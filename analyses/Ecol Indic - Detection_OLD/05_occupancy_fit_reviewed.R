# ============================================================
# 05_occupancy_fit_reviewed.R
# Species-by-species single-season occupancy models
# Reviewer-ready version for Pitfall + GPD analyses
#
# Key changes compared with the previous script:
# 1. Species inclusion is based on the number of sites with >=1 detection,
#    not on the number of sites with at least one available replicate.
# 2. Models with missing/very large SEs or extreme coefficients are classified
#    as "unstable" and are not used for p-hat prediction tables.
# 3. Empty outputs are handled safely, avoiding crashes such as
#    "object 'method' not found" when no stable model is available.
# 4. The script reports available sites, detected sites, replicates and reasons
#    for exclusion/failure for reviewer transparency.
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

tidy_unmarked_mat <- function(mat, component = c("state", "det")) {
  component <- match.arg(component)
  if (is.null(mat)) return(tibble::tibble())
  mat <- as.data.frame(mat)
  mat$term <- rownames(mat)
  rownames(mat) <- NULL
  mat %>%
    tibble::as_tibble() %>%
    dplyr::relocate(term) %>%
    dplyr::rename(
      estimate = Estimate,
      se       = SE,
      z        = z,
      p_value  = `P(>|z|)`
    ) %>%
    dplyr::mutate(component = component)
}

method_eff_medians <- function(eff_z_matrix) {
  apply(eff_z_matrix, 2, function(col) {
    if (all(is.na(col))) 0 else stats::median(col, na.rm = TRUE)
  })
}

check_unstable_fit <- function(fm,
                               se_max = get0("unstable_se_max", ifnotfound = 50),
                               coef_abs_max = get0("unstable_coef_abs_max", ifnotfound = 25)) {
  sumobj <- tryCatch(suppressWarnings(summary(fm)), error = function(e) NULL)
  
  if (is.null(sumobj)) {
    return(list(
      usable_for_prediction = FALSE,
      fit_quality = "failed_summary",
      reason = "summary_failed"
    ))
  }
  
  se_state <- suppressWarnings(as.numeric(sumobj$state[, "SE"]))
  se_det   <- suppressWarnings(as.numeric(sumobj$det[, "SE"]))
  est_state <- suppressWarnings(as.numeric(sumobj$state[, "Estimate"]))
  est_det   <- suppressWarnings(as.numeric(sumobj$det[, "Estimate"]))
  
  se_all  <- c(se_state, se_det)
  est_all <- c(est_state, est_det)
  
  bad_se <- any(!is.finite(se_all) | is.na(se_all) | se_all > se_max)
  bad_coef <- any(!is.finite(est_all) | is.na(est_all) | abs(est_all) > coef_abs_max)
  
  if (bad_se && bad_coef) {
    return(list(
      usable_for_prediction = FALSE,
      fit_quality = "failed_unstable",
      reason = "large_or_missing_se_and_extreme_coef"
    ))
  }
  
  if (bad_se) {
    return(list(
      usable_for_prediction = TRUE,
      fit_quality = "usable_large_se",
      reason = "large_or_missing_se"
    ))
  }
  
  if (bad_coef) {
    return(list(
      usable_for_prediction = TRUE,
      fit_quality = "usable_extreme_coef",
      reason = "extreme_coef"
    ))
  }
  
  list(
    usable_for_prediction = TRUE,
    fit_quality = "stable",
    reason = NA_character_
  )
}


# ------------------------------------------------------------
# Core fitting function
# ------------------------------------------------------------

fit_one_species_occu <- function(sp_name,
                                 Y,
                                 eff_z,
                                 site_cov,
                                 form_state_terms = c("DOY_z"),
                                 min_sites = min_sites,
                                 require_contrast = TRUE,
                                 quiet = TRUE,
                                 maxit = 400) {

  k <- which(dimnames(Y)$species == sp_name)
  if (length(k) != 1) {
    return(list(status = "skip", species = sp_name, reason = "species_not_found"))
  }

  Yk0 <- Y[, , k, drop = FALSE][, , 1]
  if (is.null(rownames(Yk0)) || is.null(colnames(Yk0))) {
    return(list(status = "fail", species = sp_name, reason = "Yk_missing_dimnames"))
  }

  # Available sites = at least one replicate/protocol with usable data
  keep_site <- apply(Yk0, 1, function(a) any(!is.na(a)))
  Yk <- Yk0[keep_site, , drop = FALSE]

  n_sites_available <- nrow(Yk)
  n_sites_detected <- sum(rowSums(Yk == 1, na.rm = TRUE) > 0)

  # Species inclusion should be based on sites where the species was detected,
  # not on all sites where the protocol was available.
  if (n_sites_detected < min_sites) {
    return(list(
      status = "skip", species = sp_name,
      n_sites_available = n_sites_available,
      n_sites_detected = n_sites_detected,
      reason = paste0("n_detected_sites<", min_sites)
    ))
  }

  # Keep only reps with at least one available observation across retained sites
  keep_rep <- apply(Yk, 2, function(a) any(!is.na(a)))
  Yk <- Yk[, keep_rep, drop = FALSE]
  rep_names_k <- colnames(Yk)
  n_rep <- ncol(Yk)

  if (n_rep < 2) {
    return(list(status = "skip", species = sp_name,
                n_sites_available = n_sites_available,
                n_sites_detected = n_sites_detected,
                n_rep = n_rep,
                reason = "n_rep<2"))
  }

  if (require_contrast) {
    yvec <- as.vector(Yk[!is.na(Yk)])
    if (!(any(yvec == 1) && any(yvec == 0))) {
      return(list(status = "skip", species = sp_name,
                  n_sites_available = n_sites_available,
                  n_sites_detected = n_sites_detected,
                  n_rep = n_rep,
                  reason = "no_detection_contrast"))
    }
  }

  # ----------------------------------------------------------
  # Observation covariates: method indicators + effort
  # ----------------------------------------------------------
  mkI <- function(tag) {
    matrix(as.numeric(rep_names_k == tag),
           nrow = nrow(Yk), ncol = ncol(Yk), byrow = TRUE,
           dimnames = dimnames(Yk))
  }
  I_list <- stats::setNames(lapply(rep_names_k, mkI), rep_names_k)

  eff_sp <- eff_z[rownames(Yk), rep_names_k, drop = FALSE]

  # obsCovs must be NA exactly where y is NA
  for (nm in rep_names_k) I_list[[nm]][is.na(Yk)] <- NA
  eff_sp[is.na(Yk)] <- NA

  # ----------------------------------------------------------
  # Site covariates
  # ----------------------------------------------------------
  sc <- site_cov %>%
    dplyr::filter(site_id %in% rownames(Yk)) %>%
    dplyr::arrange(match(site_id, rownames(Yk)))

  if (nrow(sc) != nrow(Yk)) {
    return(list(status = "fail", species = sp_name,
                n_sites_available = n_sites_available,
                n_sites_detected = n_sites_detected,
                n_rep = n_rep,
                reason = "site_cov_alignment_mismatch"))
  }

  sc_umf <- sc %>% dplyr::select(-site_id) %>% as.data.frame()
  rownames(sc_umf) <- sc$site_id

  # Keep only valid/non-constant state covariates
  available_state <- intersect(form_state_terms, colnames(sc_umf))
  if (length(available_state)) {
    v_ok <- vapply(available_state, function(v) {
      x <- sc_umf[[v]]
      all(is.finite(x) | is.na(x)) && stats::sd(x, na.rm = TRUE) > 0
    }, logical(1))
    available_state <- available_state[v_ok]
  }
  state_rhs <- if (length(available_state)) paste(available_state, collapse = " + ") else "1"

  # Detection formula. Use Pitfall10 as reference when available.
  ref_method <- if ("Pitfall10" %in% rep_names_k) "Pitfall10" else rep_names_k[1]

  eff_var_ok <- stats::sd(as.vector(eff_sp), na.rm = TRUE) > 0
  det_terms <- setdiff(rep_names_k, ref_method)
  if (eff_var_ok) det_terms <- c(det_terms, "eff_z")

  det_rhs <- if (length(det_terms)) paste(det_terms, collapse = " + ") else "1"
  form_occu <- stats::as.formula(paste("~", det_rhs, "~", state_rhs))

  obsCovs <- c(I_list, list(eff_z = eff_sp))

  umf_k <- unmarked::unmarkedFrameOccu(
    y        = Yk,
    siteCovs = sc_umf,
    obsCovs  = obsCovs
  )

  fm <- tryCatch(
    {
      if (quiet) {
        suppressWarnings(suppressMessages(
          unmarked::occu(form_occu, data = umf_k, control = list(maxit = maxit))
        ))
      } else {
        unmarked::occu(form_occu, data = umf_k, control = list(maxit = maxit))
      }
    },
    error = function(e) e
  )

  if (inherits(fm, "error")) {
    return(list(status = "fail", species = sp_name,
                n_sites_available = n_sites_available,
                n_sites_detected = n_sites_detected,
                n_rep = n_rep,
                reason = conditionMessage(fm)))
  }

  AIC_fm <- tryCatch(stats::AIC(fm), error = function(e) NA_real_)
  
  fit_check <- check_unstable_fit(fm)
  
  sumobj <- tryCatch(suppressWarnings(summary(fm)), error = function(e) NULL)
  beta_tab <- tibble::tibble()
  
  if (!is.null(sumobj)) {
    beta_tab <- dplyr::bind_rows(
      tidy_unmarked_mat(sumobj$state, component = "state"),
      tidy_unmarked_mat(sumobj$det,   component = "det")
    ) %>%
      dplyr::mutate(species = sp_name) %>%
      dplyr::rename(p = p_value) %>%
      dplyr::relocate(species, component, term, estimate, se, z, p)
  }
  
  if (!isTRUE(fit_check$usable_for_prediction)) {
    return(list(
      status = "unstable",
      species = sp_name,
      n_sites_available = n_sites_available,
      n_sites_detected = n_sites_detected,
      n_rep = n_rep,
      AIC = AIC_fm,
      reason = fit_check$reason,
      fit_quality = fit_check$fit_quality,
      beta = beta_tab,
      fit = fm
    ))
  }

  # ----------------------------------------------------------
  # Detection predictions by method at method-specific median effort
  # ----------------------------------------------------------
  eff_med_all <- method_eff_medians(eff_z)
  eff_for <- function(m) if (m %in% names(eff_med_all)) eff_med_all[[m]] else 0

  mk_pred <- function(tag) {
    nd <- data.frame(row_id = 1)

    # Detection covariates used or possibly requested by predict()
    for (nm in rep_names_k) nd[[nm]] <- as.numeric(nm == tag)
    nd$eff_z <- eff_for(tag)

    # State covariates set to mean standardized value = 0
    for (v in available_state) nd[[v]] <- 0

    nd$row_id <- NULL
    pr <- tryCatch(
      suppressWarnings(as.data.frame(unmarked::predict(fm, type = "det", newdata = nd))),
      error = function(e) e
    )
    pr
  }

  pred_list <- lapply(rep_names_k, mk_pred)
  pred_errors <- vapply(pred_list, inherits, logical(1), what = "error")
  if (any(pred_errors)) {
    return(list(status = "fail", species = sp_name,
                n_sites_available = n_sites_available,
                n_sites_detected = n_sites_detected,
                n_rep = n_rep,
                AIC = AIC_fm,
                reason = paste0("prediction_failed: ", conditionMessage(pred_list[[which(pred_errors)[1]]])),
                beta = beta_tab,
                fit = fm))
  }

  tab_p <- tibble::tibble(
    species = sp_name,
    method  = rep_names_k,
    p_hat   = vapply(pred_list, function(z) z$Predicted, numeric(1)),
    lcl     = vapply(pred_list, function(z) z$lower,     numeric(1)),
    ucl     = vapply(pred_list, function(z) z$upper,     numeric(1)),
    n_sites_available = n_sites_available,
    n_sites_detected  = n_sites_detected,
    n_rep   = n_rep,
    AIC     = AIC_fm,
    any_se_bad = FALSE
  ) %>%
    dplyr::mutate(dplyr::across(c(p_hat, lcl, ucl), ~ ifelse(is.nan(.x), NA_real_, .x)))

  # Occupancy predictions by site
  psi_pred <- tryCatch(
    suppressWarnings(as.data.frame(unmarked::predict(fm, type = "state"))),
    error = function(e) NULL
  )

  psi_tab <- tibble::tibble()
  if (!is.null(psi_pred)) {
    psi_tab <- tibble::tibble(
      site_id = rownames(Yk),
      species = sp_name,
      psi_hat = psi_pred$Predicted
    )
  }

  list(
    status     = ifelse(fit_check$fit_quality == "stable", "ok", "usable_unstable"),
    species    = sp_name,
    n_sites_available = n_sites_available,
    n_sites_detected  = n_sites_detected,
    n_rep      = n_rep,
    AIC        = AIC_fm,
    reason     = fit_check$reason,
    fit_quality = fit_check$fit_quality,
    any_se_bad = fit_check$fit_quality != "stable",
    p_table    = tab_p,
    beta       = beta_tab,
    psi_hat    = psi_tab,
    fit        = fm
  )
}

# ------------------------------------------------------------
# Run all species
# ------------------------------------------------------------

spp_all <- dimnames(Y)$species

safe_fit <- purrr::safely(
  ~ fit_one_species_occu(.x, Y = Y, eff_z = eff_z, site_cov = site_cov,
                         min_sites = min_sites,
                         require_contrast = TRUE, quiet = TRUE, maxit = 400),
  otherwise = list(status = "fail", species = NA_character_, reason = "safely_otherwise")
)

res <- tibble::tibble(species = spp_all) %>%
  dplyr::mutate(
    out    = purrr::map(species, safe_fit),
    result = purrr::map(out, "result"),
    error  = purrr::map(out, "error"),
    status = purrr::map_chr(result, ~ if (is.null(.x)) "fail" else .x$status),
    reason = purrr::map_chr(out, ~ {
      if (!is.null(.x$result)) return(.x$result$reason %||% NA_character_)
      if (!is.null(.x$error))  return(conditionMessage(.x$error))
      NA_character_
    }),
    n_sites_available = purrr::map_int(result, ~ .x$n_sites_available %||% NA_integer_),
    n_sites_detected  = purrr::map_int(result, ~ .x$n_sites_detected  %||% NA_integer_),
    n_rep             = purrr::map_int(result, ~ .x$n_rep %||% NA_integer_)
  )

cat("\nOccupancy fit status summary:\n")
print(res %>% dplyr::count(status, reason, sort = TRUE), n = Inf)

readr::write_csv(res %>% dplyr::select(species, status, reason, n_sites_available, n_sites_detected, n_rep),
                 file.path(out_dir, "nb_esp_modelized.csv"))

# ------------------------------------------------------------
# Species selection / representativeness summary
# ------------------------------------------------------------

rank_freq <- dat0 %>%
  dplyr::select(site_id, species, abund) %>%
  dplyr::group_by(site_id, species) %>%
  dplyr::summarise(sumAb = sum(abund, na.rm = TRUE), .groups = "drop")

sp_freq <- rank_freq %>%
  dplyr::mutate(present = sumAb > 0) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_sites = sum(present, na.rm = TRUE),
    total_abundance = sum(sumAb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    res %>% dplyr::select(species, status, reason),
    by = "species"
  ) %>%
  dplyr::mutate(
    occupancy_fit_status = dplyr::case_when(
      status == "ok"              ~ "stable",
      status == "usable_unstable" ~ "usable_unstable",
      status == "unstable"        ~ "unstable",
      status == "skip"            ~ "skipped",
      TRUE                        ~ "not_modelled"
    ),
    final_inclusion = dplyr::case_when(
      status %in% c("ok", "usable_unstable") ~ "Included",
      TRUE                                   ~ "Excluded"
    )
  ) %>%
  dplyr::arrange(dplyr::desc(n_sites)) %>%
  dplyr::mutate(rank = dplyr::row_number())

readr::write_csv(
  sp_freq,
  file.path(out_dir, "species_selection/species_frequency_model_status.csv")
)

sp_freq_summary <- sp_freq %>%
  dplyr::count(occupancy_fit_status, final_inclusion, name = "n_species")

readr::write_csv(
  sp_freq_summary,
  file.path(out_dir, "species_selection/species_frequency_model_status_summary.csv")
)
print(sp_freq_summary)

# Figures only if ggplot2 is available and directories exist
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_rank_freq <- ggplot2::ggplot(sp_freq, ggplot2::aes(x = rank, y = n_sites, colour = model_status)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::geom_hline(yintercept = min_sites, linetype = "dashed") +
    ggplot2::scale_y_continuous(trans = "log10", breaks = c(1, 2, 5, 10, 20, 50), name = "Number of sites occupied (log scale)") +
    ggplot2::labs(x = "Species rank (by decreasing site frequency)", colour = "Species status") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "species_selection/rang_freq_selection.png"), p_rank_freq, width = 12, height = 12, dpi = 600)

  p_hist_sites <- ggplot2::ggplot(sp_freq, ggplot2::aes(x = n_sites, fill = model_status)) +
    ggplot2::geom_histogram(binwidth = 1, position = "identity", alpha = 0.6) +
    ggplot2::geom_vline(xintercept = min_sites, linetype = "dashed") +
    ggplot2::labs(x = "Number of sites occupied", y = "Number of species", fill = "Species status") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "species_selection/hist_sites_selection.png"), p_hist_sites, width = 12, height = 12, dpi = 600)
}

# ------------------------------------------------------------
# Extract outputs safely
# ------------------------------------------------------------

res_ok <- res %>%
  dplyr::filter(status %in% c("ok", "usable_unstable"))

res_unstable <- res %>%
  dplyr::filter(status == "unstable")

if (nrow(res_ok) > 0) {
  tab_p <- res_ok %>%
    dplyr::transmute(tbl = purrr::map(result, "p_table")) %>%
    tidyr::unnest(tbl)
  tab_beta <- res_ok %>%
    dplyr::transmute(tbl = purrr::map(result, "beta")) %>%
    tidyr::unnest(tbl)
  tab_psi <- res_ok %>%
    dplyr::transmute(tbl = purrr::map(result, "psi_hat")) %>%
    tidyr::unnest(tbl)
} else {
  tab_p <- tibble::tibble(species = character(), method = character(), p_hat = numeric(), lcl = numeric(), ucl = numeric(),
                          n_sites_available = integer(), n_sites_detected = integer(), n_rep = integer(), AIC = numeric(), any_se_bad = logical())
  tab_beta <- tibble::tibble(species = character(), component = character(), term = character(), estimate = numeric(), se = numeric(), z = numeric(), p = numeric())
  tab_psi <- tibble::tibble(site_id = character(), species = character(), psi_hat = numeric())
}

# Store beta for unstable models separately, useful for diagnostics but not for inference
if (nrow(res_unstable) > 0) {
  tab_beta_unstable <- res_unstable %>%
    dplyr::transmute(tbl = purrr::map(result, "beta")) %>%
    tidyr::unnest(tbl)
} else {
  tab_beta_unstable <- tibble::tibble()
}

readr::write_csv(tab_p,    file.path(out_dir, "occ_model_output/p_hat_by_method_unmarked_full.csv"))
readr::write_csv(tab_beta, file.path(out_dir, "occ_model_output/beta_by_method_unmarked_full.csv"))
readr::write_csv(tab_psi,  file.path(out_dir, "occ_model_output/psi_hat_by_method_unmarked_full.csv"))
readr::write_csv(tab_beta_unstable, file.path(out_dir, "occ_model_output/beta_by_method_unmarked_unstable_diagnostics.csv"))

cat("\nStable occupancy models:", nrow(res_ok), "species\n")
cat("Unstable fitted models excluded from p-hat tables:", nrow(res_unstable), "species\n")

# ------------------------------------------------------------
# Optional gt annexes, safely skipped if no stable model is available
# ------------------------------------------------------------

if (requireNamespace("gt", quietly = TRUE) && nrow(res_ok) > 0) {
  tab_status <- res_ok %>%
    dplyr::transmute(
      species,
      n_sites_available = purrr::map_int(result, "n_sites_available"),
      n_sites_detected  = purrr::map_int(result, "n_sites_detected"),
      n_rep   = purrr::map_int(result, "n_rep"),
      AIC     = purrr::map_dbl(result, ~ .x$AIC %||% NA_real_),
      any_se_bad = purrr::map_lgl(result, ~ .x$any_se_bad %||% FALSE)
    ) %>%
    dplyr::arrange(n_sites_detected)

  gt_status <- tab_status %>%
    gt::gt() %>%
    gt::tab_header(title = "Annex — Occupancy models (unmarked::occu): diagnostic per species") %>%
    gt::fmt_number(columns = c(AIC), decimals = 2) %>%
    gt::cols_label(
      species = "Species",
      n_sites_available = "Available sites",
      n_sites_detected = "Detected sites",
      n_rep = "Replicates",
      any_se_bad = "Unstable SE?",
      AIC = "AIC"
    )
  gt::gtsave(gt_status, filename = file.path(out_dir, "occ_model_output/Annex_occu_status.html"))

  gt_p <- tab_p %>%
    dplyr::mutate(method = factor(method, levels = unique(method))) %>%
    gt::gt(groupname_col = "species") %>%
    gt::tab_header(
      title = "Annex — Detection probabilities by method",
      subtitle = "Predictions at method-specific median effort, 95% CI"
    ) %>%
    gt::fmt_number(columns = c(p_hat, lcl, ucl), decimals = 3) %>%
    gt::cols_label(method = "Method", p_hat = gt::html("p&#770;"), lcl = "LCL", ucl = "UCL") %>%
    gt::tab_options(row_group.as_column = TRUE)
  gt::gtsave(gt_p, filename = file.path(out_dir, "occ_model_output/Annex_p_hat_by_method.html"))

  gt_beta <- tab_beta %>%
    gt::gt(groupname_col = "species") %>%
    gt::tab_header(
      title = "Annex — Model coefficients (logit scale)",
      subtitle = "State = occupancy, Det = detection"
    ) %>%
    gt::fmt_number(columns = c(estimate, se, z, p), decimals = 3) %>%
    gt::tab_options(row_group.as_column = TRUE)
  gt::gtsave(gt_beta, filename = file.path(out_dir, "occ_model_output/Annex_occu_coefficients.html"))
} else {
  cat("\nNo gt annex generated: either gt is not installed or no stable occupancy model is available.\n")
}


# ------------------------------------------------------------
# Missing data summary in the detection matrix
# ------------------------------------------------------------

missing_global <- tibble::tibble(
  n_cells_total   = length(Y),
  n_cells_missing = sum(is.na(Y)),
  prop_missing    = mean(is.na(Y)),
  perc_missing    = 100 * mean(is.na(Y))
)

print(missing_global)

readr::write_csv(
  missing_global,
  file.path(out_dir, "occ_model_output/missing_data_global_summary.csv")
)

missing_by_protocol <- tibble::tibble(
  protocol = dimnames(Y)[[2]],
  n_cells_total = apply(Y, 2, length),
  n_cells_missing = apply(Y, 2, function(x) sum(is.na(x))),
  prop_missing = apply(Y, 2, function(x) mean(is.na(x))),
  perc_missing = 100 * prop_missing
)

print(missing_by_protocol)

readr::write_csv(
  missing_by_protocol,
  file.path(out_dir, "occ_model_output/missing_data_by_protocol.csv")
)

availability_site_protocol <- apply(Y, c(1, 2), function(x) any(!is.na(x)))

missing_site_protocol <- tibble::tibble(
  n_site_protocol_total = length(availability_site_protocol),
  n_site_protocol_missing = sum(!availability_site_protocol),
  prop_missing = mean(!availability_site_protocol),
  perc_missing = 100 * mean(!availability_site_protocol)
)

print(missing_site_protocol)

readr::write_csv(
  missing_site_protocol,
  file.path(out_dir, "occ_model_output/missing_site_protocol_summary.csv")
)

missing_site_protocol_by_protocol <- tibble::tibble(
  protocol = dimnames(Y)[[2]],
  n_sites_total = nrow(availability_site_protocol),
  n_sites_missing = colSums(!availability_site_protocol),
  prop_missing = n_sites_missing / n_sites_total,
  perc_missing = 100 * prop_missing
)

print(missing_site_protocol_by_protocol)

readr::write_csv(
  missing_site_protocol_by_protocol,
  file.path(out_dir, "occ_model_output/missing_site_protocol_by_protocol.csv")
)



tab_p %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    median_p = median(p_hat, na.rm = TRUE),
    q25 = quantile(p_hat, 0.25, na.rm = TRUE),
    q75 = quantile(p_hat, 0.75, na.rm = TRUE),
    min_p = min(p_hat, na.rm = TRUE),
    max_p = max(p_hat, na.rm = TRUE),
    .groups = "drop"
  )