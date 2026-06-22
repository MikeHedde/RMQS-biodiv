# ============================================================
# 15_gpd_IL_sensitivity_common_z.R
# Sensitivity of conditional GPD-versus-pitfall contrasts to the
# assumed GPD interception-length proxy.
#
# IMPORTANT:
# This script estimates p-hat for ALL protocols at the SAME
# standardised effort value (eff_z = 0) under every IL scenario.
# It therefore differs from the manuscript's descriptive predictions
# at method-specific median effort, which are necessarily almost
# invariant to a method-specific rescaling of GPD effort.
#
# Prerequisite: apply the small patch to 05_occupancy_fit_reviewed.R
# so fit_one_species_occu() accepts `pred_eff_z`.
# ============================================================

stopifnot(
  exists("Y"), exists("eff_mat"), exists("visit_mask"),
  exists("site_cov"), exists("fit_one_species_occu"),
  exists("gpd_avail_site"), exists("out_dir")
)

sensitivity_dir <- file.path(out_dir, "reviewer_response", "gpd_IL_sensitivity_common_z")
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

# Broad, plausible range around the geometric baseline:
# 1.16 m = 1 m drift element + pi * 0.05 m pitfall circumference.
gpd_il_grid_m <- get0(
  "gpd_il_sensitivity_m",
  ifnotfound = c(0.80, 1.00, 1.16, 1.30, 1.50)
)
gpd_il_grid_m <- sort(unique(as.numeric(gpd_il_grid_m)))
stopifnot(all(is.finite(gpd_il_grid_m)), all(gpd_il_grid_m > 0))

methods_compare <- c("Pitfall2", "Pitfall4", "Pitfall6", "Pitfall8", "Pitfall10", "GPD")
methods_compare <- methods_compare[methods_compare %in% colnames(eff_mat)]
stopifnot("GPD" %in% methods_compare)

# The common prediction point: pooled mean of log-transformed effort
# after standardisation, i.e. eff_z = 0, for every method and scenario.
# This estimates conditional protocol contrasts at the same standardised
# effort, rather than observed p-hat at each method's own median effort.
common_eff_z <- 0

# Extract site-level GPD components produced by 04_detection_effort_reviewed.R
gpd_components <- tibble::tibble(site_id = rownames(eff_mat)) %>%
  dplyr::left_join(
    gpd_avail_site %>%
      dplyr::select(site_id, n_gpd_eff, days_gpd, visit_available),
    by = "site_id"
  ) %>%
  dplyr::mutate(
    visit_available = dplyr::coalesce(as.logical(visit_available), FALSE),
    n_gpd_eff = dplyr::coalesce(as.numeric(n_gpd_eff), 0),
    days_gpd = as.numeric(days_gpd)
  )

if (!identical(
  as.logical(gpd_components$visit_available),
  as.logical(visit_mask[rownames(eff_mat), "GPD"])
)) {
  stop("Mismatch between GPD availability in gpd_avail_site and visit_mask.")
}

if (any(gpd_components$visit_available &
        (!is.finite(gpd_components$days_gpd) | gpd_components$days_gpd <= 0))) {
  stop("At least one available GPD visit has invalid trapping duration.")
}
if (any(gpd_components$visit_available & gpd_components$n_gpd_eff <= 0)) {
  stop("At least one available GPD visit has zero usable sub-units.")
}

standardise_effort_matrix <- function(effort_matrix, availability_mask) {
  out <- effort_matrix
  out[!availability_mask] <- NA_real_
  log_eff <- log1p(out)
  centre <- mean(log_eff, na.rm = TRUE)
  scale_ <- stats::sd(log_eff, na.rm = TRUE)
  if (!is.finite(scale_) || is.na(scale_) || scale_ == 0) scale_ <- 1
  
  z <- (log_eff - centre) / scale_
  z[!availability_mask] <- NA_real_
  z
}

summarise_effort_support <- function(eff_z, il_value) {
  purrr::map_dfr(methods_compare, function(m) {
    z <- eff_z[, m]
    z <- z[is.finite(z)]
    tibble::tibble(
      IL_GPD_m = il_value,
      method = m,
      n_visits = length(z),
      min_eff_z = min(z),
      q25_eff_z = stats::quantile(z, 0.25, names = FALSE),
      median_eff_z = stats::median(z),
      q75_eff_z = stats::quantile(z, 0.75, names = FALSE),
      max_eff_z = max(z),
      common_reference_within_observed_range = min(z) <= common_eff_z && max(z) >= common_eff_z
    )
  })
}

fit_scenario <- function(il_per_unit_m) {
  # Replace ONLY GPD effort. IL is per active GPD sub-unit.
  eff_alt <- eff_mat
  gpd_available <- gpd_components$visit_available
  eff_alt[, "GPD"] <- NA_real_
  eff_alt[gpd_available, "GPD"] <-
    il_per_unit_m *
    gpd_components$n_gpd_eff[gpd_available] *
    gpd_components$days_gpd[gpd_available]
  
  eff_z_alt <- standardise_effort_matrix(eff_alt, visit_mask)
  
  spp_all_local <- dimnames(Y)$species
  fit_list <- stats::setNames(
    lapply(spp_all_local, function(sp) {
      tryCatch(
        fit_one_species_occu(
          sp_name = sp,
          Y = Y,
          eff_z = eff_z_alt,
          site_cov = site_cov,
          min_sites = min_sites,
          require_contrast = require_contrast,
          quiet = TRUE,
          maxit = 400,
          pred_eff_z = common_eff_z
        ),
        error = function(e) list(
          status = "fail", species = sp,
          reason = paste0("error: ", conditionMessage(e))
        )
      )
    }),
    spp_all_local
  )
  
  status_tbl <- purrr::imap_dfr(fit_list, function(x, sp) {
    tibble::tibble(
      IL_GPD_m = il_per_unit_m,
      species = sp,
      status = x$status %||% "fail",
      reason = x$reason %||% NA_character_,
      n_sites_available = x$n_sites_available %||% NA_integer_,
      n_sites_detected = x$n_sites_detected %||% NA_integer_,
      n_rep = x$n_rep %||% NA_integer_,
      fit_quality = x$fit_quality %||% NA_character_
    )
  })
  
  p_tbl <- purrr::map_dfr(fit_list, function(x) {
    if (is.null(x$p_table)) return(tibble::tibble())
    x$p_table
  }) %>%
    dplyr::mutate(
      IL_GPD_m = il_per_unit_m,
      prediction_eff_z = common_eff_z,
      .before = 1
    )
  
  list(
    p = p_tbl,
    status = status_tbl,
    support = summarise_effort_support(eff_z_alt, il_per_unit_m)
  )
}

message("Running common-effort GPD IL sensitivity refits: ",
        paste(gpd_il_grid_m, collapse = ", "), " m")
scenario_fits <- lapply(gpd_il_grid_m, fit_scenario)
names(scenario_fits) <- sprintf("IL_%0.2f", gpd_il_grid_m)

p_all <- dplyr::bind_rows(lapply(scenario_fits, `[[`, "p"))
status_all <- dplyr::bind_rows(lapply(scenario_fits, `[[`, "status"))
support_all <- dplyr::bind_rows(lapply(scenario_fits, `[[`, "support"))

readr::write_csv(status_all, file.path(sensitivity_dir, "gpd_IL_common_z_model_status.csv"))
readr::write_csv(support_all, file.path(sensitivity_dir, "gpd_IL_common_z_effort_support.csv"))

# Keep the same species in every IL scenario and for every protocol.
valid_species <- p_all %>%
  dplyr::filter(method %in% methods_compare) %>%
  dplyr::group_by(IL_GPD_m, species) %>%
  dplyr::summarise(
    all_protocols_predicted = dplyr::n_distinct(method) == length(methods_compare),
    .groups = "drop"
  ) %>%
  dplyr::filter(all_protocols_predicted) %>%
  dplyr::count(species, name = "n_IL") %>%
  dplyr::filter(n_IL == length(gpd_il_grid_m)) %>%
  dplyr::pull(species)

if (length(valid_species) < 5) {
  stop("Fewer than five common species have valid p-hat estimates across all IL scenarios.")
}

p_common <- p_all %>%
  dplyr::filter(species %in% valid_species, method %in% methods_compare) %>%
  dplyr::mutate(
    method = factor(method, levels = methods_compare),
    p_clip = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    logit_p = stats::qlogis(p_clip),
    se_p = pmax((ucl - lcl) / (2 * 1.96), 1e-6),
    se_logit = pmax(se_p / (p_clip * (1 - p_clip)), 1e-6),
    weight = pmin(1 / se_logit^2, 1e6)
  )

readr::write_csv(tibble::tibble(species = valid_species),
                 file.path(sensitivity_dir, "gpd_IL_common_z_species.csv"))
readr::write_csv(p_common,
                 file.path(sensitivity_dir, "gpd_IL_common_z_p_hat_by_species.csv"))

p_summary <- p_common %>%
  dplyr::group_by(IL_GPD_m, method) %>%
  dplyr::summarise(
    n_species = dplyr::n_distinct(species),
    median_p = stats::median(p_hat, na.rm = TRUE),
    q25_p = stats::quantile(p_hat, 0.25, na.rm = TRUE),
    q75_p = stats::quantile(p_hat, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

fit_contrast_model <- function(dat) {
  mod <- lme4::lmer(
    logit_p ~ method + (1 | species),
    data = dat,
    weights = weight,
    REML = FALSE
  )
  
  em <- emmeans::emmeans(mod, ~ method)
  out <- emmeans::contrast(
    em,
    method = list(
      "GPD - Pitfall4" = c(0, -1, 0, 0, 0, 1),
      "GPD - Pitfall6" = c(0, 0, -1, 0, 0, 1),
      "GPD - Pitfall8" = c(0, 0, 0, -1, 0, 1)
    ),
    adjust = "none"
  )
  
  as.data.frame(summary(out, infer = TRUE)) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      contrast,
      estimate_logit = estimate,
      lcl_logit = lower.CL,
      ucl_logit = upper.CL,
      odds_ratio = exp(estimate),
      lcl_OR = exp(lower.CL),
      ucl_OR = exp(upper.CL),
      p_value = p.value
    )
}

contrast_tbl <- p_common %>%
  dplyr::group_by(IL_GPD_m) %>%
  dplyr::group_modify(~ fit_contrast_model(.x)) %>%
  dplyr::ungroup()

pivot_p <- p_summary %>%
  dplyr::select(IL_GPD_m, method, median_p) %>%
  tidyr::pivot_wider(names_from = method, values_from = median_p)

position_tbl <- pivot_p %>%
  dplyr::mutate(
    GPD_position_at_common_eff_z = dplyr::case_when(
      GPD < Pitfall4 ~ "below Pitfall4",
      GPD <= Pitfall6 ~ "between Pitfall4 and Pitfall6",
      GPD <= Pitfall8 ~ "between Pitfall6 and Pitfall8",
      TRUE ~ "above Pitfall8"
    )
  )

readr::write_csv(p_summary, file.path(sensitivity_dir, "gpd_IL_common_z_protocol_summary.csv"))
readr::write_csv(position_tbl, file.path(sensitivity_dir, "gpd_IL_common_z_position.csv"))
readr::write_csv(contrast_tbl, file.path(sensitivity_dir, "gpd_IL_common_z_contrasts.csv"))

p_protocol <- ggplot2::ggplot(
  p_summary,
  ggplot2::aes(x = IL_GPD_m, y = median_p, group = method, shape = method, linetype = method)
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = q25_p, ymax = q75_p),
    width = 0.025,
    alpha = 0.65
  ) +
  ggplot2::geom_vline(xintercept = 1.16, linetype = "dashed") +
  ggplot2::labs(
    x = "Assumed GPD interception-length proxy (m per active sub-unit)",
    y = expression("Median predicted detection probability at common " * z[effort] * " = 0"),
    shape = "Protocol",
    linetype = "Protocol"
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  file.path(sensitivity_dir, "FigS_GPD_IL_common_z_protocol_p.png"),
  p_protocol, width = 8.5, height = 5.5, dpi = 400
)

p_contrast <- ggplot2::ggplot(
  contrast_tbl,
  ggplot2::aes(x = IL_GPD_m, y = odds_ratio, group = contrast, shape = contrast, linetype = contrast)
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lcl_OR, ymax = ucl_OR),
    width = 0.025,
    alpha = 0.65
  ) +
  ggplot2::geom_vline(xintercept = 1.16, linetype = "dashed") +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = "Assumed GPD interception-length proxy (m per active sub-unit)",
    y = "GPD / pitfall detection odds ratio at common effort",
    shape = "Contrast",
    linetype = "Contrast"
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  file.path(sensitivity_dir, "FigS_GPD_IL_common_z_contrasts.png"),
  p_contrast, width = 8.5, height = 5.5, dpi = 400
)

message("\nCommon-effort GPD IL sensitivity completed.")
message("Common species across scenarios: ", length(valid_species))
message("Inspect gpd_IL_common_z_effort_support.csv before interpreting predictions.")
print(position_tbl)
print(contrast_tbl)
