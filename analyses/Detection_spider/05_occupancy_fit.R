# ============================================================
# MODÈLES unmarked::occu ESPÈCE-PAR-ESPÈCE
# Sorties fichiers: tab_p, tab_beta, tab_psi + CSV
# ============================================================


#--- utilitaire: extraire matrice "Estimate/SE/z/P"
tidy_unmarked_mat <- function(mat, component = c("state", "det")) {
  component <- match.arg(component)
  if (is.null(mat)) return(tibble())
  mat <- as.data.frame(mat)
  mat$term <- rownames(mat)
  rownames(mat) <- NULL
  mat %>%
    as_tibble() %>%
    relocate(term) %>%
    rename(
      estimate = Estimate,
      se       = SE,
      z        = z,
      p_value  = `P(>|z|)`
    ) %>%
    mutate(component = component)
}

#--- diagnostic simple
flag_bad_fit <- function(fm) {
  co <- summary(fm)
  mats <- list(co$state, co$det)
  any_se_bad <- any(map_lgl(mats, ~ any(is.na(.x[,"SE"]) | is.nan(.x[,"SE"]) | is.infinite(.x[,"SE"]))))
  list(any_se_bad = any_se_bad)
}

method_eff_medians <- function(eff_z_matrix) {
  apply(eff_z_matrix, 2, function(col) if (all(is.na(col))) 0 else stats::median(col, na.rm = TRUE))
}
  
fit_one_species_occu <- function(sp_name,
                                 Y,
                                 eff_z,
                                 site_cov,
                                 form_state_terms = c("DOY_z"), # "ALTITUDE_z", "NDVI_m_z", "NDVI_sd_z"),
                                 max_state_terms = 1,
                                 min_sites = min_sites,
                                 require_contrast = TRUE,
                                 quiet = TRUE,
                                 maxit = 200) {
  
  #--------------------------------------------------
  # 0) Espèce
  #--------------------------------------------------
  k <- which(dimnames(Y)$species == sp_name)
  if (length(k) != 1) {
    return(list(status = "skip", species = sp_name, reason = "species_not_found"))
  }
  
  Yk0 <- Y[, , k, drop = FALSE][, , 1] # sites x reps
  if (is.null(rownames(Yk0)) || is.null(colnames(Yk0))) {
    return(list(status = "fail", species = sp_name, reason = "Yk missing row/col names"))
  }
  
  #--------------------------------------------------
  # 1) Filtrer sites / reps informatifs
  #--------------------------------------------------
  keep_site <- apply(Yk0, 1, function(a) any(!is.na(a)))
  Yk <- Yk0[keep_site, , drop = FALSE]
  n_sites <- nrow(Yk)
  
  if (n_sites < min_sites) {
    return(list(status = "skip", species = sp_name, n_sites = n_sites,
                reason = paste0("n_sites<", min_sites)))
  }
  
  keep_rep <- apply(Yk, 2, function(a) any(!is.na(a)))
  Yk <- Yk[, keep_rep, drop = FALSE]
  n_rep <- ncol(Yk)
  
  if (n_rep < 2) {
    return(list(status = "skip", species = sp_name, n_sites = n_sites, n_rep = n_rep,
                reason = "n_rep<2"))
  }
  
  if (require_contrast) {
    yvec <- as.vector(Yk[!is.na(Yk)])
    if (!(any(yvec == 1) && any(yvec == 0))) {
      return(list(status = "skip", species = sp_name, n_sites = n_sites, n_rep = n_rep,
                  reason = "no_det_contrast"))
    }
  }
  
  rep_names_k <- colnames(Yk)
  
  #--------------------------------------------------
  # 2) obsCovs : indicateurs de méthode + effort
  #--------------------------------------------------
  mkI <- function(tag) {
    matrix(as.numeric(rep_names_k == tag),
           nrow = nrow(Yk), ncol = ncol(Yk), byrow = TRUE,
           dimnames = dimnames(Yk))
  }
  I_list <- setNames(lapply(rep_names_k, mkI), rep_names_k)
  
  # effort aligné
  eff_sp <- eff_z[rownames(Yk), rep_names_k, drop = FALSE]
  
  # masquer là où Y est NA
  for (nm in rep_names_k) I_list[[nm]][is.na(Yk)] <- NA
  eff_sp[is.na(Yk)] <- NA
  
  #--------------------------------------------------
  # 3) siteCovs (alignement strict sur Yk)
  #--------------------------------------------------
  sc <- site_cov %>%
    dplyr::filter(site_id %in% rownames(Yk)) %>%
    dplyr::arrange(match(site_id, rownames(Yk)))
  
  if (nrow(sc) != nrow(Yk)) {
    return(list(status = "fail", species = sp_name, n_sites = n_sites, n_rep = n_rep,
                reason = "site_cov alignment mismatch"))
  }
  
  # IMPORTANT : as.data.frame() pour éviter les warnings "rownames on tibble"
  sc_umf <- sc %>% dplyr::select(-site_id) %>% as.data.frame()
  rownames(sc_umf) <- sc$site_id
  
  #--------------------------------------------------
  # 4) Formules
  #--------------------------------------------------
  available_state <- intersect(form_state_terms, colnames(sc_umf))
  if (length(available_state)) {
    v_ok <- vapply(available_state, function(v) stats::sd(sc_umf[[v]], na.rm = TRUE) > 0, logical(1))
    available_state <- available_state[v_ok]
  }
  state_rhs <- if (length(available_state)) paste(available_state, collapse = " + ") else "1"
 
  det_terms <- c(rep_names_k, "eff_z")
  
  # drop 1 méthode pour éviter la colinéarité avec l'intercept
  if (length(rep_names_k) >= 2) {
    ref_method <- rep_names_k[1]              # ou fixe: "Pitfall" etc.
    det_terms  <- setdiff(det_terms, ref_method)
  }
  det_rhs <- paste(det_terms, collapse = " + ")
  form_occu <- as.formula(paste("~", det_rhs, "~", state_rhs))
  
  #--------------------------------------------------
  # 5) Construire umf_k et fitter
  #--------------------------------------------------
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
    return(list(status = "fail", species = sp_name, n_sites = n_sites, n_rep = n_rep,
                reason = conditionMessage(fm)))
  }
  
  # AIC (APRÈS fit)
  AIC_fm <- tryCatch(stats::AIC(fm), error = function(e) NA_real_)
  
  # Flag SE "bizarres" (APRÈS fit)
  sumobj <- tryCatch(summary(fm), error = function(e) NULL)
  any_se_bad <- NA
  if (!is.null(sumobj)) {
    se_state <- suppressWarnings(as.numeric(sumobj$state[, "SE"]))
    se_det   <- suppressWarnings(as.numeric(sumobj$det[, "SE"]))
    se_all   <- c(se_state, se_det)
    any_se_bad <- any(!is.finite(se_all) | is.na(se_all) | se_all > 50) # seuil ajustable
  }
  
  #--------------------------------------------------
  # 6) p̂ par méthode (effort médian global)
  #--------------------------------------------------
  eff_med_all <- method_eff_medians(eff_z)
  eff_for <- function(m) if (m %in% names(eff_med_all)) eff_med_all[[m]] else 0
  
  mk_pred <- function(tag) {
    nd <- data.frame(eff_z = eff_for(tag))
    
    # indicateurs méthode
    for (nm in rep_names_k) nd[[nm]] <- as.numeric(nm == tag)
    
    # ces colonnes n'appartiennent PAS à la formule de détection, mais predict()
    # peut réclamer des colonnes -> on les met à 0 par sécurité
    for (v in available_state) nd[[v]] <- 0
    
    as.data.frame(unmarked::predict(fm, type = "det", newdata = nd))
  }
  
  pred_list <- lapply(rep_names_k, mk_pred)
  
  tab_p <- tibble::tibble(
    species = sp_name,
    method  = rep_names_k,
    p_hat   = vapply(pred_list, function(z) z$Predicted, numeric(1)),
    lcl     = vapply(pred_list, function(z) z$lower,     numeric(1)),
    ucl     = vapply(pred_list, function(z) z$upper,     numeric(1)),
    n_sites = n_sites,
    n_rep   = n_rep,
    AIC     = AIC_fm,
    any_se_bad = any_se_bad
  ) %>%
    dplyr::mutate(
      dplyr::across(c(lcl, ucl), ~ ifelse(is.nan(.x), NA_real_, .x))
    )
  
  #--------------------------------------------------
  # 7) Coefficients (state + det)
  #--------------------------------------------------
  sumobj <- summary(fm)
  
  beta_state <- tidy_unmarked_mat(sumobj$state, component = "state")
  beta_det   <- tidy_unmarked_mat(sumobj$det,   component = "det")
  
  beta_tab <- bind_rows(beta_state, beta_det) %>%
    mutate(species = sp_name) %>%
    rename(p = p_value) %>%
    relocate(species, component, term, estimate, se, z, p)
  
  
  #--------------------------------------------------
  # 8) psi_hat par site
  #--------------------------------------------------
  psi_pred <- as.data.frame(unmarked::predict(fm, type = "state"))
  psi_tab <- tibble::tibble(
    site_id = rownames(Yk),
    species = sp_name,
    psi_hat = psi_pred$Predicted
  )
  
  #--------------------------------------------------
  # 9) Sortie
  #--------------------------------------------------
  list(
    status     = "ok",
    species    = sp_name,
    n_sites    = n_sites,
    n_rep      = n_rep,
    AIC        = AIC_fm,
    any_se_bad = any_se_bad,
    p_table    = tab_p,
    beta       = beta_tab,
    psi_hat    = psi_tab,
    fit        = fm
  )
}

  spp_all <- dimnames(Y)$species
  
  safe_fit <- purrr::safely(
    ~ fit_one_species_occu(.x, Y = Y, eff_z = eff_z, site_cov = site_cov,
                           max_state_terms = 1, min_sites = min_sites,
                           require_contrast = TRUE, quiet = TRUE, maxit = 400), 
    otherwise = list(status="fail", species=NA_character_, reason="safely_otherwise")
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
      })
    )
  
  # Combien d'espèces modélisées
  res %>% dplyr::count(status, reason, sort = TRUE)
  readr::write_csv(res, file.path(out_dir, "nb_esp_modelized.csv"))
  
  # diagramme rang-fréquence
  rank_freq <- dat0 %>%
    select(site_id, species, abund) %>%
    group_by(site_id, species) %>%
    summarise(sumAb = sum(abund), .groups = "drop")
  
  sp_freq <- rank_freq %>%
    mutate(present = sumAb > 0) %>%
    group_by(species) %>%
    summarise(
      n_sites = sum(present, na.rm = TRUE),
      total_abundance = sum(sumAb, na.rm = TRUE),
      .groups = "drop"
    )
  
  status_df <- res %>%
    select(species, status)
  
  sp_freq <- sp_freq %>%
    left_join(status_df, by = "species") %>%
    mutate(
      model_status = case_when(
        status == "ok"   ~ "Modelled",
        TRUE             ~ "Excluded"
      )
    )
  
  sp_freq <- sp_freq %>%
    arrange(desc(n_sites)) %>%
    mutate(rank = row_number())

  p_rank_freq <- ggplot(
    sp_freq,
    aes(x = rank, y = n_sites, colour = model_status)
  ) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = min_sites, linetype = "dashed") +
    scale_y_continuous(
      trans = "log10",
      breaks = c(1, 2, 5, 10, 20, 50),
      name = "Number of sites occupied (log scale)"
    ) +
    scale_colour_manual(
      values = c(
        "Modelled" = "#1b9e77",
        "Excluded" = "grey70"
      )
    ) +
    labs(
      x = "Species rank (by decreasing site frequency)",
      colour = "Species status"
    ) +
    theme_minimal()
  
    ggsave(file.path(out_dir, "species_selection/rang_freq_selection.png"), p_rank_freq, width = 12, height = 12, dpi = 600)
         
  p_hist_sites <- ggplot(
    sp_freq,
    aes(x = n_sites, fill = model_status)
  ) +
    geom_histogram(binwidth = 1, position = "identity", alpha = 0.6) +
    geom_vline(xintercept = min_sites, linetype = "dashed") +
    scale_fill_manual(
      values = c(
        "Modelled" = "#1b9e77",
        "Excluded" = "grey70"
      )
    ) +
    labs(
      x = "Number of sites occupied",
      y = "Number of species",
      fill = "Species status"
    ) +
    theme_minimal()
  ggsave(file.path(out_dir, "species_selection/hist_sites_selection.png"), p_hist_sites, width = 12, height = 12, dpi = 600)
  
  sp_freq %>%
    summarise(
      n_species_total = n(),
      n_modelled      = sum(model_status == "Modelled"),
      n_excluded      = sum(model_status == "Excluded"),
      prop_modelled   = mean(model_status == "Modelled")
    )
  
  
  res_ok <- res %>% dplyr::filter(status == "ok")
  
  tab_p    <- res_ok %>% dplyr::transmute(tbl = purrr::map(result, "p_table")) %>% tidyr::unnest(tbl)
  tab_beta <- res_ok %>% dplyr::transmute(tbl = purrr::map(result, "beta"))    %>% tidyr::unnest(tbl)
  tab_psi  <- res_ok %>% dplyr::transmute(tbl = purrr::map(result, "psi_hat")) %>% tidyr::unnest(tbl)
  
  readr::write_csv(tab_p,    file.path(out_dir, "occ_model_output/p_hat_by_method_unmarked_full.csv"))
  readr::write_csv(tab_beta, file.path(out_dir, "occ_model_output/beta_by_method_unmarked_full.csv"))
  readr::write_csv(tab_psi,  file.path(out_dir, "occ_model_output/psi_hat_by_method_unmarked_full.csv"))
  
  
# option: un tableau de statut pour annexe
  tab_status <- res_ok %>%
    transmute(
      species,
      n_sites = map_int(result, "n_sites"),
      n_rep   = map_int(result, "n_rep"),
      AIC     = map_dbl(result, ~ .x$AIC %||% NA_real_),
      any_se_bad = map_lgl(result, ~ .x$any_se_bad %||% NA)
    ) %>%
    arrange(desc(any_se_bad), n_sites)
  

gt_status <- tab_status %>%
  gt() %>%
  tab_header(
    title = "Annexe — Occupancy models (unmarked::occu) : diagnostic per species" ) %>%
  fmt_number(columns = c(AIC), decimals = 2) %>%
  cols_label(
    species = "Species",
    n_sites = "Sites",
    n_rep = "Replicates",
    any_se_bad = "Unstable SE?",
    AIC = "AIC"
  ) %>%
  tab_style(
    style = cell_fill(),
    locations = cells_body(rows = any_se_bad == TRUE, columns = any_se_bad)
  )

gt::gtsave(gt_status, filename = file.path(out_dir, "occ_model_output/Annex_occu_status.html"))


gt_p <- tab_p %>%
  mutate(method = factor(method, levels = unique(method))) %>%
  gt(groupname_col = "species") %>%
  tab_header(
    title = "Annexe — Probabilités de détection (p̂) par méthode",
    subtitle = "Prédictions à effort médian (eff_z), IC 95%"
  ) %>%
  fmt_number(columns = c(p_hat, lcl, ucl), decimals = 3) %>%
  cols_label(
    method = "Method",
    p_hat  = html("p&#770;"),
    lcl    = "LCL",
    ucl    = "UCL",
    n_sites = "Sites",
    n_rep   = "Reps"
  ) %>%
  cols_hide(columns = c(n_sites, n_rep)) %>%
  tab_options(row_group.as_column = TRUE)

gt::gtsave(gt_p, filename = file.path(out_dir, "occ_model_output/Annex_p_hat_by_method.html"))


gt_beta <- tab_beta %>%
  # option: virer les termes de dummies si tu veux une annexe plus courte :
  # filter(!(component == "det" & term %in% unique(tab_p$method)))
  gt(groupname_col = "species") %>%
  tab_header(
    title = "Annexe — Coefficients des modèles (logit-scale)",
    subtitle = "State = ψ (occupation), Det = p (détection)"
  ) %>%
  fmt_number(columns = c(estimate, se, z, p), decimals = 3) %>%
  tab_options(row_group.as_column = TRUE)

gt::gtsave(gt_beta, filename = file.path(out_dir, "occ_model_output/Annex_occu_coefficients.html"))

