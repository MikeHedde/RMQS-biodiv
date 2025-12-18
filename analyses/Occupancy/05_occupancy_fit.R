# ============================================================
# MODÈLES unmarked::occu ESPÈCE-PAR-ESPÈCE
# Sorties fichiers: tab_p, tab_beta, tab_psi + CSV
# ============================================================

method_eff_medians <- function(eff_z_matrix)
  apply(eff_z_matrix, 2, function(col) if (all(is.na(col))) 0 else median(col, na.rm=TRUE))

fit_one_species_occu <- function(sp_name,
                                 min_sites = 2,
                                 require_contrast = TRUE) {
  #--------------------------------------------------
  # 0) Identifier l'espèce dans Y
  #--------------------------------------------------
  k <- which(dimnames(Y)$species == sp_name)
  if (length(k) != 1) return(NULL)
  
  Yk <- Y[,,k, drop = FALSE][,,1]   # matrice J × R
  # rownames(Yk) = sites ; colnames(Yk) = réplicats
  
  #--------------------------------------------------
  # 1) Filtrer sites et réplicats informatifs
  #--------------------------------------------------
  # Sites où au moins une obs non-NA
  keep_site <- apply(Yk, 1, function(a) any(!is.na(a)))
  Yk <- Yk[keep_site, , drop = FALSE]
  if (nrow(Yk) < min_sites) return(NULL)
  
  # Réplicats où au moins une obs non-NA
  keep_rep <- apply(Yk, 2, function(a) any(!is.na(a)))
  Yk <- Yk[, keep_rep, drop = FALSE]
  if (ncol(Yk) < 2) return(NULL)   # au moins 2 réplicats pour séparer ψ / p
  
  # Contraste (au moins une détection et une non-détection)
  yvec <- as.vector(Yk[!is.na(Yk)])
  if (require_contrast && !(any(yvec == 1) && any(yvec == 0))) return(NULL)
  
  rep_names_k <- colnames(Yk)
  
  #--------------------------------------------------
  # 2) Matrices d’indicateurs de méthode + effort
  #    (obsCovs pour la détection)
  #--------------------------------------------------
  mkI <- function(tag) {
    matrix(as.numeric(rep_names_k == tag),
           nrow = nrow(Yk), ncol = ncol(Yk), byrow = TRUE)
  }
  
  # Une matrice par réplicat (Pitfall10, Pitfall8, GPD, DVAC, TM, etc.)
  I_list <- setNames(lapply(rep_names_k, mkI), rep_names_k)
  
  # Effort aligné sur Yk
  eff_sp <- eff_z[rownames(Yk), rep_names_k, drop = FALSE]
  
  # Masquer là où Y est NA
  for (nm in rep_names_k) {
    I_list[[nm]][is.na(Yk)] <- NA
  }
  eff_sp[is.na(Yk)] <- NA
  
  #--------------------------------------------------
  # 3) Covariables de site (ψ et éventuellement p)
  #--------------------------------------------------
  sc <- site_cov %>%
    dplyr::filter(site_id %in% rownames(Yk)) %>%
    dplyr::arrange(match(site_id, rownames(Yk)))
  
  # Site covariates passés à unmarked (sans site_id)
  sc_umf <- sc %>% dplyr::select(-site_id)
  rownames(sc_umf) <- sc$site_id
  
  #--------------------------------------------------
  # 4) Formules : ψ ~ DOY_z ; p ~ méthode + effort_z (+ HAB si dispo)
  #--------------------------------------------------
  # Occupancy (ψ) : DOY_z ou intercept
  state_terms <- c()
  state_terms <- c(state_terms, "DOY_z", "ALTITUDE_z", "NDVI_m_z", "NDVI_sd_z")
  state_rhs <- if (length(state_terms)) paste(state_terms, collapse = " + ") else "1"
  
  # Détection (p) : indicateurs de méthode + effort + éventuellement HAB
  det_terms <- c(rep_names_k, "eff_z")

  det_rhs <- paste(det_terms, collapse = " + ")
  
  form_occu <- as.formula(paste("~", det_rhs, "~", state_rhs))
  
  #--------------------------------------------------
  # 5) unmarkedFrameOccu + fit du modèle
  #--------------------------------------------------
  obsCovs <- c(I_list, list(eff_z = eff_sp))
  
  umf <- unmarked::unmarkedFrameOccu(
    y       = Yk,
    siteCovs = sc_umf,
    obsCovs  = obsCovs
  )
  
  fm <- try(
    unmarked::occu(form_occu, data = umf,
                   control = list(maxit = 200)),
    silent = TRUE
  )
  if (inherits(fm, "try-error")) return(NULL)
  
  #--------------------------------------------------
  # 6) Prédictions de p̂ par méthode
  #     à effort médian global, DOY_z = 0, HAB = ref
  #--------------------------------------------------
  method_eff_medians <- function(eff_z_matrix) {
    apply(eff_z_matrix, 2, function(col) {
      if (all(is.na(col))) 0 else stats::median(col, na.rm = TRUE)
    })
  }
  eff_med_all <- method_eff_medians(eff_z)
  
  eff_for <- function(m) {
    if (m %in% names(eff_med_all)) eff_med_all[[m]] else 0
  }

  mk_pred <- function(tag) {
    nd <- data.frame(eff_z = eff_for(tag))
    # indicateurs de méthode
    for (nm in rep_names_k) {
      nd[[nm]] <- as.numeric(nm == tag)
    }
    # DOY_z fixé à 0 si présent dans le modèle
    if (has_doy) nd$DOY_z <- 0
    if (has_altitude) nd$ALTITUDE_z <- 0
    if (has_ndvi_m) nd$NDVI_m_z <- 0
    if (has_ndvi_sd) nd$NDVI_sd_z <- 0
    
    as.data.frame(unmarked::predict(fm, type = "det", newdata = nd))
  }
  
  pred_list <- lapply(rep_names_k, mk_pred)
  
  tab_p <- tibble::tibble(
    species = sp_name,
    method  = rep_names_k,
    p_hat   = vapply(pred_list, function(z) z$Predicted, numeric(1)),
    lcl     = vapply(pred_list, function(z) z$lower,     numeric(1)),
    ucl     = vapply(pred_list, function(z) z$upper,     numeric(1)),
    n_sites = nrow(Yk),
    n_rep   = ncol(Yk)
  )
  
  #--------------------------------------------------
  # 7) Coefficients d’occupation (ψ) pour DOY_z
  #--------------------------------------------------
  co <- summary(fm)
  co_state <- co$state
  
  get_or_na <- function(mat, par, col) {
    if (!(is.matrix(mat) || is.data.frame(mat))) return(NA_real_)
    if (!(par %in% rownames(mat))) return(NA_real_)
    mat[par, col]
  }
  
  beta_DOY <- get_or_na(co_state, "DOY_z", "Estimate")
  se_DOY   <- get_or_na(co_state, "DOY_z", "SE")
  
  beta_ALTITUDE <- get_or_na(co_state, "ALTITUDE_z", "Estimate")
  se_ALTITUDE   <- get_or_na(co_state, "ALTITUDE_z", "SE")
  
  beta_NDVI_m <- get_or_na(co_state, "NDVI_m_z", "Estimate")
  se_NDVI_m   <- get_or_na(co_state, "NDVI_m_z", "SE")
  
  beta_NDVI_sd <- get_or_na(co_state, "NDVI_sd_z", "Estimate")
  se_NDVI_sd   <- get_or_na(co_state, "NDVI_sd_z", "SE")
  
  beta_tab <- tibble::tibble(
    species  = sp_name,
    beta_doy = beta_DOY,
    se_doy   = se_DOY,
    beta_altitude = beta_ALTITUDE,
    se_altitude = se_ALTITUDE,
    beta_ndvi_m = beta_NDVI_m,
    se_ndvi_m = se_NDVI_m,
    beta_ndvi_sd = beta_NDVI_sd,
    se_ndvi_sd = se_NDVI_sd
    )
  
  #--------------------------------------------------
  # 8) ψ̂ par site (richesse attendue ensuite via Σ ψ̂)
  #--------------------------------------------------
  psi_pred <- as.data.frame(unmarked::predict(fm, type = "state"))
  psi_hat  <- psi_pred$Predicted
  names(psi_hat) <- rownames(Yk)
  
  psi_tab <- tibble::tibble(
    site_id = rownames(Yk),
    species = sp_name,
    psi_hat = psi_hat
  )
  
  #--------------------------------------------------
  # 9) Sortie
  #--------------------------------------------------
  list(
    p_table = tab_p,
    beta    = beta_tab,
    psi_hat = psi_tab
  )
}

## ---------- Boucle sur toutes les espèces + exports ----------
spp_all  <- dimnames(Y)$species

res_list <- purrr::map(spp_all, ~ try(fit_one_species_occu(.x), silent = FALSE))

tab_p    <- res_list %>%
  purrr::keep(~ !inherits(.x,"try-error") && !is.null(.x$p_table)) %>%
  purrr::map("p_table") %>% dplyr::bind_rows()

tab_beta <- res_list %>%
  purrr::keep(~ !inherits(.x,"try-error") && !is.null(.x$beta)) %>%
  purrr::map("beta") %>% dplyr::bind_rows()

tab_psi  <- res_list %>%
  purrr::keep(~ !inherits(.x,"try-error") && !is.null(.x$psi_hat)) %>%
  purrr::map("psi_hat") %>% dplyr::bind_rows()

readr::write_csv(tab_p,    file.path(out_dir, "p_hat_by_method_unmarked_full.csv"))
readr::write_csv(tab_beta, file.path(out_dir, "psi_coef_cov_by_species.csv"))
readr::write_csv(tab_psi,  file.path(out_dir, "psi_hat_by_site_species.csv"))
