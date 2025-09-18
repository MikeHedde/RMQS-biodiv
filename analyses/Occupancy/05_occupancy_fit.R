# ============================================================
# MODÈLES unmarked::occu ESPÈCE-PAR-ESPÈCE
# Sorties fichiers: tab_p, tab_beta, tab_psi + CSV
# ============================================================

method_eff_medians <- function(eff_z_matrix)
  apply(eff_z_matrix, 2, function(col) if (all(is.na(col))) 0 else median(col, na.rm=TRUE))

fit_one_species_occu <- function(sp_name) {
  k <- which(dimnames(Y)$species == sp_name); if (length(k) != 1) return(NULL)
  Yk <- Y[,,k, drop=FALSE][,,1]
  keep_site <- apply(Yk, 1, function(a) any(!is.na(a))); Yk <- Yk[keep_site,, drop=FALSE]
  if (nrow(Yk) < min_sites) return(NULL)
  keep_rep  <- apply(Yk, 2, function(a) any(!is.na(a))); Yk <- Yk[, keep_rep, drop=FALSE]
  if (ncol(Yk) < 2) return(NULL)
  yvec <- as.vector(Yk[!is.na(Yk)]); if (require_contrast && !(any(yvec==1) && any(yvec==0))) return(NULL)
  
  rep_names_k <- colnames(Yk)
  mkI <- function(tag) matrix(as.numeric(rep_names_k == tag), nrow(Yk), ncol(Yk), byrow=TRUE)
  I_list <- setNames(lapply(rep_names_k, mkI), rep_names_k)
  eff_sp <- eff_z[rownames(Yk), rep_names_k, drop=FALSE]
  for (nm in rep_names_k) I_list[[nm]][is.na(Yk)] <- NA
  eff_sp[is.na(Yk)] <- NA
  
  sc <- site_cov %>% filter(site_id %in% rownames(Yk)) %>% arrange(match(site_id, rownames(Yk))) %>% select(ALTITUDE_z, DOY_z)
  has_var_alt <- sd(sc$ALTITUDE_z, na.rm=TRUE) > 0; has_var_doy <- sd(sc$DOY_z, na.rm=TRUE) > 0
  state_terms <- c(); if (has_var_alt) state_terms <- c(state_terms,"ALTITUDE_z"); if (has_var_doy) state_terms <- c(state_terms,"DOY_z")
  state_rhs <- if (length(state_terms)) paste(state_terms, collapse=" + ") else "1"
  
  obsCovs <- c(I_list, list(eff_z = eff_sp))
  umf <- unmarked::unmarkedFrameOccu(y = Yk, siteCovs = sc, obsCovs = obsCovs)
  
  det_terms <- paste(c(rep_names_k, "eff_z"), collapse = " + ")
  form_occu <- as.formula(paste("~", det_terms, "~", state_rhs))
  
  fm <- try(unmarked::occu(form_occu, data = umf, control = list(maxit = 200)), silent = TRUE)
  if (inherits(fm, "try-error")) return(NULL)
  
  eff_med_all <- method_eff_medians(eff_z)
  eff_for <- function(m) if (m %in% names(eff_med_all)) eff_med_all[[m]] else 0
  mk_pred <- function(tag) {
    nd <- data.frame(eff_z = eff_for(tag))
    for (nm in rep_names_k) nd[[nm]] <- as.numeric(nm == tag)
    if (has_var_alt) nd$ALTITUDE_z <- 0
    if (has_var_doy) nd$DOY_z      <- 0
    as.data.frame(unmarked::predict(fm, type = "det", newdata = nd))
  }
  pred_list <- lapply(rep_names_k, mk_pred)
  
  co_state <- summary(fm)$state
  getcoef <- function(nm) if (nm %in% rownames(co_state)) co_state[nm,"Estimate"] else NA_real_
  getse   <- function(nm) if (nm %in% rownames(co_state)) co_state[nm,"SE"]       else NA_real_
  
  list(
    p_table = tibble::tibble(
      species = sp_name, method = rep_names_k,
      p_hat = vapply(pred_list, function(z) z$Predicted, numeric(1)),
      lcl   = vapply(pred_list, function(z) z$lower,     numeric(1)),
      ucl   = vapply(pred_list, function(z) z$upper,     numeric(1)),
      n_sites = nrow(Yk), n_rep = ncol(Yk)
    ),
    beta = tibble::tibble(
      species  = sp_name,
      beta_alt = getcoef("ALTITUDE_z"), se_alt = getse("ALTITUDE_z"),
      beta_doy = getcoef("DOY_z"),      se_doy = getse("DOY_z")
    ),
    psi_hat = tibble::tibble(
      site_id = rownames(Yk), species = sp_name,
      psi_hat = as.data.frame(unmarked::predict(fm, type="state"))$Predicted
    )
  )
}

spp_all  <- dimnames(Y)$species
res_list <- purrr::map(spp_all, ~try(fit_one_species_occu(.x), silent = TRUE))

tab_p   <- res_list %>% keep(~!inherits(.x,"try-error") && !is.null(.x$p_table)) %>% map("p_table") %>% bind_rows()
tab_beta<- res_list %>% keep(~!inherits(.x,"try-error") && !is.null(.x$beta))     %>% map("beta")     %>% bind_rows()
tab_psi <- res_list %>% keep(~!inherits(.x,"try-error") && !is.null(.x$psi_hat))  %>% map("psi_hat")  %>% bind_rows()

readr::write_csv(tab_p,   file.path(out_dir, "p_hat_by_method_unmarked_full.csv"))
readr::write_csv(tab_beta,file.path(out_dir, "psi_coef_ALT_DOY_by_species.csv"))
readr::write_csv(tab_psi, file.path(out_dir, "psi_hat_by_site_species.csv"))
