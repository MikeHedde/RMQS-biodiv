# ============================================================
# ANNEXE — DIAGNOSTICS LMM/GLMM (robuste)
# Requiert: m1, tab_p2, out_dir
# ============================================================

pkgs <- c("DHARMa", "performance", "broom.mixed", "dplyr", "purrr", "readr")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(DHARMa)
library(performance)
library(broom.mixed)
library(dplyr)
library(purrr)
library(readr)

set.seed(1)

safe_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

# ---- 1) DHARMa diagnostics: PDF (base) + table tests ----
run_dharma_bundle <- function(model, model_name, out_dir, nsim = 1000,
                              width = 10, height = 7) {
  
  model_slug <- safe_slug(model_name)
  
  sim <- DHARMa::simulateResiduals(fittedModel = model, n = nsim, plot = FALSE)
  
  # Tests principaux
  t_uni <- DHARMa::testUniformity(sim)
  t_dis <- DHARMa::testDispersion(sim)
  t_out <- DHARMa::testOutliers(sim)
  
  tab_tests <- dplyr::bind_rows(
    tibble(test = "uniformity", statistic = unname(t_uni$statistic), p.value = unname(t_uni$p.value)),
    tibble(test = "dispersion", statistic = unname(t_dis$statistic), p.value = unname(t_dis$p.value)),
    tibble(test = "outliers",   statistic = unname(t_out$statistic), p.value = unname(t_out$p.value))
  ) %>%
    mutate(model = model_name) %>%
    select(model, test, statistic, p.value)
  
  # Export tests
  csv_path <- file.path(out_dir, paste0("DARMHa/TabS_DHARMa_tests_", model_slug, ".csv"))
  readr::write_csv(tab_tests, csv_path)
  
  # Export figure (base plots)
  pdf_path <- file.path(out_dir, paste0("DARMHa/FigS_DHARMa_", model_slug, ".pdf"))
  grDevices::pdf(pdf_path, width = width, height = height)
  op <- par(mfrow = c(2, 2), mar = c(4,4,2,1))
  
  # 1) plot global (remplace plotSimulatedResiduals, désormais deprecated)
  plot(sim, main = "DHARMa: simulated residual diagnostics")
  
  # 2) QQ-uniform
  DHARMa::plotQQunif(sim, main = "DHARMa: QQ-uniform")
  
  # 3) residuals vs predicted
  DHARMa::plotResiduals(sim, main = "DHARMa: residuals vs predicted")
  
  # 4) un panneau texte avec les p-values
  plot.new()
  txt <- paste0(
    "DHARMa tests (", nsim, " simulations)\n\n",
    "Uniformity p = ", signif(unname(t_uni$p.value), 3), "\n",
    "Dispersion p = ", signif(unname(t_dis$p.value), 3), "\n",
    "Outliers   p = ", signif(unname(t_out$p.value), 3), "\n"
  )
  text(0, 1, adj = c(0,1), labels = txt, cex = 1)
  
  par(op)
  grDevices::dev.off()
  
  message("✅ DHARMa figure exported: ", pdf_path)
  message("✅ DHARMa tests exported:  ", csv_path)
  
  invisible(list(sim = sim, tests = tab_tests, fig_path = pdf_path, tab_path = csv_path))
}

# ---- 2) performance bundle: table robuste ----
run_performance_bundle <- function(model, model_name, out_dir) {
  
  model_slug <- safe_slug(model_name)
  
  sing <- performance::check_singularity(model)
  singular_flag <- if (is.list(sing) && "singular" %in% names(sing)) isTRUE(sing$singular) else isTRUE(sing)
  
  r2v  <- performance::r2(model)
  iccv <- performance::icc(model)
  
  # calculs hors tibble (évite le masquage)
  aic  <- AIC(model)
  bic  <- BIC(model)
  ll   <- as.numeric(stats::logLik(model))
  
  tab_perf <- tibble::tibble(
    model_id = model_name,
    singular = singular_flag,
    r2_marginal    = if ("R2_marginal" %in% names(r2v)) r2v$R2_marginal else NA_real_,
    r2_conditional = if ("R2_conditional" %in% names(r2v)) r2v$R2_conditional else NA_real_,
    icc = dplyr::coalesce(
      iccv$ICC_adjusted %||% NA_real_,
      iccv$ICC %||% NA_real_
    ),
    AIC = aic,
    BIC = bic,
    logLik = ll
  )
  
  csv_path <- file.path(out_dir, paste0("DARMHa/TabS_performance_", model_slug, ".csv"))
  readr::write_csv(tab_perf, csv_path)
  
  message("✅ performance table exported: ", csv_path)
  invisible(tab_perf)
}


`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- 3) Leave-one-species-out sensitivity (fix capture) ----
run_leave_one_species_out <- function(data, formula,
                                      weights_col = NULL,
                                      species_var = "species",
                                      model_name = "m1",
                                      out_dir,
                                      lmer_fun = lme4::lmer) {
  
  sv <- species_var
  model_slug <- safe_slug(paste0(model_name, "DARMHa/_LOO"))
  sp_levels <- sort(unique(data[[sv]]))
  
  # gabarit de sortie en cas d'erreur (colonnes attendues)
  err_row <- function(sp, msg) {
    tibble::tibble(
      term      = NA_character_,
      estimate  = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value   = NA_real_,
      left_out  = sp,
      ok        = FALSE,
      message   = msg
    )
  }
  
  results <- vector("list", length(sp_levels))
  names(results) <- sp_levels
  
  for (sp in sp_levels) {
    
    dat <- dplyr::filter(data, .data[[sv]] != sp)
    dat <- droplevels(dat)
    
    w_vec <- if (!is.null(weights_col)) dat[[weights_col]] else NULL
    
    fit <- tryCatch({
      
      mod <- if (is.null(w_vec)) {
        lmer_fun(formula, data = dat)
      } else {
        lmer_fun(formula, data = dat, weights = w_vec)
      }
      
      td <- broom.mixed::tidy(mod, effects = "fixed") |>
        dplyr::filter(term != "(Intercept)") |>
        dplyr::mutate(
          left_out = sp,
          ok = TRUE,
          message = NA_character_
        )
      
      # au cas où tidy renverrait 0 ligne (rare mais possible)
      if (nrow(td) == 0) err_row(sp, "tidy() returned 0 fixed-effect rows") else td
      
    }, error = function(e) {
      err_row(sp, conditionMessage(e))
    })
    
    results[[sp]] <- fit
  }
  
  loo <- dplyr::bind_rows(results)
  
  # Log des échecs
  log_df <- loo |>
    dplyr::filter(ok == FALSE) |>
    dplyr::select(left_out, message)
  
  log_path <- file.path(out_dir, paste0("DARMHA/TabS_LOO_log_", model_slug, ".csv"))
  readr::write_csv(log_df, log_path)
  
  # Garder les fits OK
  loo_ok <- loo |>
    dplyr::filter(ok == TRUE)
  
  # Si tout a échoué, on sort proprement
  if (nrow(loo_ok) == 0) {
    warning("LOO: no successful refits. See log: ", log_path)
    return(invisible(list(all = loo_ok, summary = NULL, log = log_df, log_path = log_path)))
  }
  
  loo_summary <- loo_ok |>
    dplyr::group_by(term) |>
    dplyr::summarise(
      est_min   = min(estimate, na.rm = TRUE),
      est_max   = max(estimate, na.rm = TRUE),
      est_range = est_max - est_min,
      .groups   = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(est_range))
  
  # Export
  csv_all <- file.path(out_dir, paste0("DARMHa/TabS_LOO_all_", model_slug, ".csv"))
  csv_sum <- file.path(out_dir, paste0("DARMHa/TabS_LOO_summary_", model_slug, ".csv"))
  readr::write_csv(loo_ok, csv_all)
  readr::write_csv(loo_summary, csv_sum)
  
  message("✅ LOO OK refits: ", length(unique(loo_ok$left_out)))
  message("✅ LOO exported:  ", csv_all)
  message("✅ LOO summary:   ", csv_sum)
  message("ℹ️  LOO log:      ", log_path)
  
  invisible(list(all = loo_ok, summary = loo_summary, log = log_df, log_path = log_path))
}


# ============================================================
# RUN (pour ton m1)
# ============================================================

d_m1 <- run_dharma_bundle(
  model = m1,
  model_name = "m1_logit_p_by_method_weighted_random_species",
  out_dir = out_dir,
  nsim = 1000
)

p_m1 <- run_performance_bundle(m1, "m1_logit_p_by_method_weighted_random_species", out_dir)

loo_m1 <- run_leave_one_species_out(
  data = tab_p2,
  formula = logit_p ~ method + (1 | species),
  weights_col = "weight",
  species_var = "species",
  model_name = "m1",
  out_dir = out_dir
)


# ============================================================
# OPTIONAL: do the same for m_trend and m_gpd if you want
# ============================================================
# run_dharma_bundle(m_trend, "m_trend_ntrap", out_dir)
# run_performance_bundle(m_trend, "m_trend_ntrap", out_dir)
# run_dharma_bundle(m_gpd, "m_gpd_vs_pitfall6", out_dir)
# run_performance_bundle(m_gpd, "m_gpd_vs_pitfall6", out_dir)
