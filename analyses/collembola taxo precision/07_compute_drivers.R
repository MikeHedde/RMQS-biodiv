# ============================================================
# 07_compute_drivers — PERMANOVA et modèles de drivers alpha
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "C:/Users/Hedde/Documents/R/RMQS-biodiv/analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances
load_required_step("01_import_clean")
load_required_step("03_prepare_environment_drivers")
load_required_step("04_build_scenarios")
load_required_step("05_compute_alpha_metrics")

STEP_ID <- "07_compute_drivers"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("07_compute_drivers — PERMANOVA et modèles de drivers alpha : cache chargé")
} else {
  # 8. DRIVERS ENVIRONNEMENTAUX : PERMANOVA UNIVARIEE
  # -----------------------------

  message_header("Bloc drivers environnementaux")

  permanova_by_iter <- tibble()
  permanova_summary <- tibble()
  driver_rank_stability <- tibble()

  if (RUN_DRIVER_BLOCK && length(selected_drivers) > 0) {
  
    meta_unit <- if (ANALYSIS_UNIT == "station") {
      meta_station %>% transmute(unit = station, across(all_of(selected_drivers)))
    } else {
      col_dat %>%
        distinct(sample_id, station) %>%
        left_join(meta_station, by = "station") %>%
        transmute(unit = sample_id, across(all_of(selected_drivers)))
    }
  
    num_to_scale <- intersect(selected_numeric_drivers_driverblock, names(meta_unit))
    meta_unit <- meta_unit %>%
      mutate(across(all_of(num_to_scale), ~ as.numeric(scale(.x))))
  
    run_univariate_permanova_one <- function(df_one, scenario_name, scenario_family_name, baseline_name, iter_id) {
      mat <- make_comm_matrix(df_one, unit_col = ANALYSIS_UNIT)
      common <- intersect(rownames(mat), meta_unit$unit)
      mat <- mat[common, , drop = FALSE]
      meta <- meta_unit %>% filter(unit %in% common) %>% arrange(match(unit, rownames(mat)))
    
      purrr::map_dfr(selected_drivers, function(term_i) {
        complete <- complete.cases(meta %>% select(all_of(term_i)))
        mat_i <- mat[complete, , drop = FALSE]
        meta_i <- meta[complete, , drop = FALSE]
      
        keep_cols <- if (ncol(mat_i) > 0) colSums(mat_i) > 0 else logical(0)
        mat_i <- mat_i[, keep_cols, drop = FALSE]
      
        if (
          nrow(mat_i) < 10 ||
          ncol(mat_i) < 2 ||
          dplyr::n_distinct(meta_i[[term_i]], na.rm = TRUE) < 2
        ) {
          return(tibble(
            scenario_family = scenario_family_name,
            scenario = scenario_name,
            baseline_scenario = baseline_name,
            iter = iter_id,
            term = term_i,
            df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
            n_units = nrow(mat_i), n_taxa = ncol(mat_i)
          ))
        }
      
        meta_tmp <- tibble(x = meta_i[[term_i]])
        fit <- tryCatch(
          vegan::adonis2(mat_i ~ x, data = meta_tmp, method = "bray", permutations = PERMANOVA_N_PERM),
          error = function(e) NULL
        )
      
        if (is.null(fit)) {
          return(tibble(
            scenario_family = scenario_family_name,
            scenario = scenario_name,
            baseline_scenario = baseline_name,
            iter = iter_id,
            term = term_i,
            df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
            n_units = nrow(mat_i), n_taxa = ncol(mat_i)
          ))
        }
      
        as.data.frame(fit) %>%
          rownames_to_column("adonis_term") %>%
          filter(adonis_term == "x") %>%
          janitor::clean_names() %>%
          transmute(
            scenario_family = scenario_family_name,
            scenario = scenario_name,
            baseline_scenario = baseline_name,
            iter = iter_id,
            term = term_i,
            df = df,
            sumsq = sum_of_sqs,
            r2 = r2,
            f = f,
            p = pr_f,
            n_units = nrow(mat_i),
            n_taxa = ncol(mat_i)
          )
      })
    }
  
    deterministic_scenarios <- scenario_data %>%
      group_by(scenario_family, scenario) %>%
      summarise(n_iter = n_distinct(iter), .groups = "drop") %>%
      filter(n_iter == 1) %>%
      pull(scenario)
  
    driver_scenario_data <- scenario_data %>%
      filter(
        scenario %in% deterministic_scenarios |
          iter <= PERMANOVA_MAX_STOCHASTIC_ITER
      )
  
    permanova_by_iter <- driver_scenario_data %>%
      distinct(scenario_family, scenario, baseline_scenario, iter) %>%
      group_by(scenario_family, scenario, baseline_scenario, iter) %>%
      group_split() %>%
      purrr::map_dfr(function(key) {
        df <- driver_scenario_data %>%
          filter(
            scenario_family == key$scenario_family[1],
            scenario == key$scenario[1],
            baseline_scenario == key$baseline_scenario[1],
            iter == key$iter[1]
          )
        run_univariate_permanova_one(
          df_one = df,
          scenario_name = key$scenario[1],
          scenario_family_name = key$scenario_family[1],
          baseline_name = key$baseline_scenario[1],
          iter_id = key$iter[1]
        )
      })
  
    permanova_summary <- permanova_by_iter %>%
      group_by(scenario_family, scenario, baseline_scenario, term) %>%
      summarise(
        r2_mean = safe_mean(r2),
        r2_sd = safe_sd(r2),
        p_median = safe_median(p),
        p_prop_lt_005 = if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE),
        n_iter = n(),
        n_valid = sum(!is.na(r2)),
        .groups = "drop"
      )
  
    baseline_rank <- permanova_summary %>%
      filter(scenario == baseline_scenario, !is.na(r2_mean)) %>%
      group_by(scenario_family, baseline_scenario) %>%
      mutate(rank_baseline = min_rank(desc(r2_mean))) %>%
      ungroup() %>%
      select(scenario_family, baseline_scenario, term, rank_baseline, r2_baseline = r2_mean)
  
    driver_rank_stability <- permanova_summary %>%
      group_by(scenario_family, scenario, baseline_scenario) %>%
      mutate(rank_scenario = min_rank(desc(r2_mean))) %>%
      ungroup() %>%
      left_join(baseline_rank, by = c("scenario_family", "baseline_scenario", "term")) %>%
      group_by(scenario_family, scenario, baseline_scenario) %>%
      summarise(
        n_complete_rank_pairs = sum(complete.cases(rank_baseline, rank_scenario)),
        driver_rank_spearman = safe_spearman(rank_baseline, rank_scenario),
        top_driver = safe_top_driver(term, r2_mean),
        baseline_top_driver = safe_top_driver(term, r2_baseline),
        top_driver_changed = case_when(
          is.na(top_driver) | is.na(baseline_top_driver) ~ NA,
          TRUE ~ top_driver != baseline_top_driver
        ),
        .groups = "drop"
      )
  
    # PERMANOVA multivariée marginale : T360_mean + MOS + pH dans le même modèle.
    run_multivariate_permanova_one <- function(df_one, scenario_name, scenario_family_name, baseline_name, iter_id) {
      mat <- make_comm_matrix(df_one, unit_col = ANALYSIS_UNIT)
      common <- intersect(rownames(mat), meta_unit$unit)
      mat <- mat[common, , drop = FALSE]
      meta <- meta_unit %>% filter(unit %in% common) %>% arrange(match(unit, rownames(mat)))
    
      complete <- complete.cases(meta %>% select(all_of(selected_drivers)))
      mat_i <- mat[complete, , drop = FALSE]
      meta_i <- meta[complete, , drop = FALSE]
    
      keep_cols <- if (ncol(mat_i) > 0) colSums(mat_i) > 0 else logical(0)
      mat_i <- mat_i[, keep_cols, drop = FALSE]
    
      if (nrow(mat_i) < 10 || ncol(mat_i) < 2 || length(selected_drivers) < 1) {
        return(tibble(
          scenario_family = scenario_family_name,
          scenario = scenario_name,
          baseline_scenario = baseline_name,
          iter = iter_id,
          term = selected_drivers,
          df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
          n_units = nrow(mat_i), n_taxa = ncol(mat_i), model = "multivariate_margin"
        ))
      }
    
      form <- as.formula(paste("mat_i ~", paste(selected_drivers, collapse = " + ")))
      fit <- tryCatch(
        vegan::adonis2(form, data = meta_i, method = "bray", permutations = PERMANOVA_N_PERM, by = "margin"),
        error = function(e) NULL
      )
    
      if (is.null(fit)) {
        return(tibble(
          scenario_family = scenario_family_name,
          scenario = scenario_name,
          baseline_scenario = baseline_name,
          iter = iter_id,
          term = selected_drivers,
          df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
          n_units = nrow(mat_i), n_taxa = ncol(mat_i), model = "multivariate_margin"
        ))
      }
    
      as.data.frame(fit) %>%
        rownames_to_column("adonis_term") %>%
        filter(adonis_term %in% selected_drivers) %>%
        janitor::clean_names() %>%
        transmute(
          scenario_family = scenario_family_name,
          scenario = scenario_name,
          baseline_scenario = baseline_name,
          iter = iter_id,
          term = adonis_term,
          df = df,
          sumsq = sum_of_sqs,
          r2 = r2,
          f = f,
          p = pr_f,
          n_units = nrow(mat_i),
          n_taxa = ncol(mat_i),
          model = "multivariate_margin"
        )
    }
  
    permanova_multivar_by_iter <- driver_scenario_data %>%
      distinct(scenario_family, scenario, baseline_scenario, iter) %>%
      group_by(scenario_family, scenario, baseline_scenario, iter) %>%
      group_split() %>%
      purrr::map_dfr(function(key) {
        df <- driver_scenario_data %>%
          filter(
            scenario_family == key$scenario_family[1],
            scenario == key$scenario[1],
            baseline_scenario == key$baseline_scenario[1],
            iter == key$iter[1]
          )
        run_multivariate_permanova_one(
          df_one = df,
          scenario_name = key$scenario[1],
          scenario_family_name = key$scenario_family[1],
          baseline_name = key$baseline_scenario[1],
          iter_id = key$iter[1]
        )
      })
  
    permanova_multivar_summary <- permanova_multivar_by_iter %>%
      group_by(scenario_family, scenario, baseline_scenario, term) %>%
      summarise(
        r2_mean = safe_mean(r2),
        r2_sd = safe_sd(r2),
        p_median = safe_median(p),
        p_prop_lt_005 = if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE),
        n_iter = n(),
        n_valid = sum(!is.na(r2)),
        .groups = "drop"
      )
  
    # Modèles alpha : stabilité des coefficients directionnels pour q0, q1, q2 et couverture.
    alpha_driver_input <- alpha_diversity_by_unit %>%
      left_join(meta_unit, by = "unit") %>%
      pivot_longer(
        cols = c(q0, q1, q2, coverage_chao),
        names_to = "alpha_metric",
        values_to = "alpha_value"
      )
  
    fit_alpha_driver_one <- function(df) {
      fam <- df$scenario_family[1]
      scen <- df$scenario[1]
      base <- df$baseline_scenario[1]
      iter_i <- df$iter[1]
      metric_i <- df$alpha_metric[1]
    
      dd <- df %>%
        select(alpha_value, all_of(selected_drivers)) %>%
        mutate(y = log1p(alpha_value)) %>%
        filter(complete.cases(.))
    
      if (nrow(dd) < 10 || n_distinct(dd$y) < 3 || length(selected_drivers) < 1) {
        return(tibble(
          scenario_family = fam, scenario = scen, baseline_scenario = base,
          iter = iter_i, alpha_metric = metric_i,
          term = selected_drivers, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p = NA_real_,
          adj_r2 = NA_real_, n_units = nrow(dd)
        ))
      }
    
      form <- as.formula(paste("y ~", paste(selected_drivers, collapse = " + ")))
      fit <- tryCatch(lm(form, data = dd), error = function(e) NULL)
      if (is.null(fit)) {
        return(tibble(
          scenario_family = fam, scenario = scen, baseline_scenario = base,
          iter = iter_i, alpha_metric = metric_i,
          term = selected_drivers, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p = NA_real_,
          adj_r2 = NA_real_, n_units = nrow(dd)
        ))
      }
    
      adj_r2 <- broom::glance(fit)$adj.r.squared[1]
      broom::tidy(fit) %>%
        filter(term %in% selected_drivers) %>%
        transmute(
          scenario_family = fam, scenario = scen, baseline_scenario = base,
          iter = iter_i, alpha_metric = metric_i,
          term = term,
          estimate = estimate,
          std_error = std.error,
          statistic = statistic,
          p = p.value,
          adj_r2 = adj_r2,
          n_units = nrow(dd)
        )
    }
  
    alpha_driver_by_iter <- alpha_driver_input %>%
      filter(scenario %in% deterministic_scenarios | iter <= PERMANOVA_MAX_STOCHASTIC_ITER) %>%
      group_by(scenario_family, scenario, baseline_scenario, iter, alpha_metric) %>%
      group_split() %>%
      purrr::map_dfr(fit_alpha_driver_one)
  
    alpha_driver_summary <- alpha_driver_by_iter %>%
      group_by(scenario_family, scenario, baseline_scenario, alpha_metric, term) %>%
      summarise(
        estimate_mean = safe_mean(estimate),
        estimate_sd = safe_sd(estimate),
        p_median = safe_median(p),
        p_prop_lt_005 = if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE),
        adj_r2_mean = safe_mean(adj_r2),
        n_iter = n(),
        n_valid = sum(!is.na(estimate)),
        .groups = "drop"
      )
  
    alpha_driver_baseline <- alpha_driver_summary %>%
      filter(scenario == baseline_scenario) %>%
      transmute(
        scenario_family, baseline_scenario, alpha_metric, term,
        baseline_estimate = estimate_mean,
        baseline_abs_estimate_rank = min_rank(desc(abs(estimate_mean)))
      )
  
    alpha_driver_stability <- alpha_driver_summary %>%
      left_join(alpha_driver_baseline, by = c("scenario_family", "baseline_scenario", "alpha_metric", "term")) %>%
      mutate(
        estimate_sign_changed = case_when(
          is.na(estimate_mean) | is.na(baseline_estimate) ~ NA,
          estimate_mean == 0 | baseline_estimate == 0 ~ FALSE,
          TRUE ~ sign(estimate_mean) != sign(baseline_estimate)
        ),
        estimate_ratio = if_else(!is.na(baseline_estimate) & baseline_estimate != 0, estimate_mean / baseline_estimate, NA_real_)
      )
  
    readr::write_csv(permanova_multivar_by_iter, file.path(OUT_DIR, "permanova_multivar_by_iter.csv"))
    readr::write_csv(permanova_multivar_summary, file.path(OUT_DIR, "permanova_multivar_summary.csv"))
    readr::write_csv(alpha_driver_by_iter, file.path(OUT_DIR, "alpha_driver_by_iter.csv"))
    readr::write_csv(alpha_driver_summary, file.path(OUT_DIR, "alpha_driver_summary.csv"))
    readr::write_csv(alpha_driver_stability, file.path(OUT_DIR, "alpha_driver_stability.csv"))
  
    readr::write_csv(permanova_by_iter, file.path(OUT_DIR, "permanova_by_iter.csv"))
    readr::write_csv(permanova_summary, file.path(OUT_DIR, "permanova_summary.csv"))
    readr::write_csv(driver_rank_stability, file.path(OUT_DIR, "driver_rank_stability.csv"))
  
    message("PERMANOVA univariée : ", nrow(permanova_by_iter), " lignes écrites.")
    message("PERMANOVA multivariée : ", nrow(permanova_multivar_by_iter), " lignes écrites.")
    message("Modèles alpha : ", nrow(alpha_driver_by_iter), " lignes écrites.")
  
  } else {
    warning("Bloc drivers ignoré : aucun driver sélectionné ou RUN_DRIVER_BLOCK = FALSE.")
    readr::write_csv(permanova_by_iter, file.path(OUT_DIR, "permanova_by_iter.csv"))
    readr::write_csv(permanova_summary, file.path(OUT_DIR, "permanova_summary.csv"))
    readr::write_csv(driver_rank_stability, file.path(OUT_DIR, "driver_rank_stability.csv"))
    readr::write_csv(tibble(), file.path(OUT_DIR, "permanova_multivar_by_iter.csv"))
    readr::write_csv(tibble(), file.path(OUT_DIR, "permanova_multivar_summary.csv"))
    readr::write_csv(tibble(), file.path(OUT_DIR, "alpha_driver_by_iter.csv"))
    readr::write_csv(tibble(), file.path(OUT_DIR, "alpha_driver_summary.csv"))
    readr::write_csv(tibble(), file.path(OUT_DIR, "alpha_driver_stability.csv"))
  }


  # Sécurise les objets optionnels si le bloc drivers est ignoré.
  .optional_driver_objects <- c(
    "meta_unit",
    "deterministic_scenarios",
    "driver_scenario_data",
    "permanova_multivar_by_iter",
    "permanova_multivar_summary",
    "alpha_driver_input",
    "alpha_driver_by_iter",
    "alpha_driver_summary",
    "alpha_driver_baseline",
    "alpha_driver_stability"
  )
  for (.obj in .optional_driver_objects) {
    if (!exists(.obj, inherits = FALSE)) assign(.obj, tibble())
  }

  save_step(STEP_ID, c(
    "permanova_by_iter",
    "permanova_summary",
    "driver_rank_stability",
    "meta_unit",
    "deterministic_scenarios",
    "driver_scenario_data",
    "permanova_multivar_by_iter",
    "permanova_multivar_summary",
    "alpha_driver_input",
    "alpha_driver_by_iter",
    "alpha_driver_summary",
    "alpha_driver_baseline",
    "alpha_driver_stability"
  ))
}
