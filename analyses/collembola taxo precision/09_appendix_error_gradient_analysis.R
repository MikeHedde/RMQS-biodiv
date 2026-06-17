# ============================================================
# 09_appendix_error_gradient_analysis — gradient 1% à 20% d'erreur d'identification
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "C:/Users/Hedde/Documents/R/RMQS-biodiv/analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances minimales
load_required_step("01_import_clean")
load_required_step("02_prepare_taxonomy_pools")

STEP_ID <- "09_appendix_error_gradient_analysis"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("09_appendix_error_gradient_analysis — gradient 1% à 20% : cache chargé")
} else {
  message_header("Analyse annexe — gradient d'erreur d'identification de 1% à 20%")

  appendix_dir <- file.path(OUT_DIR, "appendix_error_gradient")
  dir.create(appendix_dir, showWarnings = FALSE, recursive = TRUE)

  appendix_pool_name <- match.arg(ID_ERROR_APPENDIX_POOL, choices = c("observed", "taxref_mainland"))
  appendix_pool_tbl <- switch(
    appendix_pool_name,
    observed = observed_same_genus_pool,
    taxref_mainland = taxref_same_genus_pool
  )

  appendix_error_grid <- tibble(
    error_rate = ID_ERROR_APPENDIX_RATES,
    error_pct = round(100 * error_rate),
    scenario = paste0("appendix_species_same_genus_", appendix_pool_name, "_", sprintf("%02d", round(100 * error_rate)), "pct")
  )

  base_cols <- c("station", "repetition", "sample_id", "abundance", "taxo_level", "binomial", "species_key", "species_label", "genus", "famille", "ordre")

  choose_same_genus_candidate <- function(genus, current_species_key, pool_tbl) {
    cand <- pool_tbl %>%
      filter(.data$genus == !!genus, .data$candidate_key != !!current_species_key) %>%
      distinct(candidate_key, candidate_label, genus)

    if (nrow(cand) == 0) {
      return(tibble(candidate_key = NA_character_, candidate_label = NA_character_, genus = NA_character_))
    }

    cand %>% slice_sample(n = 1)
  }

  simulate_same_genus_error <- function(dat_species, pool_tbl, scenario_name, iter_id, fixed_rate) {
    dat0 <- dat_species %>%
      filter(taxo_level == "species", !is.na(species_key), !is.na(genus)) %>%
      select(all_of(base_cols)) %>%
      mutate(p_error = fixed_rate)

    purrr::pmap_dfr(dat0, function(station, repetition, sample_id, abundance, taxo_level,
                                   binomial, species_key, species_label, genus, famille, ordre, p_error) {
      n <- as.integer(round(abundance))
      if (is.na(p_error) || p_error <= 0 || n <= 0) {
        return(tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = abundance, taxo_level = taxo_level, binomial = binomial,
          species_key = species_key, species_label = species_label,
          genus = genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", species_key)
        ))
      }

      target <- choose_same_genus_candidate(genus, species_key, pool_tbl)
      if (nrow(target) == 0 || is.na(target$candidate_key[1])) {
        return(tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = abundance, taxo_level = taxo_level, binomial = binomial,
          species_key = species_key, species_label = species_label,
          genus = genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", species_key)
        ))
      }

      target_key <- target$candidate_key[1]
      target_label <- target$candidate_label[1]
      target_genus <- target$genus[1]

      n_swap <- rbinom(1, size = n, prob = p_error)
      if (n_swap == 0) {
        return(tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = abundance, taxo_level = taxo_level, binomial = binomial,
          species_key = species_key, species_label = species_label,
          genus = genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", species_key)
        ))
      }

      bind_rows(
        if (n - n_swap > 0) {
          tibble(
            station = station, repetition = repetition, sample_id = sample_id,
            abundance = n - n_swap, taxo_level = taxo_level, binomial = binomial,
            species_key = species_key, species_label = species_label,
            genus = genus, famille = famille, ordre = ordre,
            taxon_unit = paste0("species:", species_key)
          )
        },
        tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = n_swap, taxo_level = taxo_level, binomial = target_label,
          species_key = target_key, species_label = target_label,
          genus = target_genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", target_key)
        )
      )
    }) %>%
      mutate(
        scenario = scenario_name,
        scenario_family = "appendix_error_gradient",
        baseline_scenario = "appendix_species_baseline",
        iter = iter_id,
        unit_type = "species_units"
      )
  }

  species_records <- col_used %>% filter(taxo_level == "species")

  scenario_baseline_appendix <- species_records %>%
    mutate(
      taxon_unit = paste0("species:", species_key),
      scenario = "appendix_species_baseline",
      scenario_family = "appendix_error_gradient",
      baseline_scenario = "appendix_species_baseline",
      iter = 1,
      unit_type = "species_units"
    ) %>%
    select(all_of(base_cols), taxon_unit, scenario, scenario_family, baseline_scenario, iter, unit_type)

  scenario_gradient_appendix <- purrr::map_dfr(seq_len(ID_ERROR_APPENDIX_N_SIM), function(iter_i) {
    set.seed(9000 + iter_i)
    purrr::pmap_dfr(
      appendix_error_grid,
      function(error_rate, error_pct, scenario) {
        simulate_same_genus_error(
          dat_species = species_records,
          pool_tbl = appendix_pool_tbl,
          scenario_name = scenario,
          iter_id = iter_i,
          fixed_rate = error_rate
        )
      }
    )
  })

  appendix_scenario_data <- bind_rows(scenario_baseline_appendix, scenario_gradient_appendix) %>%
    left_join(appendix_error_grid %>% select(scenario, error_rate, error_pct), by = "scenario") %>%
    mutate(
      error_rate = if_else(is.na(error_rate) & scenario == "appendix_species_baseline", 0, error_rate),
      error_pct = if_else(is.na(error_pct) & scenario == "appendix_species_baseline", 0, error_pct),
      pool_source = appendix_pool_name
    )

  appendix_scenario_definitions <- appendix_scenario_data %>%
    distinct(scenario, scenario_family, baseline_scenario, unit_type, error_rate, error_pct, pool_source) %>%
    arrange(error_rate, scenario)

  readr::write_csv(appendix_scenario_definitions, file.path(appendix_dir, "appendix_scenario_definitions.csv"))

  # --- Alpha/gamma
  appendix_alpha_diversity_by_unit <- appendix_scenario_data %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, pool_source) %>%
    group_split() %>%
    purrr::map_dfr(function(df) {
      mat <- make_comm_matrix(df, unit_col = ANALYSIS_UNIT)
      hill_metrics_from_matrix(mat) %>%
        mutate(
          scenario_family = df$scenario_family[1],
          scenario = df$scenario[1],
          baseline_scenario = df$baseline_scenario[1],
          unit_type = df$unit_type[1],
          iter = df$iter[1],
          error_rate = df$error_rate[1],
          error_pct = df$error_pct[1],
          pool_source = df$pool_source[1]
        )
    }) %>%
    relocate(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, pool_source, unit)

  readr::write_csv(appendix_alpha_diversity_by_unit, file.path(appendix_dir, "appendix_alpha_diversity_by_unit.csv"))

  appendix_gamma_by_scenario_iter <- appendix_scenario_data %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, pool_source) %>%
    summarise(gamma_taxon_units = n_distinct(taxon_unit), .groups = "drop")

  appendix_scenario_summary_by_iter <- appendix_alpha_diversity_by_unit %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, pool_source) %>%
    summarise(
      n_units = n(),
      total_abundance = sum(total_abundance, na.rm = TRUE),
      mean_local_q0 = mean(q0, na.rm = TRUE),
      mean_local_q1 = mean(q1, na.rm = TRUE),
      mean_local_q2 = mean(q2, na.rm = TRUE),
      mean_coverage_chao = mean(coverage_chao, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(appendix_gamma_by_scenario_iter,
              by = c("scenario_family", "scenario", "baseline_scenario", "unit_type", "iter", "error_rate", "error_pct", "pool_source"))

  readr::write_csv(appendix_scenario_summary_by_iter, file.path(appendix_dir, "appendix_scenario_summary_by_iter.csv"))

  appendix_metric_key <- tribble(
    ~metric, ~metric_label,
    "mean_local_q0", "Alpha q0 richness",
    "mean_local_q1", "Alpha Hill q1",
    "mean_local_q2", "Alpha Hill q2",
    "mean_coverage_chao", "Mean sample coverage",
    "gamma_taxon_units", "Gamma richness"
  )

  appendix_baseline_metrics <- appendix_scenario_summary_by_iter %>%
    filter(scenario == "appendix_species_baseline") %>%
    summarise(
      total_abundance = mean(total_abundance, na.rm = TRUE),
      mean_local_q0 = mean(mean_local_q0, na.rm = TRUE),
      mean_local_q1 = mean(mean_local_q1, na.rm = TRUE),
      mean_local_q2 = mean(mean_local_q2, na.rm = TRUE),
      mean_coverage_chao = mean(mean_coverage_chao, na.rm = TRUE),
      gamma_taxon_units = mean(gamma_taxon_units, na.rm = TRUE)
    ) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "baseline_value")

  appendix_biodiversity_curves <- appendix_scenario_summary_by_iter %>%
    filter(scenario != "appendix_species_baseline") %>%
    select(error_rate, error_pct, iter, total_abundance, mean_local_q0, mean_local_q1, mean_local_q2, mean_coverage_chao, gamma_taxon_units) %>%
    pivot_longer(
      cols = c(total_abundance, mean_local_q0, mean_local_q1, mean_local_q2, mean_coverage_chao, gamma_taxon_units),
      names_to = "metric",
      values_to = "value"
    ) %>%
    left_join(appendix_baseline_metrics, by = "metric") %>%
    mutate(
      response_ratio = value / baseline_value,
      percent_change = 100 * (response_ratio - 1)
    )

  appendix_biodiversity_curve_summary <- appendix_biodiversity_curves %>%
    group_by(metric, error_rate, error_pct) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      p10_value = quantile(value, 0.10, na.rm = TRUE, names = FALSE),
      p90_value = quantile(value, 0.90, na.rm = TRUE, names = FALSE),
      mean_response_ratio = mean(response_ratio, na.rm = TRUE),
      p10_response_ratio = quantile(response_ratio, 0.10, na.rm = TRUE, names = FALSE),
      p90_response_ratio = quantile(response_ratio, 0.90, na.rm = TRUE, names = FALSE),
      mean_percent_change = mean(percent_change, na.rm = TRUE),
      p10_percent_change = quantile(percent_change, 0.10, na.rm = TRUE, names = FALSE),
      p90_percent_change = quantile(percent_change, 0.90, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    left_join(appendix_metric_key, by = "metric")

  readr::write_csv(appendix_biodiversity_curve_summary, file.path(appendix_dir, "appendix_biodiversity_curve_summary.csv"))

  # --- Stabilité
  appendix_matrix_index <- appendix_scenario_data %>%
    distinct(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, pool_source) %>%
    arrange(error_rate, iter)

  make_matrix_for_key <- function(scen, iter_i) {
    appendix_scenario_data %>%
      filter(scenario == scen, iter == iter_i) %>%
      make_comm_matrix(unit_col = ANALYSIS_UNIT)
  }

  appendix_stability_by_iter <- appendix_matrix_index %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, pool_source) %>%
    group_split() %>%
    purrr::map_dfr(function(key) {
      scen <- key$scenario[1]
      iter_i <- key$iter[1]
      mat_base <- make_matrix_for_key("appendix_species_baseline", 1)
      mat_scen <- make_matrix_for_key(scen, iter_i)

      stability_one(mat_base, mat_scen) %>%
        mutate(
          scenario_family = key$scenario_family[1],
          scenario = scen,
          baseline_scenario = key$baseline_scenario[1],
          unit_type = key$unit_type[1],
          iter = iter_i,
          error_rate = key$error_rate[1],
          error_pct = key$error_pct[1],
          pool_source = key$pool_source[1]
        )
    }) %>%
    relocate(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, pool_source)

  readr::write_csv(appendix_stability_by_iter, file.path(appendix_dir, "appendix_stability_by_iter.csv"))

  appendix_stability_long <- appendix_stability_by_iter %>%
    filter(scenario != "appendix_species_baseline") %>%
    select(error_rate, error_pct, iter,
           q0_cor, q1_cor, q2_cor, coverage_chao_cor,
           bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2) %>%
    pivot_longer(
      cols = c(q0_cor, q1_cor, q2_cor, coverage_chao_cor,
               bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2),
      names_to = "metric",
      values_to = "stability"
    ) %>%
    mutate(
      metric_label = recode(
        metric,
        q0_cor = "Alpha q0 richness",
        q1_cor = "Alpha Hill q1",
        q2_cor = "Alpha Hill q2",
        coverage_chao_cor = "Mean sample coverage",
        bray_curtis_cor = "Bray-Curtis",
        jaccard_pa_cor = "Jaccard",
        sorensen_pa_cor = "Sørensen",
        procrustes_r2 = "Ordination (Procrustes R2)"
      )
    )

  appendix_stability_curve_summary <- appendix_stability_long %>%
    group_by(metric, metric_label, error_rate, error_pct) %>%
    summarise(
      mean_stability = mean(stability, na.rm = TRUE),
      p10_stability = quantile(stability, 0.10, na.rm = TRUE, names = FALSE),
      p90_stability = quantile(stability, 0.90, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )

  readr::write_csv(appendix_stability_curve_summary, file.path(appendix_dir, "appendix_stability_curve_summary.csv"))

  appendix_threshold_biodiv <- appendix_biodiversity_curve_summary %>%
    filter(metric %in% appendix_metric_key$metric) %>%
    group_by(metric, metric_label) %>%
    summarise(
      critical_definition = paste0("|relative change| >= ", round(100 * ID_ERROR_APPENDIX_CRITICAL_REL_CHANGE), "%"),
      critical_error_pct = suppressWarnings(min(error_pct[abs(mean_percent_change) >= (100 * ID_ERROR_APPENDIX_CRITICAL_REL_CHANGE)], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(critical_error_pct = ifelse(is.infinite(critical_error_pct), NA_real_, critical_error_pct))

  appendix_threshold_stability <- appendix_stability_curve_summary %>%
    group_by(metric, metric_label) %>%
    summarise(
      critical_definition = paste0("mean stability <= ", ID_ERROR_APPENDIX_CRITICAL_STABILITY),
      critical_error_pct = suppressWarnings(min(error_pct[mean_stability <= ID_ERROR_APPENDIX_CRITICAL_STABILITY], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(critical_error_pct = ifelse(is.infinite(critical_error_pct), NA_real_, critical_error_pct))

  appendix_critical_thresholds <- bind_rows(
    appendix_threshold_biodiv %>% mutate(threshold_type = "relative_change"),
    appendix_threshold_stability %>% mutate(threshold_type = "stability")
  ) %>%
    relocate(threshold_type, metric, metric_label, critical_definition, critical_error_pct)

  readr::write_csv(appendix_critical_thresholds, file.path(appendix_dir, "appendix_critical_thresholds.csv"))

  save_step(STEP_ID, c(
    "appendix_dir",
    "appendix_pool_name",
    "appendix_pool_tbl",
    "appendix_error_grid",
    "appendix_scenario_data",
    "appendix_scenario_definitions",
    "appendix_alpha_diversity_by_unit",
    "appendix_gamma_by_scenario_iter",
    "appendix_scenario_summary_by_iter",
    "appendix_biodiversity_curves",
    "appendix_biodiversity_curve_summary",
    "appendix_matrix_index",
    "appendix_stability_by_iter",
    "appendix_stability_long",
    "appendix_stability_curve_summary",
    "appendix_critical_thresholds"
  ))
}
