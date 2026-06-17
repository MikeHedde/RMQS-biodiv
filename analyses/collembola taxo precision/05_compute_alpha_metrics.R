# ============================================================
# 05_compute_alpha_metrics — alpha diversité, couverture et gamma
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances
load_required_step("04_build_scenarios")

STEP_ID <- "05_compute_alpha_metrics"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("05_compute_alpha_metrics — alpha diversité, couverture et gamma : cache chargé")
} else {
  # -----------------------------
  # 6. METRIQUES ALPHA + COUVERTURE
  # -----------------------------

  message_header("Calcul diversité alpha et couverture avec divent")

  alpha_diversity_by_unit <- scenario_data %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
    group_split() %>%
    purrr::map_dfr(function(df) {
      mat <- make_comm_matrix(df, unit_col = ANALYSIS_UNIT)
      hill_metrics_from_matrix(mat) %>%
        mutate(
          scenario_family = df$scenario_family[1],
          scenario = df$scenario[1],
          baseline_scenario = df$baseline_scenario[1],
          unit_type = df$unit_type[1],
          iter = df$iter[1]
        )
    }) %>%
    relocate(scenario_family, scenario, baseline_scenario, unit_type, iter, unit)

  readr::write_csv(alpha_diversity_by_unit, file.path(OUT_DIR, "alpha_diversity_and_coverage_by_unit.csv"))

  gamma_by_scenario_iter <- scenario_data %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
    summarise(gamma_taxon_units = n_distinct(taxon_unit), .groups = "drop")

  scenario_summary_by_iter <- alpha_diversity_by_unit %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
    summarise(
      n_units = n(),
      total_abundance = sum(total_abundance, na.rm = TRUE),
      mean_local_q0 = mean(q0, na.rm = TRUE),
      mean_local_q1 = mean(q1, na.rm = TRUE),
      mean_local_q2 = mean(q2, na.rm = TRUE),
      mean_coverage_chao = mean(coverage_chao, na.rm = TRUE),
      p10_coverage_chao = quantile(coverage_chao, probs = 0.10, na.rm = TRUE, names = FALSE),
      median_coverage_chao = median(coverage_chao, na.rm = TRUE),
      p90_coverage_chao = quantile(coverage_chao, probs = 0.90, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    left_join(gamma_by_scenario_iter, by = c("scenario_family", "scenario", "baseline_scenario", "unit_type", "iter"))

  readr::write_csv(scenario_summary_by_iter, file.path(OUT_DIR, "scenario_summary_by_iter.csv"))

  coverage_summary_by_iter <- alpha_diversity_by_unit %>%
    select(scenario_family, scenario, baseline_scenario, unit_type, iter, unit, starts_with("coverage_")) %>%
    pivot_longer(starts_with("coverage_"), names_to = "coverage_estimator", values_to = "coverage") %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, coverage_estimator) %>%
    summarise(
      mean = mean(coverage, na.rm = TRUE),
      p10 = quantile(coverage, 0.10, na.rm = TRUE, names = FALSE),
      median = median(coverage, na.rm = TRUE),
      p90 = quantile(coverage, 0.90, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
  readr::write_csv(coverage_summary_by_iter, file.path(OUT_DIR, "coverage_summary_by_iter.csv"))

  # -----------------------------

  save_step(STEP_ID, c(
    "alpha_diversity_by_unit",
    "gamma_by_scenario_iter",
    "scenario_summary_by_iter",
    "coverage_summary_by_iter"
  ))
}
