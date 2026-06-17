# ============================================================
# 06_compute_stability — stabilité vs baseline
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances
load_required_step("04_build_scenarios")

STEP_ID <- "06_compute_stability"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("06_compute_stability — stabilité vs baseline : cache chargé")
} else {
  # 7. STABILITE VS BASELINE PROPRE A CHAQUE FAMILLE DE SCENARIOS
  # -----------------------------

  message_header("Calcul stabilité des inférences")

  # Matrices par scénario/itération
  matrix_index <- scenario_data %>%
    distinct(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
    arrange(scenario_family, scenario, iter)

  make_matrix_for_key <- function(fam, scen, iter_i) {
    scenario_data %>%
      filter(scenario_family == fam, scenario == scen, iter == iter_i) %>%
      make_comm_matrix(unit_col = ANALYSIS_UNIT)
  }

  stability_by_iter <- matrix_index %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
    group_split() %>%
    purrr::map_dfr(function(key) {
      fam <- key$scenario_family[1]
      scen <- key$scenario[1]
      base <- key$baseline_scenario[1]
      iter_i <- key$iter[1]
    
      # Baseline déterministe : iter = 1.
      mat_base <- make_matrix_for_key(fam, base, 1)
      mat_scen <- make_matrix_for_key(fam, scen, iter_i)
    
      stability_one(mat_base, mat_scen) %>%
        mutate(
          scenario_family = fam,
          scenario = scen,
          baseline_scenario = base,
          unit_type = key$unit_type[1],
          iter = iter_i
        )
    }) %>%
    relocate(scenario_family, scenario, baseline_scenario, unit_type, iter)

  readr::write_csv(stability_by_iter, file.path(OUT_DIR, "stability_by_iter.csv"))

  stability_summary <- stability_by_iter %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type) %>%
    summarise(
      across(c(abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
               bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor,
               bray_mean_ratio, jaccard_mean_ratio, sorensen_mean_ratio,
               procrustes_r2, n_common_units),
             safe_mean, .names = "{.col}_mean"),
      across(c(abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
               bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor,
               bray_mean_ratio, jaccard_mean_ratio, sorensen_mean_ratio,
               procrustes_r2, n_common_units),
             safe_sd, .names = "{.col}_sd"),
      n_iter = n(),
      .groups = "drop"
    )
  readr::write_csv(stability_summary, file.path(OUT_DIR, "stability_summary.csv"))

  # Format long pour figures
  stability_long <- stability_by_iter %>%
    select(scenario_family, scenario, baseline_scenario, unit_type, iter,
           abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
           bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2) %>%
    pivot_longer(
      cols = c(abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
               bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2),
      names_to = "inference",
      values_to = "stability"
    ) %>%
    mutate(
      inference_label = recode(
        inference,
        abundance_cor = "total abundance",
        q0_cor = "q0 richness / number of units",
        q1_cor = "Hill q1 / effective units",
        q2_cor = "Hill q2 / dominant units",
        coverage_chao_cor = "sample coverage Chao",
        bray_curtis_cor = "Bray-Curtis abundance distance",
        jaccard_pa_cor = "Jaccard presence-absence distance",
        sorensen_pa_cor = "Sørensen presence-absence distance",
        procrustes_r2 = "ordination Procrustes R2"
      )
    )
  readr::write_csv(stability_long, file.path(OUT_DIR, "stability_long.csv"))

  # -----------------------------

  save_step(STEP_ID, c(
    "matrix_index",
    "stability_by_iter",
    "stability_summary",
    "stability_long"
  ))
}
