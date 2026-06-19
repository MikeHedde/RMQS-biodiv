# ============================================================
# 002_run_controlled_one_factor.R
# Controlled one-factor-at-a-time simulation experiment
# ============================================================
#
# This is the next analytical stage after the exploratory random-world screen.
# It varies one mechanism at a time while holding the other three at reference
# values, allowing direct and interpretable general predictions.

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)

source(file.path(THEORY_DIR, "00_theory_config.R"))
source(file.path(THEORY_DIR, "01_taxonomic_architecture.R"))
source(file.path(THEORY_DIR, "02_community_generator.R"))
source(file.path(THEORY_DIR, "03_observation_process.R"))
source(file.path(THEORY_DIR, "04_metrics.R"))
source(file.path(THEORY_DIR, "05_simulation_engine.R"))
source(file.path(THEORY_DIR, "10_controlled_one_factor_experiment.R"))
source(file.path(THEORY_DIR, "11_controlled_one_factor_plots.R"))

# ----------------------------------------------------------------
# User-facing configuration
# ----------------------------------------------------------------
# First full test: 5 replicates. For final production, use 20 or more.
N_REPLICATES_PER_WORLD <- 5

ID_ERROR_RATES <- c(0.01, 0.03, 0.05, 0.10, 0.20)
REPORTING_RATES <- c(0.01, 0.03, 0.05, 0.10, 0.20)

reference <- controlled_reference_parameters()
levels <- controlled_factor_levels()

design <- make_controlled_one_factor_design(
  reference = reference,
  factor_levels = levels
)

scenario_grid <- make_controlled_scenario_grid(
  id_error_rates = ID_ERROR_RATES,
  reporting_rates = REPORTING_RATES,
  rare_quantile = 0.10,
  difficult_multiplier = 4
)

experiment <- run_controlled_one_factor_experiment(
  design = design,
  scenario_grid = scenario_grid,
  n_replicates_per_world = N_REPLICATES_PER_WORLD,
  seed = THEORY_SEED,
  progress = TRUE
)

summary_out <- summarise_controlled_one_factor_experiment(experiment)

# ----------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------
saveRDS(experiment, file.path(THEORY_OUT_DIR, "controlled_ofat_experiment.rds"))
write_theory_csv(design, "controlled_ofat_design.csv")
write_theory_csv(scenario_grid %>% select(-baseline_process, -target_process), "controlled_ofat_scenario_grid.csv")
write_theory_csv(experiment$world_profiles, "controlled_ofat_world_profiles.csv")
write_theory_csv(experiment$blowes, "controlled_ofat_blowes_by_replicate.csv")
write_theory_csv(experiment$effects, "controlled_ofat_effects_by_replicate.csv")
write_theory_csv(experiment$stability, "controlled_ofat_stability_by_replicate.csv")
write_theory_csv(experiment$process, "controlled_ofat_process_by_replicate.csv")
write_theory_csv(summary_out$blowes_summary, "controlled_ofat_blowes_summary.csv")
write_theory_csv(summary_out$effects_summary, "controlled_ofat_effects_summary.csv")
write_theory_csv(summary_out$stability_summary, "controlled_ofat_stability_summary.csv")
write_theory_csv(summary_out$process_summary, "controlled_ofat_process_summary.csv")
write_theory_csv(summary_out$endpoint_contrasts, "controlled_ofat_endpoint_contrasts.csv")

# Identification-error mechanisms.
fig_id_gamma <- controlled_response_plot(
  summary_out$blowes_summary,
  scenario_domain = "identification_error",
  response = "delta_gamma",
  title = "Controlled effects of identification-error mechanisms on regional richness"
)
fig_id_occ <- controlled_response_plot(
  summary_out$blowes_summary,
  scenario_domain = "identification_error",
  response = "delta_occupancy",
  title = "Controlled effects of identification-error mechanisms on mean occupancy"
)

# Reporting workflows with a variable unresolved-record rate.
fig_reporting_gamma <- controlled_response_plot(
  summary_out$blowes_summary,
  scenario_domain = "reporting_workflow",
  response = "delta_gamma",
  title = "Controlled effects of reporting workflows on regional richness"
)
fig_reporting_occ <- controlled_response_plot(
  summary_out$blowes_summary,
  scenario_domain = "reporting_workflow",
  response = "delta_occupancy",
  title = "Controlled effects of reporting workflows on mean occupancy"
)

# Fixed coarsening/reporting workflows.
fig_fixed_gamma <- controlled_fixed_workflow_plot(
  summary_out$blowes_summary,
  response = "delta_gamma",
  title = "Controlled effects of fixed taxonomic workflows on regional richness"
)
fig_fixed_occ <- controlled_fixed_workflow_plot(
  summary_out$blowes_summary,
  response = "delta_occupancy",
  title = "Controlled effects of fixed taxonomic workflows on mean occupancy"
)

fig_endpoint_gamma <- controlled_endpoint_heatmap(
  summary_out$endpoint_contrasts,
  response = "delta_gamma_med"
)
fig_endpoint_occ <- controlled_endpoint_heatmap(
  summary_out$endpoint_contrasts,
  response = "delta_occupancy_med"
)

figures <- list(
  FigF1_controlled_ID_delta_gamma = fig_id_gamma,
  FigF2_controlled_ID_delta_occupancy = fig_id_occ,
  FigF3_controlled_reporting_delta_gamma = fig_reporting_gamma,
  FigF4_controlled_reporting_delta_occupancy = fig_reporting_occ,
  FigF5_controlled_fixed_workflows_delta_gamma = fig_fixed_gamma,
  FigF6_controlled_fixed_workflows_delta_occupancy = fig_fixed_occ,
  FigF7_controlled_endpoint_gamma = fig_endpoint_gamma,
  FigF8_controlled_endpoint_occupancy = fig_endpoint_occ
)

purrr::iwalk(figures, function(fig, name) {
  ggsave(
    filename = file.path(THEORY_OUT_DIR, paste0(name, ".pdf")),
    plot = fig,
    width = 11,
    height = 8,
    units = "in"
  )
  ggsave(
    filename = file.path(THEORY_OUT_DIR, paste0(name, ".png")),
    plot = fig,
    width = 11,
    height = 8,
    units = "in",
    dpi = 300
  )
})

message("Controlled OFAT experiment completed. Outputs: ", THEORY_OUT_DIR)
