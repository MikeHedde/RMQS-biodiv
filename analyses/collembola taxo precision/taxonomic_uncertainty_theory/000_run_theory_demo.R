# ============================================================
# 000_run_theory_demo.R
# Demonstration run of the stand-alone theoretical model
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
THEORY_DIR <- "analyses/collembola taxo precision/taxonomic_uncertainty_theory"

source(file.path(THEORY_DIR, "00_theory_config.R"))
source(file.path(THEORY_DIR, "01_taxonomic_architecture.R"))
source(file.path(THEORY_DIR, "02_community_generator.R"))
source(file.path(THEORY_DIR, "03_observation_process.R"))
source(file.path(THEORY_DIR, "04_metrics.R"))
source(file.path(THEORY_DIR, "05_simulation_engine.R"))
source(file.path(THEORY_DIR, "06_experiment_design.R"))
source(file.path(THEORY_DIR, "07_empirical_bridge.R"))
source(file.path(THEORY_DIR, "08_plotting.R"))

# This demo explores theoretical worlds. It does not use Collembola,
# earthworm or ant data, and does not calibrate any parameter on them.
design <- make_world_design(
  n_worlds = 60,
  seed = THEORY_SEED
)

scenarios <- default_scenario_library(
  id_error_rate = 0.10,
  workflow_unresolved_rate = 0.20,
  rare_quantile = 0.10
)

experiment <- run_theoretical_experiment(
  world_design = design,
  scenario_library = scenarios,
  n_replicates_per_world = 5,
  seed = THEORY_SEED
)

theory_summary <- summarise_theoretical_experiment(experiment)

saveRDS(experiment, file.path(THEORY_OUT_DIR, "theory_experiment.rds"))
write_theory_csv(experiment$world_profiles, "theory_world_profiles.csv")
write_theory_csv(experiment$blowes, "theory_blowes_by_replicate.csv")
write_theory_csv(experiment$effects, "theory_effects_by_replicate.csv")
write_theory_csv(experiment$stability, "theory_stability_by_replicate.csv")
write_theory_csv(experiment$process, "theory_process_by_replicate.csv")
write_theory_csv(theory_summary, "theory_blowes_summary.csv")

p_family <- plot_theory_blowes(theory_summary, scenario = "family_level")
ggsave(
  filename = file.path(THEORY_OUT_DIR, "theory_blowes_family_level.png"),
  plot = p_family,
  width = 7.5,
  height = 6,
  dpi = 300
)

p_regional <- plot_theory_blowes(theory_summary, scenario = "regional_congeneric_error")
ggsave(
  filename = file.path(THEORY_OUT_DIR, "theory_blowes_regional_error.png"),
  plot = p_regional,
  width = 7.5,
  height = 6,
  dpi = 300
)

message("Theoretical simulation finished. Outputs: ", THEORY_OUT_DIR)
