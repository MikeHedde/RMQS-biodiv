# ============================================================
# 003_run_taxon_profiles.R
# Build empirical profiles, then run cross-taxon theory simulations
# ============================================================
#
# Stage A: reads standardised inputs and creates descriptive profiles.
# Stage B: optionally simulates generic communities using each group's real
# regional taxonomy. Empirical scenario outcomes are never used as targets.

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
THEORY_DIR <- if (!is.na(.workflow_file)) dirname(.workflow_file) else getwd()

source(file.path(THEORY_DIR, "00_theory_config.R"))
source(file.path(THEORY_DIR, "01_taxonomic_architecture.R"))
source(file.path(THEORY_DIR, "02_community_generator.R"))
source(file.path(THEORY_DIR, "03_observation_process.R"))
source(file.path(THEORY_DIR, "04_metrics.R"))
source(file.path(THEORY_DIR, "05_simulation_engine.R"))
source(file.path(THEORY_DIR, "10_controlled_one_factor_experiment.R"))
source(file.path(THEORY_DIR, "12_taxon_profile_builder.R"))
source(file.path(THEORY_DIR, "13_taxon_specific_simulation.R"))
source(file.path(THEORY_DIR, "14_taxon_profile_plots.R"))

INPUT_CONFIG <- file.path(THEORY_DIR, "inputs", "taxon_profile_inputs.R")
if (!file.exists(INPUT_CONFIG)) {
  stop("Missing input configuration: ", INPUT_CONFIG, call. = FALSE)
}
source(INPUT_CONFIG)

# Resolve relative paths from THEORY_DIR rather than from the R working directory.
resolve_theory_input <- function(x) {
  if (grepl("^([A-Za-z]:)?[\\\\/]", x)) x else file.path(THEORY_DIR, x)
}

# ----------------------------------------------------------------
# Stage A. Build profiles
# ----------------------------------------------------------------
profiles <- purrr::imap(TAXON_PROFILE_INPUTS, function(spec, taxon) {
  community_file <- resolve_theory_input(spec$community_file)
  taxonomy_file <- resolve_theory_input(spec$taxonomy_file)

  if (!file.exists(community_file) || !file.exists(taxonomy_file)) {
    stop(
      sprintf("Missing standard inputs for %s: %s / %s", taxon, community_file, taxonomy_file),
      call. = FALSE
    )
  }

  inputs <- read_standard_taxon_inputs(
    community_file = community_file,
    taxonomy_file = taxonomy_file,
    community_delim = spec$community_delim %||% ",",
    taxonomy_delim = spec$taxonomy_delim %||% ","
  )

  build_empirical_taxon_profile(
    taxon = taxon,
    community_long = inputs$community_long,
    regional_taxonomy = inputs$regional_taxonomy
  )
})

purrr::walk(profiles, write_empirical_taxon_profile)
profile_summaries <- bind_rows(purrr::map(profiles, "summary"))
write_theory_csv(profile_summaries, "empirical_taxon_profiles.csv")

fig_profiles <- plot_taxon_profile_axes(profile_summaries)
ggsave(file.path(THEORY_OUT_DIR, "FigP1_empirical_taxon_profiles.pdf"), fig_profiles, width = 10, height = 4.5, units = "in")
ggsave(file.path(THEORY_OUT_DIR, "FigP1_empirical_taxon_profiles.png"), fig_profiles, width = 10, height = 4.5, units = "in", dpi = 300)

# ----------------------------------------------------------------
# Stage B. Generic taxon-specific simulations
# ----------------------------------------------------------------
# Set FALSE for profile-only runs. Use 10 initially; increase to >= 50 for
# stable cross-taxon uncertainty envelopes once all groups are available.
RUN_TAXON_SPECIFIC_SIMULATIONS <- TRUE
N_REPLICATES_PER_TAXON <- 10

if (RUN_TAXON_SPECIFIC_SIMULATIONS) {
  scenario_grid <- make_controlled_scenario_grid(
    id_error_rates = c(0.01, 0.03, 0.05, 0.10, 0.20),
    reporting_rates = c(0.01, 0.03, 0.05, 0.10, 0.20),
    rare_quantile = 0.10,
    difficult_multiplier = 4
  )

  experiment <- run_taxon_profile_experiment(
    profiles = profiles,
    scenario_grid = scenario_grid,
    n_replicates_per_taxon = N_REPLICATES_PER_TAXON,
    seed = THEORY_SEED,
    progress = TRUE
  )
  summary_out <- summarise_taxon_profile_experiment(experiment)

  saveRDS(experiment, file.path(THEORY_OUT_DIR, "taxon_profile_experiment.rds"))
  write_theory_csv(experiment$profiles, "taxon_profile_simulated_worlds.csv")
  write_theory_csv(experiment$results, "taxon_profile_simulation_by_replicate.csv")
  write_theory_csv(summary_out$blowes_summary, "taxon_profile_blowes_summary.csv")
  write_theory_csv(summary_out$stability_summary, "taxon_profile_stability_summary.csv")

  fig_gamma <- plot_cross_taxon_scenarios(summary_out$blowes_summary, response = "delta_gamma")
  fig_occupancy <- plot_cross_taxon_scenarios(summary_out$blowes_summary, response = "delta_occupancy")
  fig_blowes <- plot_cross_taxon_blowes(summary_out$blowes_summary)

  ggsave(file.path(THEORY_OUT_DIR, "FigP2_cross_taxon_delta_gamma.pdf"), fig_gamma, width = 10, height = 6, units = "in")
  ggsave(file.path(THEORY_OUT_DIR, "FigP3_cross_taxon_delta_occupancy.pdf"), fig_occupancy, width = 10, height = 6, units = "in")
  ggsave(file.path(THEORY_OUT_DIR, "FigP4_cross_taxon_blowes.pdf"), fig_blowes, width = 9, height = 7, units = "in")

  ggsave(file.path(THEORY_OUT_DIR, "FigP2_cross_taxon_delta_gamma.png"), fig_gamma, width = 10, height = 6, units = "in", dpi = 300)
  ggsave(file.path(THEORY_OUT_DIR, "FigP3_cross_taxon_delta_occupancy.png"), fig_occupancy, width = 10, height = 6, units = "in", dpi = 300)
  ggsave(file.path(THEORY_OUT_DIR, "FigP4_cross_taxon_blowes.png"), fig_blowes, width = 9, height = 7, units = "in", dpi = 300)
}

message("Taxon-profile stage completed. Outputs: ", THEORY_OUT_DIR)
