# ============================================================
# 09_explore_parameter_sensitivity.R
# Exploratory explanation of the theoretical alpha-gamma space
# ============================================================
#
# Purpose
# -------
# This script does NOT change the simulator and does NOT calibrate it
# on empirical taxa. It reads a completed theoretical experiment and asks:
#
#   Which properties of the regional taxonomic architecture and of the
#   true community covary with the direction and magnitude of the
#   simulated artefact?
#
# The present script is deliberately exploratory. It uses rank
# associations across independently generated theoretical worlds, rather
# than treating the simulated worlds as an inferential sample from nature.
#
# Use it after 000_run_theory_demo.R.
#
# Outputs
# -------
# - theory_parameter_space.csv
# - theory_parameter_associations.csv
# - theory_parameter_associations_top.csv
# - FigT3_regional_error_parameter_screen.png / .pdf
# - FigT4_family_coarsening_parameter_screen.png / .pdf
# - FigT5_regional_error_response_curves.png / .pdf
# - FigT6_family_coarsening_response_curves.png / .pdf
#
# Interpretation
# --------------
# The results identify candidate mechanisms for a subsequent, structured
# factorial experiment. They should not yet be interpreted as unique or
# causal effects because several simulated descriptors remain correlated
# by construction (e.g. regional richness, numbers of genera and families).

.theory_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)


source(file.path(THEORY_DIR, "00_theory_config.R"))
source(file.path(THEORY_DIR, "01_taxonomic_architecture.R"))
source(file.path(THEORY_DIR, "02_community_generator.R"))
source(file.path(THEORY_DIR, "03_observation_process.R"))
source(file.path(THEORY_DIR, "04_metrics.R"))
source(file.path(THEORY_DIR, "05_simulation_engine.R"))
source(file.path(THEORY_DIR, "08_plotting.R"))

EXPERIMENT_FILE <- file.path(THEORY_OUT_DIR, "theory_experiment.rds")
if (!file.exists(EXPERIMENT_FILE)) {
  stop(
    "No theoretical experiment found. Run `000_run_theory_demo.R` first, ",
    "or set `taxo.theory_dir` to the directory containing outputs/theory_experiment.rds.",
    call. = FALSE
  )
}

experiment <- readRDS(EXPERIMENT_FILE)
theory_summary <- summarise_theoretical_experiment(experiment)

# -------------------------------------------------------------------------
# 1. Build a transparent, derived descriptor table for each theoretical world
# -------------------------------------------------------------------------

# These design parameters are not necessarily retained in world_profiles.
# We join only variables that are not duplicated by the profile table.
design_keep <- experiment$design %>%
  distinct(world_id, .keep_all = TRUE) %>%
  select(
    world_id,
    n_sites,
    species_per_genus_concentration,
    genera_per_family_concentration,
    difficult_genus_fraction,
    occupancy_mean,
    occupancy_precision,
    mean_abundance,
    abundance_nb_size,
    occupancy_abundance_slope,
    observed_selection_bias,
    environmental_gradient_strength
  )

# Process-realisation descriptors are useful because the number of eligible
# species differs among worlds even when the requested error rate is fixed.
process_summary <- experiment$process %>%
  group_by(world_id, scenario) %>%
  summarise(
    realised_mean_error = safe_mean(realised_mean_error),
    realised_unresolved_rate = safe_mean(realised_unresolved_rate),
    n_eligible_source_species = safe_mean(n_eligible_source_species),
    .groups = "drop"
  )

parameter_space <- theory_summary %>%
  left_join(design_keep, by = "world_id", suffix = c("", "_design")) %>%
  left_join(process_summary, by = c("world_id", "scenario")) %>%
  mutate(
    # Pool size available for new labels relative to the sampled species pool.
    regional_to_observed_ratio = n_regional_species / pmax(n_observed_species, 1),
    observed_fraction = n_observed_species / pmax(n_regional_species, 1),

    # Compression potential when all species are merged into a higher rank.
    species_per_family = n_regional_species / pmax(n_families, 1),
    taxonomic_compression_genus = n_regional_species / pmax(n_genera, 1),
    taxonomic_compression_family = n_regional_species / pmax(n_families, 1),

    # Log transforms are supplied for modelling/screening, while figures retain
    # the original values whenever possible.
    log_regional_to_observed = log10(regional_to_observed_ratio),
    log_mean_species_per_genus = log10(mean_species_per_genus),
    log_mean_genera_per_family = log10(mean_genera_per_family),
    log_species_per_family = log10(species_per_family),
    log_n_sites = log10(n_sites)
  )

write_theory_csv(parameter_space, "theory_parameter_space.csv")

# -------------------------------------------------------------------------
# 2. Rank-association screen
# -------------------------------------------------------------------------

screen_predictors <- c(
  "log_regional_to_observed",
  "log_mean_species_per_genus",
  "cv_species_per_genus",
  "prop_monotypic_genera",
  "log_mean_genera_per_family",
  "cv_genera_per_family",
  "species_per_genus_concentration",
  "genera_per_family_concentration",
  "taxonomic_compression_family",
  "occupancy_mean",
  "occupancy_precision",
  "true_mean_occupancy",
  "true_alpha_q0",
  "true_gamma",
  "log_n_sites",
  "realised_mean_error",
  "n_eligible_source_species"
)

predictor_labels <- c(
  log_regional_to_observed = "Regional / observed species ratio (log10)",
  log_mean_species_per_genus = "Mean species per genus (log10)",
  cv_species_per_genus = "CV of species per genus",
  prop_monotypic_genera = "Proportion of monotypic genera",
  log_mean_genera_per_family = "Mean genera per family (log10)",
  cv_genera_per_family = "CV of genera per family",
  species_per_genus_concentration = "Species-per-genus concentration parameter",
  genera_per_family_concentration = "Genera-per-family concentration parameter",
  taxonomic_compression_family = "Species per family",
  occupancy_mean = "Input mean occupancy",
  occupancy_precision = "Occupancy heterogeneity parameter",
  true_mean_occupancy = "Realised mean occupancy",
  true_alpha_q0 = "True mean local richness",
  true_gamma = "True gamma richness",
  log_n_sites = "Number of sites (log10)",
  realised_mean_error = "Realised identification-error rate",
  n_eligible_source_species = "Eligible source species"
)

screen_associations <- function(data, scenario_name, outcome) {
  dat <- data %>% filter(scenario == scenario_name)

  purrr::map_dfr(screen_predictors, function(predictor) {
    if (!predictor %in% names(dat) || !outcome %in% names(dat)) {
      return(tibble(
        scenario = scenario_name,
        outcome = outcome,
        predictor = predictor,
        n = NA_integer_,
        spearman_rho = NA_real_,
        spearman_p = NA_real_
      ))
    }

    dd <- dat %>% select(all_of(c(predictor, outcome))) %>% tidyr::drop_na()
    if (nrow(dd) < 5 || dplyr::n_distinct(dd[[predictor]]) < 3 || dplyr::n_distinct(dd[[outcome]]) < 3) {
      return(tibble(
        scenario = scenario_name,
        outcome = outcome,
        predictor = predictor,
        n = nrow(dd),
        spearman_rho = NA_real_,
        spearman_p = NA_real_
      ))
    }

    test <- suppressWarnings(stats::cor.test(
      dd[[predictor]], dd[[outcome]], method = "spearman", exact = FALSE
    ))

    tibble(
      scenario = scenario_name,
      outcome = outcome,
      predictor = predictor,
      n = nrow(dd),
      spearman_rho = unname(test$estimate),
      spearman_p = test$p.value
    )
  }) %>%
    mutate(
      predictor_label = unname(predictor_labels[predictor]),
      abs_rho = abs(spearman_rho)
    ) %>%
    arrange(desc(abs_rho))
}

# Scenario-specific outcomes reflect the mechanisms illustrated by the
# two first theoretical figures.
association_targets <- tribble(
  ~scenario, ~outcome, ~screen_label,
  "regional_congeneric_error", "delta_gamma_med", "Regional-pool error: change in gamma richness",
  "regional_congeneric_error", "delta_occupancy_med", "Regional-pool error: change in mean occupancy",
  "family_level", "delta_occupancy_med", "Family coarsening: change in mean occupancy",
  "family_level", "delta_gamma_med", "Family coarsening: change in gamma richness"
)

parameter_associations <- association_targets %>%
  pmap_dfr(function(scenario, outcome, screen_label) {
    screen_associations(parameter_space, scenario_name = scenario, outcome = outcome) %>%
      mutate(screen_label = screen_label)
  })

write_theory_csv(parameter_associations, "theory_parameter_associations.csv")

# Top associations are simply descriptive candidates for the subsequent
# designed factorial experiment; they are not p-value-selected conclusions.
parameter_associations_top <- parameter_associations %>%
  group_by(scenario, outcome) %>%
  slice_max(abs_rho, n = 6, with_ties = FALSE) %>%
  ungroup()

write_theory_csv(parameter_associations_top, "theory_parameter_associations_top.csv")

# -------------------------------------------------------------------------
# 3. Screening figures: rank associations
# -------------------------------------------------------------------------

plot_rank_screen <- function(association_table, scenario_name, outcomes, title) {
  dat <- association_table %>%
    filter(scenario == scenario_name, outcome %in% outcomes) %>%
    group_by(outcome) %>%
    slice_max(abs_rho, n = 8, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      predictor_label = factor(predictor_label, levels = rev(unique(predictor_label))),
      outcome = recode(
        outcome,
        delta_gamma_med = "Change in gamma richness",
        delta_occupancy_med = "Change in mean occupancy"
      )
    )

  ggplot(dat, aes(x = predictor_label, y = spearman_rho)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_col(width = 0.72) +
    coord_flip() +
    facet_wrap(~ outcome, scales = "free_y") +
    labs(
      title = title,
      subtitle = "Exploratory Spearman associations across theoretical worlds; used to identify candidate mechanisms",
      x = NULL,
      y = "Spearman rho"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )
}

p_regional_screen <- plot_rank_screen(
  association_table = parameter_associations,
  scenario_name = "regional_congeneric_error",
  outcomes = c("delta_gamma_med", "delta_occupancy_med"),
  title = "Theoretical drivers of regional-pool congeneric-error artefacts"
)

p_family_screen <- plot_rank_screen(
  association_table = parameter_associations,
  scenario_name = "family_level",
  outcomes = c("delta_gamma_med", "delta_occupancy_med"),
  title = "Theoretical drivers of family-level coarsening artefacts"
)

ggsave(
  file.path(THEORY_OUT_DIR, "FigT3_regional_error_parameter_screen.png"),
  p_regional_screen, width = 9.2, height = 6.4, dpi = 300
)
ggsave(
  file.path(THEORY_OUT_DIR, "FigT3_regional_error_parameter_screen.pdf"),
  p_regional_screen, width = 9.2, height = 6.4
)
ggsave(
  file.path(THEORY_OUT_DIR, "FigT4_family_coarsening_parameter_screen.png"),
  p_family_screen, width = 9.2, height = 6.4, dpi = 300
)
ggsave(
  file.path(THEORY_OUT_DIR, "FigT4_family_coarsening_parameter_screen.pdf"),
  p_family_screen, width = 9.2, height = 6.4
)

# -------------------------------------------------------------------------
# 4. Response curves for a pre-specified set of interpretable mechanisms
# -------------------------------------------------------------------------

make_response_data <- function(data, scenario_name, outcome, predictors) {
  data %>%
    filter(scenario == scenario_name) %>%
    select(all_of(c("world_id", outcome, predictors))) %>%
    pivot_longer(
      cols = all_of(predictors),
      names_to = "predictor",
      values_to = "x"
    ) %>%
    mutate(
      predictor_label = unname(predictor_labels[predictor])
    )
}

plot_response_curves <- function(data, scenario_name, outcome, predictors, title) {
  dat <- make_response_data(data, scenario_name, outcome, predictors) %>%
    filter(is.finite(x), is.finite(.data[[outcome]]))

  ggplot(dat, aes(x = x, y = .data[[outcome]])) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(alpha = 0.45, size = 1.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.7) +
    facet_wrap(~ predictor_label, scales = "free_x", ncol = 2) +
    labs(
      title = title,
      subtitle = "Each panel varies across independently generated theoretical worlds; curves are descriptive",
      x = NULL,
      y = if (outcome == "delta_gamma_med") expression(Delta * gamma ~ "(%)") else expression(Delta * " mean occupancy (%)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )
}

regional_predictors <- c(
  "regional_to_observed_ratio",
  "mean_species_per_genus",
  "cv_species_per_genus",
  "true_mean_occupancy"
)

family_predictors <- c(
  "species_per_family",
  "mean_species_per_genus",
  "mean_genera_per_family",
  "prop_monotypic_genera"
)

p_regional_curves <- plot_response_curves(
  parameter_space,
  scenario_name = "regional_congeneric_error",
  outcome = "delta_gamma_med",
  predictors = regional_predictors,
  title = "Mechanistic candidates for gamma inflation under regional-pool error"
)

p_family_curves <- plot_response_curves(
  parameter_space,
  scenario_name = "family_level",
  outcome = "delta_occupancy_med",
  predictors = family_predictors,
  title = "Mechanistic candidates for artificial homogenisation under family coarsening"
)

ggsave(
  file.path(THEORY_OUT_DIR, "FigT5_regional_error_response_curves.png"),
  p_regional_curves, width = 9.2, height = 6.4, dpi = 300
)
ggsave(
  file.path(THEORY_OUT_DIR, "FigT5_regional_error_response_curves.pdf"),
  p_regional_curves, width = 9.2, height = 6.4
)
ggsave(
  file.path(THEORY_OUT_DIR, "FigT6_family_coarsening_response_curves.png"),
  p_family_curves, width = 9.2, height = 6.4, dpi = 300
)
ggsave(
  file.path(THEORY_OUT_DIR, "FigT6_family_coarsening_response_curves.pdf"),
  p_family_curves, width = 9.2, height = 6.4
)

# Console summary: compact and designed to be pasted into a lab notebook.
message("\nTop exploratory associations for regional-pool congeneric error:")
print(
  parameter_associations_top %>%
    filter(scenario == "regional_congeneric_error") %>%
    select(screen_label, predictor_label, spearman_rho, n)
)

message("\nTop exploratory associations for family-level coarsening:")
print(
  parameter_associations_top %>%
    filter(scenario == "family_level") %>%
    select(screen_label, predictor_label, spearman_rho, n)
)

message("\nParameter-space exploration complete. Outputs: ", THEORY_OUT_DIR)
