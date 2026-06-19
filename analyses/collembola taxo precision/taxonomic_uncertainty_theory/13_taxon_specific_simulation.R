# ============================================================
# 13_taxon_specific_simulation.R
# Taxon-specific theoretical simulations using real taxonomy only
# ============================================================
#
# This script does not use empirical scenario outcomes as a calibration target.
# It uses a group's observed regional taxonomy and basic sampling descriptors to
# generate synthetic communities, then applies identical theoretical scenarios.
# This allows an apples-to-apples comparison of vulnerability among taxa.

run_taxon_profile_world <- function(
  profile,
  scenario_grid,
  replicate_id = 1,
  seed = THEORY_SEED,
  community_overrides = list()
) {
  if (!is.list(profile) || is.null(profile$summary) || is.null(profile$regional_taxonomy)) {
    stop("`profile` must be returned by build_empirical_taxon_profile().", call. = FALSE)
  }
  if (!is.data.frame(scenario_grid) || nrow(scenario_grid) == 0) {
    stop("`scenario_grid` must be a non-empty scenario table.", call. = FALSE)
  }

  taxon <- profile$summary$taxon[[1]]
  taxonomy <- profile$regional_taxonomy
  p <- profile$summary

  # The actual regional hierarchy is retained. Community structure remains
  # synthetic, so this is not a calibration to the observed scenario outcomes.
  default_community <- list(
    n_sites = as.integer(p$n_sites[[1]]),
    n_observed_species = as.integer(p$n_observed_species_from_community[[1]]),
    occupancy_mean = clamp(p$observed_mean_occupancy[[1]], 0.01, 0.95),
    occupancy_precision = 6,
    mean_abundance = 8,
    abundance_nb_size = 1.5,
    occupancy_abundance_slope = 0.6,
    observed_selection_bias = 0,
    environmental_gradient_strength = 0
  )
  community_parameters <- utils::modifyList(default_community, community_overrides)
  community_parameters$n_observed_species <- min(
    as.integer(community_parameters$n_observed_species),
    nrow(taxonomy)
  )

  community <- do.call(
    generate_true_community,
    c(list(taxonomy = taxonomy), community_parameters, list(seed = seed + 1L))
  )

  tax_architecture <- summarise_taxonomic_architecture(taxonomy)
  true_metrics <- community_metrics(community$comm_true)

  profile_metadata <- p %>%
    select(
      taxon,
      n_sites,
      n_observed_species_from_community,
      regional_observed_species_ratio,
      regional_mean_species_per_genus,
      regional_mean_genera_per_family,
      observed_mean_occupancy
    )

  scenario_metadata <- scenario_grid %>% select(-baseline_process, -target_process)

  results <- purrr::map_dfr(seq_len(nrow(scenario_grid)), function(i) {
    scenario_row <- scenario_grid[i, , drop = FALSE]

    pair <- run_scenario_pair(
      comm_true = community$comm_true,
      taxonomy = taxonomy,
      baseline_process = scenario_row$baseline_process[[1]],
      target_process = scenario_row$target_process[[1]],
      seed = seed + 1000L + i * 100L
    )

    metadata <- bind_cols(
      profile_metadata,
      tibble(replicate_id = replicate_id),
      scenario_metadata[i, , drop = FALSE]
    )

    bind_rows(
      pair$blowes %>% bind_cols(metadata) %>% mutate(output_type = "blowes"),
      pair$effect_long %>% bind_cols(metadata) %>% mutate(output_type = "effect"),
      pair$stability %>% bind_cols(metadata) %>% mutate(output_type = "stability"),
      bind_cols(pair$target$error_summary, pair$target$workflow_summary) %>%
        bind_cols(metadata) %>%
        mutate(output_type = "process")
    )
  })

  list(
    profile = bind_cols(
      profile_metadata,
      tax_architecture %>% rename_with(~ paste0("simulated_", .x)),
      tibble(
        replicate_id = replicate_id,
        true_alpha_q0 = true_metrics$alpha_q0,
        true_gamma = true_metrics$gamma,
        true_mean_occupancy = true_metrics$mean_occupancy
      )
    ),
    results = results
  )
}

run_taxon_profile_experiment <- function(
  profiles,
  scenario_grid = make_controlled_scenario_grid(),
  n_replicates_per_taxon = 20,
  seed = THEORY_SEED,
  community_overrides = list(),
  progress = TRUE
) {
  if (!is.list(profiles) || length(profiles) == 0) {
    stop("`profiles` must be a non-empty named list of profile objects.", call. = FALSE)
  }
  assert_scalar(n_replicates_per_taxon, "n_replicates_per_taxon", lower = 1, integer = TRUE)

  jobs <- tidyr::crossing(
    profile_index = seq_along(profiles),
    replicate_id = seq_len(n_replicates_per_taxon)
  )

  out <- vector("list", nrow(jobs))
  for (job_i in seq_len(nrow(jobs))) {
    profile_i <- profiles[[jobs$profile_index[[job_i]]]]
    taxon <- profile_i$summary$taxon[[1]]

    if (progress && (job_i == 1 || job_i %% 10 == 0 || job_i == nrow(jobs))) {
      message(sprintf("Taxon-profile simulations: %d / %d (%s)", job_i, nrow(jobs), taxon))
    }

    out[[job_i]] <- run_taxon_profile_world(
      profile = profile_i,
      scenario_grid = scenario_grid,
      replicate_id = jobs$replicate_id[[job_i]],
      seed = seed + job_i * 10000L,
      community_overrides = community_overrides
    )
  }

  list(
    profiles = bind_rows(purrr::map(out, "profile")),
    results = bind_rows(purrr::map(out, "results")),
    scenario_grid = scenario_grid
  )
}

summarise_taxon_profile_experiment <- function(experiment) {
  if (!is.list(experiment) || is.null(experiment$results)) {
    stop("`experiment` must be returned by run_taxon_profile_experiment().", call. = FALSE)
  }

  blowes_summary <- experiment$results %>%
    filter(output_type == "blowes") %>%
    group_by(
      taxon, scenario, scenario_label, scenario_domain, mechanism,
      intensity_type, requested_intensity, intensity_label
    ) %>%
    summarise(
      delta_alpha_med = safe_median(delta_alpha_pct),
      delta_alpha_p10 = safe_quantile(delta_alpha_pct, 0.10),
      delta_alpha_p90 = safe_quantile(delta_alpha_pct, 0.90),
      delta_gamma_med = safe_median(delta_gamma_pct),
      delta_gamma_p10 = safe_quantile(delta_gamma_pct, 0.10),
      delta_gamma_p90 = safe_quantile(delta_gamma_pct, 0.90),
      delta_occupancy_med = safe_median(delta_occupancy_pct),
      delta_occupancy_p10 = safe_quantile(delta_occupancy_pct, 0.10),
      delta_occupancy_p90 = safe_quantile(delta_occupancy_pct, 0.90),
      beta_signature = names(sort(table(beta_signature), decreasing = TRUE))[1],
      n_replicates = dplyr::n(),
      .groups = "drop"
    )

  # `experiment$results` is a row-bound table that also contains the long
  # effect output, which already has a column named `metric`. After filtering
  # to stability rows this inherited column is empty, but it must be removed
  # before creating the new metric label through pivot_longer().
  stability_summary <- experiment$results %>%
    filter(output_type == "stability") %>%
    select(-any_of("metric")) %>%
    pivot_longer(
      cols = c(alpha_q0_spearman, bray_spearman, jaccard_spearman, sorensen_spearman),
      names_to = "metric", values_to = "stability"
    ) %>%
    group_by(
      taxon, scenario, scenario_label, scenario_domain, mechanism,
      intensity_type, requested_intensity, intensity_label, metric
    ) %>%
    summarise(
      stability_med = safe_median(stability),
      stability_p10 = safe_quantile(stability, 0.10),
      stability_p90 = safe_quantile(stability, 0.90),
      .groups = "drop"
    )

  list(blowes_summary = blowes_summary, stability_summary = stability_summary)
}
