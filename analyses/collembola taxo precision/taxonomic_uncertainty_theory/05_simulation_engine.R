# ============================================================
# 05_simulation_engine.R
# One world, one scenario pair, and repeated simulations
# ============================================================

run_scenario_pair <- function(
  comm_true,
  taxonomy,
  baseline_process = list(),
  target_process = list(),
  seed = NULL
) {
  baseline <- apply_observation_process(
    comm_true = comm_true,
    taxonomy = taxonomy,
    process = baseline_process,
    seed = seed
  )

  target <- apply_observation_process(
    comm_true = comm_true,
    taxonomy = taxonomy,
    process = target_process,
    seed = if (is.null(seed)) NULL else seed + 10000L
  )

  baseline_metrics <- community_metrics(baseline$comm)
  target_metrics <- community_metrics(target$comm)

  comparison <- compare_to_baseline(baseline_metrics, target_metrics)
  stability <- metric_stability(baseline$comm, target$comm)

  list(
    baseline = baseline,
    target = target,
    baseline_metrics = baseline_metrics,
    target_metrics = target_metrics,
    effect_long = comparison$long,
    blowes = comparison$blowes,
    stability = stability
  )
}

default_scenario_library <- function(
  id_error_rate = 0.10,
  workflow_unresolved_rate = 0.20,
  rare_quantile = 0.10
) {
  truth <- list(workflow = "species")

  list(
    observed_congeneric_error = list(
      baseline_process = truth,
      target_process = list(
        id_error_rate = id_error_rate,
        candidate_pool = "observed",
        workflow = "species"
      )
    ),
    regional_congeneric_error = list(
      baseline_process = truth,
      target_process = list(
        id_error_rate = id_error_rate,
        candidate_pool = "regional",
        workflow = "species"
      )
    ),
    rare_weighted_error = list(
      baseline_process = truth,
      target_process = list(
        id_error_rate = id_error_rate,
        candidate_pool = "observed",
        id_rarity_exponent = 0.5,
        workflow = "species"
      )
    ),
    mixed_rtu_reporting = list(
      baseline_process = truth,
      target_process = list(
        workflow = "mixed_rtu",
        unresolved_rate = workflow_unresolved_rate
      )
    ),
    drop_unresolved = list(
      baseline_process = list(
        workflow = "mixed_rtu",
        unresolved_rate = workflow_unresolved_rate
      ),
      target_process = list(
        workflow = "drop_unresolved",
        unresolved_rate = workflow_unresolved_rate
      )
    ),
    rare_species_to_genus = list(
      baseline_process = truth,
      target_process = list(
        workflow = "rare_to_genus",
        rare_quantile = rare_quantile
      )
    ),
    difficult_genera_to_genus = list(
      baseline_process = truth,
      target_process = list(workflow = "difficult_to_genus")
    ),
    genus_level = list(
      baseline_process = truth,
      target_process = list(workflow = "genus")
    ),
    family_level = list(
      baseline_process = truth,
      target_process = list(workflow = "family")
    )
  )
}

run_one_theoretical_world <- function(
  architecture_parameters,
  community_parameters,
  scenario_library = default_scenario_library(),
  world_id = 1,
  replicate_id = 1,
  seed = NULL
) {
  assert_named_list(architecture_parameters, "architecture_parameters")
  assert_named_list(community_parameters, "community_parameters")
  assert_named_list(scenario_library, "scenario_library")

  architecture <- do.call(generate_taxonomic_architecture, c(architecture_parameters, list(seed = seed)))
  community <- do.call(
    generate_true_community,
    c(list(taxonomy = architecture), community_parameters, list(seed = if (is.null(seed)) NULL else seed + 1L))
  )

  architecture_summary <- summarise_taxonomic_architecture(architecture)
  true_metrics <- community_metrics(community$comm_true)

  scenario_results <- imap(scenario_library, function(scenario, scenario_name) {
    result <- run_scenario_pair(
      comm_true = community$comm_true,
      taxonomy = architecture,
      baseline_process = scenario$baseline_process,
      target_process = scenario$target_process,
      seed = if (is.null(seed)) NULL else seed + 1000L + match(scenario_name, names(scenario_library))
    )

    list(
      blowes = result$blowes %>%
        mutate(
          world_id = world_id,
          replicate_id = replicate_id,
          scenario = scenario_name
        ),
      effects = result$effect_long %>%
        mutate(
          world_id = world_id,
          replicate_id = replicate_id,
          scenario = scenario_name
        ),
      stability = result$stability %>%
        mutate(
          world_id = world_id,
          replicate_id = replicate_id,
          scenario = scenario_name
        ),
      process = bind_cols(
        result$target$error_summary,
        result$target$workflow_summary
      ) %>%
        mutate(
          world_id = world_id,
          replicate_id = replicate_id,
          scenario = scenario_name
        )
    )
  })

  list(
    world_profile = bind_cols(
      tibble(world_id = world_id, replicate_id = replicate_id),
      architecture_summary,
      tibble(
        n_observed_species = length(community$observed_species),
        true_alpha_q0 = true_metrics$alpha_q0,
        true_gamma = true_metrics$gamma,
        true_mean_occupancy = true_metrics$mean_occupancy,
        true_mean_bray_curtis = true_metrics$mean_bray_curtis
      )
    ),
    blowes = bind_rows(map(scenario_results, "blowes")),
    effects = bind_rows(map(scenario_results, "effects")),
    stability = bind_rows(map(scenario_results, "stability")),
    process = bind_rows(map(scenario_results, "process"))
  )
}

run_theoretical_experiment <- function(
  world_design,
  scenario_library = default_scenario_library(),
  n_replicates_per_world = 10,
  seed = THEORY_SEED,
  progress = TRUE
) {
  if (!is.data.frame(world_design) || nrow(world_design) == 0) {
    stop("`world_design` must be a non-empty data frame.", call. = FALSE)
  }
  assert_scalar(n_replicates_per_world, "n_replicates_per_world", lower = 1, integer = TRUE)

  expected_columns <- c(
    "n_regional_species", "n_genera", "n_families",
    "species_per_genus_concentration", "genera_per_family_concentration",
    "difficult_genus_fraction", "n_sites", "n_observed_species",
    "occupancy_mean", "occupancy_precision", "mean_abundance",
    "abundance_nb_size", "occupancy_abundance_slope",
    "observed_selection_bias", "environmental_gradient_strength"
  )
  missing <- setdiff(expected_columns, names(world_design))
  if (length(missing) > 0) {
    stop("`world_design` lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  jobs <- tidyr::crossing(
    world_row = seq_len(nrow(world_design)),
    replicate_id = seq_len(n_replicates_per_world)
  )

  results <- vector("list", nrow(jobs))

  for (job_index in seq_len(nrow(jobs))) {
    row_index <- jobs$world_row[[job_index]]
    replicate_id <- jobs$replicate_id[[job_index]]
    world <- world_design[row_index, , drop = FALSE]

    architecture_parameters <- list(
      n_regional_species = world$n_regional_species[[1]],
      n_genera = world$n_genera[[1]],
      n_families = world$n_families[[1]],
      species_per_genus_concentration = world$species_per_genus_concentration[[1]],
      genera_per_family_concentration = world$genera_per_family_concentration[[1]],
      difficult_genus_fraction = world$difficult_genus_fraction[[1]]
    )

    community_parameters <- list(
      n_sites = world$n_sites[[1]],
      n_observed_species = world$n_observed_species[[1]],
      occupancy_mean = world$occupancy_mean[[1]],
      occupancy_precision = world$occupancy_precision[[1]],
      mean_abundance = world$mean_abundance[[1]],
      abundance_nb_size = world$abundance_nb_size[[1]],
      occupancy_abundance_slope = world$occupancy_abundance_slope[[1]],
      observed_selection_bias = world$observed_selection_bias[[1]],
      environmental_gradient_strength = world$environmental_gradient_strength[[1]]
    )

    if (progress && (job_index == 1 || job_index %% 25 == 0 || job_index == nrow(jobs))) {
      message(sprintf("Theoretical simulations: %d / %d", job_index, nrow(jobs)))
    }

    results[[job_index]] <- run_one_theoretical_world(
      architecture_parameters = architecture_parameters,
      community_parameters = community_parameters,
      scenario_library = scenario_library,
      world_id = world$world_id[[1]],
      replicate_id = replicate_id,
      seed = seed + job_index * 100L
    )
  }

  list(
    design = world_design,
    world_profiles = bind_rows(map(results, "world_profile")),
    blowes = bind_rows(map(results, "blowes")),
    effects = bind_rows(map(results, "effects")),
    stability = bind_rows(map(results, "stability")),
    process = bind_rows(map(results, "process")),
    scenario_library = scenario_library
  )
}

summarise_theoretical_experiment <- function(experiment) {
  blowes_summary <- experiment$blowes %>%
    group_by(world_id, scenario) %>%
    summarise(
      delta_alpha_med = median(delta_alpha_pct, na.rm = TRUE),
      delta_alpha_p10 = safe_quantile(delta_alpha_pct, 0.10),
      delta_alpha_p90 = safe_quantile(delta_alpha_pct, 0.90),
      delta_gamma_med = median(delta_gamma_pct, na.rm = TRUE),
      delta_gamma_p10 = safe_quantile(delta_gamma_pct, 0.10),
      delta_gamma_p90 = safe_quantile(delta_gamma_pct, 0.90),
      delta_occupancy_med = median(delta_occupancy_pct, na.rm = TRUE),
      beta_signature = names(sort(table(beta_signature), decreasing = TRUE))[1],
      .groups = "drop"
    )

  blowes_summary %>%
    left_join(
      experiment$world_profiles %>% distinct(world_id, .keep_all = TRUE),
      by = "world_id"
    )
}
