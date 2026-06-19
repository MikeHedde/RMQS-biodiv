# ============================================================
# 10_controlled_one_factor_experiment.R
# Controlled one-factor-at-a-time (OFAT) experiment
# ============================================================
#
# Purpose
# -------
# This script creates a deliberately controlled theoretical experiment.
# It varies one candidate mechanism at a time while holding the others at a
# reference value. It is not a calibration to Collembola or any other taxon.
#
# The four controlled axes are:
#   1. regional / observed species ratio;
#   2. mean species per genus;
#   3. mean genera per family;
#   4. initial mean occupancy.
#
# The third axis is expressed as genera per family rather than species per
# family so that it remains independent of mean species per genus. Mean species
# per family is an emergent property of these two quantities.
#
# All identification-error scenarios currently use independent individual-level
# reassignment. Mechanisms differ through the candidate pool and the taxa that
# are preferentially affected. Site-consistent and species-consistent errors
# should be treated as a later extension of the observation engine.

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

controlled_reference_parameters <- function() {
  list(
    # Four focal axes.
    regional_observed_ratio = 4,
    mean_species_per_genus = 4,
    mean_genera_per_family = 4,
    occupancy_mean = 0.20,

    # Fixed architecture shape.
    n_observed_species = 80L,
    species_per_genus_concentration = 2,
    genera_per_family_concentration = 2,
    difficult_genus_fraction = 0.20,

    # Fixed ecological context.
    n_sites = 70L,
    occupancy_precision = 6,
    mean_abundance = 8,
    abundance_nb_size = 1.5,
    occupancy_abundance_slope = 0.6,
    observed_selection_bias = 0,
    environmental_gradient_strength = 0
  )
}

controlled_factor_levels <- function() {
  list(
    regional_observed_ratio = c(1.5, 3, 6, 12),
    mean_species_per_genus = c(1.5, 2, 4, 8),
    mean_genera_per_family = c(1.5, 2, 4, 8),
    occupancy_mean = c(0.05, 0.15, 0.30, 0.50)
  )
}

controlled_factor_labels <- function() {
  c(
    regional_observed_ratio = "Regional / observed species ratio",
    mean_species_per_genus = "Mean species per genus",
    mean_genera_per_family = "Mean genera per family",
    occupancy_mean = "Initial mean occupancy"
  )
}

# Convert interpretable continuous design parameters into a valid hierarchy.
derive_controlled_world_parameters <- function(parameters) {
  assert_named_list(parameters, "parameters")

  required <- c(
    "regional_observed_ratio", "mean_species_per_genus",
    "mean_genera_per_family", "occupancy_mean", "n_observed_species",
    "species_per_genus_concentration", "genera_per_family_concentration",
    "difficult_genus_fraction", "n_sites", "occupancy_precision",
    "mean_abundance", "abundance_nb_size", "occupancy_abundance_slope",
    "observed_selection_bias", "environmental_gradient_strength"
  )
  missing <- setdiff(required, names(parameters))
  if (length(missing) > 0) {
    stop("Missing controlled-world parameters: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  n_observed <- as.integer(round(parameters$n_observed_species))
  n_regional <- max(n_observed, as.integer(round(n_observed * parameters$regional_observed_ratio)))
  n_genera <- max(1L, min(n_regional, as.integer(round(n_regional / parameters$mean_species_per_genus))))
  n_families <- max(1L, min(n_genera, as.integer(round(n_genera / parameters$mean_genera_per_family))))

  tibble(
    n_regional_species = n_regional,
    n_genera = n_genera,
    n_families = n_families,
    species_per_genus_concentration = parameters$species_per_genus_concentration,
    genera_per_family_concentration = parameters$genera_per_family_concentration,
    difficult_genus_fraction = parameters$difficult_genus_fraction,
    n_sites = as.integer(round(parameters$n_sites)),
    n_observed_species = n_observed,
    occupancy_mean = parameters$occupancy_mean,
    occupancy_precision = parameters$occupancy_precision,
    mean_abundance = parameters$mean_abundance,
    abundance_nb_size = parameters$abundance_nb_size,
    occupancy_abundance_slope = parameters$occupancy_abundance_slope,
    observed_selection_bias = parameters$observed_selection_bias,
    environmental_gradient_strength = parameters$environmental_gradient_strength,

    # Inputs retained for direct interpretation in outputs.
    regional_observed_ratio_input = parameters$regional_observed_ratio,
    mean_species_per_genus_input = parameters$mean_species_per_genus,
    mean_genera_per_family_input = parameters$mean_genera_per_family,
    occupancy_mean_input = parameters$occupancy_mean,
    mean_species_per_family_input = n_regional / n_families
  )
}

make_controlled_one_factor_design <- function(
  reference = controlled_reference_parameters(),
  factor_levels = controlled_factor_levels()
) {
  assert_named_list(reference, "reference")
  assert_named_list(factor_levels, "factor_levels")

  valid_factors <- names(controlled_factor_labels())
  unknown <- setdiff(names(factor_levels), valid_factors)
  if (length(unknown) > 0) {
    stop("Unknown controlled factors: ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  design <- purrr::imap_dfr(factor_levels, function(levels, factor_name) {
    if (!is.numeric(levels) || length(levels) < 2 || any(!is.finite(levels))) {
      stop("Each factor-level vector must contain at least two finite numeric values.", call. = FALSE)
    }

    purrr::map_dfr(levels, function(level) {
      current <- reference
      current[[factor_name]] <- level

      derive_controlled_world_parameters(current) %>%
        mutate(
          varied_parameter = factor_name,
          factor_label = unname(controlled_factor_labels()[factor_name]),
          factor_value = level,
          factor_value_label = format(level, trim = TRUE, scientific = FALSE)
        )
    })
  }) %>%
    mutate(world_id = row_number()) %>%
    relocate(world_id, varied_parameter, factor_label, factor_value, factor_value_label)

  design
}

make_controlled_scenario_grid <- function(
  id_error_rates = c(0.01, 0.03, 0.05, 0.10, 0.20),
  reporting_rates = c(0.01, 0.03, 0.05, 0.10, 0.20),
  rare_quantile = 0.10,
  difficult_multiplier = 4
) {
  if (any(id_error_rates < 0 | id_error_rates > 1)) {
    stop("`id_error_rates` must be within [0, 1].", call. = FALSE)
  }
  if (any(reporting_rates < 0 | reporting_rates > 1)) {
    stop("`reporting_rates` must be within [0, 1].", call. = FALSE)
  }

  species_truth <- list(workflow = "species")

  make_row <- function(
    scenario, scenario_label, scenario_domain, mechanism,
    intensity_type, requested_intensity, baseline_process, target_process
  ) {
    tibble(
      scenario = scenario,
      scenario_label = scenario_label,
      scenario_domain = scenario_domain,
      mechanism = mechanism,
      intensity_type = intensity_type,
      requested_intensity = requested_intensity,
      intensity_label = ifelse(
        is.na(requested_intensity),
        "fixed",
        paste0(formatC(100 * requested_intensity, format = "fg", digits = 4), "%")
      ),
      baseline_process = list(baseline_process),
      target_process = list(target_process)
    )
  }

  id_rows <- purrr::map_dfr(id_error_rates, function(rate) {
    bind_rows(
      make_row(
        scenario = paste0("id_observed_uniform_", round(100 * rate), "pct"),
        scenario_label = "Observed-pool congeneric error",
        scenario_domain = "identification_error",
        mechanism = "observed_uniform",
        intensity_type = "identification_error_rate",
        requested_intensity = rate,
        baseline_process = species_truth,
        target_process = list(
          id_error_rate = rate,
          candidate_pool = "observed",
          workflow = "species"
        )
      ),
      make_row(
        scenario = paste0("id_regional_uniform_", round(100 * rate), "pct"),
        scenario_label = "Regional-pool congeneric error",
        scenario_domain = "identification_error",
        mechanism = "regional_uniform",
        intensity_type = "identification_error_rate",
        requested_intensity = rate,
        baseline_process = species_truth,
        target_process = list(
          id_error_rate = rate,
          candidate_pool = "regional",
          workflow = "species"
        )
      ),
      make_row(
        scenario = paste0("id_observed_rare_weighted_", round(100 * rate), "pct"),
        scenario_label = "Rare-weighted observed-pool error",
        scenario_domain = "identification_error",
        mechanism = "observed_rare_weighted",
        intensity_type = "identification_error_rate",
        requested_intensity = rate,
        baseline_process = species_truth,
        target_process = list(
          id_error_rate = rate,
          candidate_pool = "observed",
          id_rarity_exponent = 0.5,
          workflow = "species"
        )
      ),
      make_row(
        scenario = paste0("id_regional_difficult_bias_", round(100 * rate), "pct"),
        scenario_label = "Difficult-genus regional-pool error",
        scenario_domain = "identification_error",
        mechanism = "regional_difficult_bias",
        intensity_type = "identification_error_rate",
        requested_intensity = rate,
        baseline_process = species_truth,
        target_process = list(
          id_error_rate = rate,
          candidate_pool = "regional",
          id_difficulty_multiplier = difficult_multiplier,
          workflow = "species"
        )
      )
    )
  })

  reporting_rows <- purrr::map_dfr(reporting_rates, function(rate) {
    bind_rows(
      make_row(
        scenario = paste0("workflow_mixed_rtu_", round(100 * rate), "pct"),
        scenario_label = "Mixed RTU reporting",
        scenario_domain = "reporting_workflow",
        mechanism = "mixed_rtu",
        intensity_type = "unresolved_reporting_rate",
        requested_intensity = rate,
        baseline_process = species_truth,
        target_process = list(
          workflow = "mixed_rtu",
          unresolved_rate = rate
        )
      ),
      make_row(
        scenario = paste0("workflow_drop_unresolved_", round(100 * rate), "pct"),
        scenario_label = "Dropping unresolved records",
        scenario_domain = "reporting_workflow",
        mechanism = "drop_unresolved",
        intensity_type = "unresolved_reporting_rate",
        requested_intensity = rate,
        baseline_process = list(
          workflow = "mixed_rtu",
          unresolved_rate = rate
        ),
        target_process = list(
          workflow = "drop_unresolved",
          unresolved_rate = rate
        )
      )
    )
  })

  fixed_rows <- bind_rows(
    make_row(
      scenario = "workflow_rare_species_to_genus",
      scenario_label = "Rare taxa reported at genus",
      scenario_domain = "fixed_workflow",
      mechanism = "rare_to_genus",
      intensity_type = "fixed_workflow",
      requested_intensity = NA_real_,
      baseline_process = species_truth,
      target_process = list(workflow = "rare_to_genus", rare_quantile = rare_quantile)
    ),
    make_row(
      scenario = "workflow_difficult_genera_to_genus",
      scenario_label = "Difficult genera reported at genus",
      scenario_domain = "fixed_workflow",
      mechanism = "difficult_to_genus",
      intensity_type = "fixed_workflow",
      requested_intensity = NA_real_,
      baseline_process = species_truth,
      target_process = list(workflow = "difficult_to_genus")
    ),
    make_row(
      scenario = "workflow_genus_level",
      scenario_label = "Genus-level aggregation",
      scenario_domain = "fixed_workflow",
      mechanism = "genus_level",
      intensity_type = "fixed_workflow",
      requested_intensity = NA_real_,
      baseline_process = species_truth,
      target_process = list(workflow = "genus")
    ),
    make_row(
      scenario = "workflow_family_level",
      scenario_label = "Family-level aggregation",
      scenario_domain = "fixed_workflow",
      mechanism = "family_level",
      intensity_type = "fixed_workflow",
      requested_intensity = NA_real_,
      baseline_process = species_truth,
      target_process = list(workflow = "family")
    )
  )

  bind_rows(id_rows, reporting_rows, fixed_rows)
}

run_controlled_one_factor_world <- function(
  design_row,
  scenario_grid,
  replicate_id,
  seed
) {
  if (!is.data.frame(design_row) || nrow(design_row) != 1) {
    stop("`design_row` must be a one-row data frame.", call. = FALSE)
  }

  architecture_parameters <- list(
    n_regional_species = design_row$n_regional_species[[1]],
    n_genera = design_row$n_genera[[1]],
    n_families = design_row$n_families[[1]],
    species_per_genus_concentration = design_row$species_per_genus_concentration[[1]],
    genera_per_family_concentration = design_row$genera_per_family_concentration[[1]],
    difficult_genus_fraction = design_row$difficult_genus_fraction[[1]]
  )

  community_parameters <- list(
    n_sites = design_row$n_sites[[1]],
    n_observed_species = design_row$n_observed_species[[1]],
    occupancy_mean = design_row$occupancy_mean[[1]],
    occupancy_precision = design_row$occupancy_precision[[1]],
    mean_abundance = design_row$mean_abundance[[1]],
    abundance_nb_size = design_row$abundance_nb_size[[1]],
    occupancy_abundance_slope = design_row$occupancy_abundance_slope[[1]],
    observed_selection_bias = design_row$observed_selection_bias[[1]],
    environmental_gradient_strength = design_row$environmental_gradient_strength[[1]]
  )

  architecture <- do.call(
    generate_taxonomic_architecture,
    c(architecture_parameters, list(seed = seed))
  )
  community <- do.call(
    generate_true_community,
    c(list(taxonomy = architecture), community_parameters, list(seed = seed + 1L))
  )

  architecture_summary <- summarise_taxonomic_architecture(architecture)
  true_metrics <- community_metrics(community$comm_true)

  design_metadata <- design_row %>%
    select(-matches("^species_per_genus_concentration$|^genera_per_family_concentration$"))

  world_profile <- bind_cols(
    design_metadata,
    tibble(replicate_id = replicate_id),
    architecture_summary,
    tibble(
      true_alpha_q0 = true_metrics$alpha_q0,
      true_gamma = true_metrics$gamma,
      true_mean_occupancy = true_metrics$mean_occupancy,
      true_mean_bray_curtis = true_metrics$mean_bray_curtis
    )
  )

  scenario_metadata <- scenario_grid %>%
    select(-baseline_process, -target_process)

  results <- purrr::map_dfr(seq_len(nrow(scenario_grid)), function(scenario_index) {
    scenario_row <- scenario_grid[scenario_index, , drop = FALSE]
    scenario_seed <- seed + 1000L + scenario_index * 100L

    pair <- run_scenario_pair(
      comm_true = community$comm_true,
      taxonomy = architecture,
      baseline_process = scenario_row$baseline_process[[1]],
      target_process = scenario_row$target_process[[1]],
      seed = scenario_seed
    )

    metadata <- bind_cols(
      design_metadata,
      tibble(replicate_id = replicate_id),
      scenario_metadata[scenario_index, , drop = FALSE]
    )

    bind_rows(
      pair$blowes %>%
        bind_cols(metadata) %>%
        mutate(output_type = "blowes"),
      pair$effect_long %>%
        bind_cols(metadata) %>%
        mutate(output_type = "effect"),
      pair$stability %>%
        bind_cols(metadata) %>%
        mutate(output_type = "stability"),
      bind_cols(pair$target$error_summary, pair$target$workflow_summary) %>%
        bind_cols(metadata) %>%
        mutate(output_type = "process")
    )
  })

  list(world_profile = world_profile, combined_results = results)
}

run_controlled_one_factor_experiment <- function(
  design,
  scenario_grid,
  n_replicates_per_world = 5,
  seed = THEORY_SEED,
  progress = TRUE
) {
  if (!is.data.frame(design) || nrow(design) == 0) {
    stop("`design` must be a non-empty data frame.", call. = FALSE)
  }
  if (!is.data.frame(scenario_grid) || nrow(scenario_grid) == 0) {
    stop("`scenario_grid` must be a non-empty data frame.", call. = FALSE)
  }
  assert_scalar(n_replicates_per_world, "n_replicates_per_world", lower = 1, integer = TRUE)

  jobs <- tidyr::crossing(
    design_row = seq_len(nrow(design)),
    replicate_id = seq_len(n_replicates_per_world)
  )

  results <- vector("list", nrow(jobs))

  for (job_index in seq_len(nrow(jobs))) {
    row_index <- jobs$design_row[[job_index]]
    rep_index <- jobs$replicate_id[[job_index]]

    if (progress && (job_index == 1 || job_index %% 10 == 0 || job_index == nrow(jobs))) {
      message(sprintf("Controlled OFAT simulation: %d / %d", job_index, nrow(jobs)))
    }

    results[[job_index]] <- run_controlled_one_factor_world(
      design_row = design[row_index, , drop = FALSE],
      scenario_grid = scenario_grid,
      replicate_id = rep_index,
      seed = seed + job_index * 100000L
    )
  }

  combined <- bind_rows(purrr::map(results, "combined_results"))

  list(
    design = design,
    scenario_grid = scenario_grid,
    world_profiles = bind_rows(purrr::map(results, "world_profile")),
    blowes = combined %>% filter(output_type == "blowes") %>% select(-output_type),
    effects = combined %>% filter(output_type == "effect") %>% select(-output_type),
    stability = combined %>% filter(output_type == "stability") %>% select(-output_type),
    process = combined %>% filter(output_type == "process") %>% select(-output_type)
  )
}

summarise_controlled_one_factor_experiment <- function(experiment) {
  metadata <- c(
    "varied_parameter", "factor_label", "factor_value", "factor_value_label",
    "scenario", "scenario_label", "scenario_domain", "mechanism",
    "intensity_type", "requested_intensity", "intensity_label"
  )

  blowes_summary <- experiment$blowes %>%
    group_by(across(all_of(metadata))) %>%
    summarise(
      n_replicates = n(),
      delta_alpha_med = median(delta_alpha_pct, na.rm = TRUE),
      delta_alpha_p10 = safe_quantile(delta_alpha_pct, 0.10),
      delta_alpha_p90 = safe_quantile(delta_alpha_pct, 0.90),
      delta_gamma_med = median(delta_gamma_pct, na.rm = TRUE),
      delta_gamma_p10 = safe_quantile(delta_gamma_pct, 0.10),
      delta_gamma_p90 = safe_quantile(delta_gamma_pct, 0.90),
      delta_occupancy_med = median(delta_occupancy_pct, na.rm = TRUE),
      delta_occupancy_p10 = safe_quantile(delta_occupancy_pct, 0.10),
      delta_occupancy_p90 = safe_quantile(delta_occupancy_pct, 0.90),
      apparent_differentiation_probability = mean(beta_signature == "apparent_differentiation", na.rm = TRUE),
      apparent_homogenisation_probability = mean(beta_signature == "apparent_homogenisation", na.rm = TRUE),
      .groups = "drop"
    )

  effects_summary <- experiment$effects %>%
    group_by(across(all_of(metadata)), metric) %>%
    summarise(
      n_replicates = n(),
      relative_delta_med = median(relative_delta_pct, na.rm = TRUE),
      relative_delta_p10 = safe_quantile(relative_delta_pct, 0.10),
      relative_delta_p90 = safe_quantile(relative_delta_pct, 0.90),
      .groups = "drop"
    )

  stability_summary <- experiment$stability %>%
    group_by(across(all_of(metadata))) %>%
    summarise(
      n_replicates = n(),
      alpha_q0_spearman_med = median(alpha_q0_spearman, na.rm = TRUE),
      bray_spearman_med = median(bray_spearman, na.rm = TRUE),
      jaccard_spearman_med = median(jaccard_spearman, na.rm = TRUE),
      sorensen_spearman_med = median(sorensen_spearman, na.rm = TRUE),
      .groups = "drop"
    )

  process_summary <- experiment$process %>%
    group_by(across(all_of(metadata))) %>%
    summarise(
      realised_mean_error_med = safe_median(realised_mean_error),
      realised_unresolved_rate_med = safe_median(realised_unresolved_rate),
      eligible_source_species_med = safe_median(n_eligible_source_species),
      .groups = "drop"
    )

  endpoint_from_blowes <- blowes_summary %>%
    select(all_of(metadata), delta_alpha_med, delta_gamma_med, delta_occupancy_med) %>%
    pivot_longer(
      cols = c(delta_alpha_med, delta_gamma_med, delta_occupancy_med),
      names_to = "response",
      values_to = "response_median"
    )

  endpoint_from_effects <- effects_summary %>%
    transmute(
      across(all_of(metadata)),
      response = paste0("effect_", metric),
      response_median = relative_delta_med
    )

  endpoint_contrasts <- bind_rows(endpoint_from_blowes, endpoint_from_effects) %>%
    group_by(
      varied_parameter, factor_label, scenario, scenario_label,
      scenario_domain, mechanism, intensity_type, requested_intensity,
      intensity_label, response
    ) %>%
    arrange(factor_value, .by_group = TRUE) %>%
    summarise(
      low_factor_value = first(factor_value),
      high_factor_value = last(factor_value),
      low_response = first(response_median),
      high_response = last(response_median),
      high_minus_low = high_response - low_response,
      .groups = "drop"
    )

  list(
    blowes_summary = blowes_summary,
    effects_summary = effects_summary,
    stability_summary = stability_summary,
    process_summary = process_summary,
    endpoint_contrasts = endpoint_contrasts
  )
}
