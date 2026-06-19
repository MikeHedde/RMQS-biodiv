# ============================================================
# 07_empirical_bridge.R
# Generic bridge from an empirical scenario table to the
# alpha-gamma-occupancy space. It does NOT source empirical scripts.
# ============================================================

empirical_blowes_from_summary <- function(
  scenario_summary_by_iter,
  scenario_col = "scenario",
  family_col = "scenario_family",
  baseline_col = "baseline_scenario",
  unit_type_col = "unit_type",
  iter_col = "iter",
  alpha_col = "mean_local_q0",
  gamma_col = "gamma_taxon_units",
  tolerance_pct = 0.5
) {
  required <- c(
    scenario_col, family_col, baseline_col, unit_type_col,
    iter_col, alpha_col, gamma_col
  )
  missing <- setdiff(required, names(scenario_summary_by_iter))
  if (length(missing) > 0) {
    stop("Scenario summary lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  dat <- scenario_summary_by_iter %>%
    transmute(
      scenario_family = .data[[family_col]],
      scenario = .data[[scenario_col]],
      baseline_scenario = .data[[baseline_col]],
      unit_type = .data[[unit_type_col]],
      iter = .data[[iter_col]],
      alpha = .data[[alpha_col]],
      gamma = .data[[gamma_col]],
      occupancy = .data[[alpha_col]] / .data[[gamma_col]]
    )

  baselines <- dat %>%
    filter(scenario == baseline_scenario) %>%
    transmute(
      scenario_family,
      baseline_scenario = scenario,
      alpha_baseline = alpha,
      gamma_baseline = gamma,
      occupancy_baseline = occupancy
    ) %>%
    distinct()

  dat %>%
    left_join(baselines, by = c("scenario_family", "baseline_scenario")) %>%
    mutate(
      delta_alpha_pct = relative_change(alpha, alpha_baseline),
      delta_gamma_pct = relative_change(gamma, gamma_baseline),
      delta_occupancy_pct = relative_change(occupancy, occupancy_baseline),
      beta_signature = case_when(
        abs(delta_alpha_pct) < tolerance_pct & abs(delta_gamma_pct) < tolerance_pct ~ "no_meaningful_change",
        delta_occupancy_pct > tolerance_pct ~ "apparent_homogenisation",
        delta_occupancy_pct < -tolerance_pct ~ "apparent_differentiation",
        TRUE ~ "no_meaningful_change"
      )
    )
}

empirical_taxon_profile <- function(
  community_matrix,
  observed_taxonomy,
  regional_taxonomy = NULL
) {
  if (!is.matrix(community_matrix)) community_matrix <- as.matrix(community_matrix)

  required <- c("species", "genus", "family")
  if (!all(required %in% names(observed_taxonomy))) {
    stop("`observed_taxonomy` must contain species, genus and family.", call. = FALSE)
  }

  observed_taxonomy <- observed_taxonomy %>%
    filter(species %in% colnames(community_matrix)) %>%
    distinct(species, .keep_all = TRUE)

  regional_taxonomy <- regional_taxonomy %||% observed_taxonomy
  if (!all(required %in% names(regional_taxonomy))) {
    stop("`regional_taxonomy` must contain species, genus and family.", call. = FALSE)
  }

  architecture <- summarise_taxonomic_architecture(regional_taxonomy)
  metrics <- community_metrics(community_matrix)

  bind_cols(
    architecture,
    tibble(
      n_sites = nrow(community_matrix),
      n_observed_species = sum(colSums(community_matrix) > 0),
      observed_alpha_q0 = metrics$alpha_q0,
      observed_gamma = metrics$gamma,
      observed_mean_occupancy = metrics$mean_occupancy,
      observed_mean_bray_curtis = metrics$mean_bray_curtis
    )
  )
}

# Base-R compatibility helper for NULL coalescence.
`%||%` <- function(x, y) if (is.null(x)) y else x

select_theoretical_worlds <- function(
  theory_world_profiles,
  taxon_profile,
  relative_tolerances = list(
    n_regional_species = 0.25,
    n_genera = 0.25,
    n_families = 0.35,
    n_sites = 0.25,
    n_observed_species = 0.25
  )
) {
  if (nrow(taxon_profile) != 1) {
    stop("`taxon_profile` must have one row.", call. = FALSE)
  }

  keep <- rep(TRUE, nrow(theory_world_profiles))

  for (variable in names(relative_tolerances)) {
    if (!variable %in% names(theory_world_profiles) || !variable %in% names(taxon_profile)) next
    target <- taxon_profile[[variable]][1]
    tolerance <- relative_tolerances[[variable]]
    keep <- keep & abs(theory_world_profiles[[variable]] - target) <= tolerance * max(target, 1)
  }

  theory_world_profiles %>% filter(keep)
}
