# ============================================================
# 06_experiment_design.R
# Parameter-space design, without taxon-specific calibration
# ============================================================

draw_uniform <- function(range, integer = FALSE) {
  if (length(range) != 2 || !all(is.finite(range)) || range[1] > range[2]) {
    stop("Ranges must be numeric vectors of length two: c(min, max).", call. = FALSE)
  }
  value <- stats::runif(1, min = range[1], max = range[2])
  if (integer) value <- as.integer(round(value))
  value
}

default_parameter_ranges <- function() {
  list(
    n_regional_species = c(50, 1000),
    n_genera = c(15, 350),
    n_families = c(5, 100),
    species_per_genus_concentration = c(0.15, 8),
    genera_per_family_concentration = c(0.15, 8),
    difficult_genus_fraction = c(0, 0.40),
    n_sites = c(30, 150),
    observed_fraction = c(0.05, 0.70),
    occupancy_mean = c(0.03, 0.60),
    occupancy_precision = c(0.5, 30),
    mean_abundance = c(1.5, 50),
    abundance_nb_size = c(0.25, 10),
    occupancy_abundance_slope = c(-0.25, 1.50),
    observed_selection_bias = c(-1, 1),
    environmental_gradient_strength = c(0, 3)
  )
}

make_world_design <- function(
  n_worlds = 100,
  ranges = default_parameter_ranges(),
  seed = THEORY_SEED
) {
  if (!is.null(seed)) set.seed(seed)
  assert_scalar(n_worlds, "n_worlds", lower = 1, integer = TRUE)
  assert_named_list(ranges, "ranges")

  required <- names(default_parameter_ranges())
  missing <- setdiff(required, names(ranges))
  if (length(missing) > 0) {
    stop("`ranges` lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  worlds <- vector("list", n_worlds)

  for (i in seq_len(n_worlds)) {
    # Redraw constrained hierarchy until species >= genera >= families.
    repeat {
      n_regional_species <- max(2L, draw_uniform(ranges$n_regional_species, integer = TRUE))
      n_genera <- max(1L, min(
        n_regional_species,
        draw_uniform(ranges$n_genera, integer = TRUE)
      ))
      n_families <- max(1L, min(
        n_genera,
        draw_uniform(ranges$n_families, integer = TRUE)
      ))
      if (n_regional_species >= n_genera && n_genera >= n_families) break
    }

    observed_fraction <- draw_uniform(ranges$observed_fraction)
    n_observed_species <- max(2L, min(
      n_regional_species,
      as.integer(round(n_regional_species * observed_fraction))
    ))

    worlds[[i]] <- tibble(
      world_id = i,
      n_regional_species = n_regional_species,
      n_genera = n_genera,
      n_families = n_families,
      species_per_genus_concentration = draw_uniform(ranges$species_per_genus_concentration),
      genera_per_family_concentration = draw_uniform(ranges$genera_per_family_concentration),
      difficult_genus_fraction = draw_uniform(ranges$difficult_genus_fraction),
      n_sites = max(2L, draw_uniform(ranges$n_sites, integer = TRUE)),
      n_observed_species = n_observed_species,
      occupancy_mean = draw_uniform(ranges$occupancy_mean),
      occupancy_precision = draw_uniform(ranges$occupancy_precision),
      mean_abundance = draw_uniform(ranges$mean_abundance),
      abundance_nb_size = draw_uniform(ranges$abundance_nb_size),
      occupancy_abundance_slope = draw_uniform(ranges$occupancy_abundance_slope),
      observed_selection_bias = draw_uniform(ranges$observed_selection_bias),
      environmental_gradient_strength = draw_uniform(ranges$environmental_gradient_strength)
    )
  }

  bind_rows(worlds)
}

make_targeted_world_design <- function(
  parameter_sets
) {
  if (!is.data.frame(parameter_sets)) {
    stop("`parameter_sets` must be a data frame.", call. = FALSE)
  }

  required <- c(
    "n_regional_species", "n_genera", "n_families",
    "species_per_genus_concentration", "genera_per_family_concentration",
    "difficult_genus_fraction", "n_sites", "n_observed_species",
    "occupancy_mean", "occupancy_precision", "mean_abundance",
    "abundance_nb_size", "occupancy_abundance_slope",
    "observed_selection_bias", "environmental_gradient_strength"
  )

  missing <- setdiff(required, names(parameter_sets))
  if (length(missing) > 0) {
    stop("`parameter_sets` lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  parameter_sets %>%
    mutate(world_id = dplyr::row_number()) %>%
    select(world_id, all_of(required))
}
