# ============================================================
# 02_community_generator.R
# True ecological community before taxonomic error
# ============================================================

draw_observed_species <- function(
  taxonomy,
  n_observed_species,
  selection_bias = 0
) {
  n_regional <- nrow(taxonomy)
  assert_scalar(n_observed_species, "n_observed_species", lower = 1, upper = n_regional, integer = TRUE)

  latent_accessibility <- stats::rnorm(n_regional)
  weights <- exp(selection_bias * latent_accessibility)
  sample(taxonomy$species, n_observed_species, replace = FALSE, prob = weights)
}

generate_true_community <- function(
  taxonomy,
  n_sites = 70,
  n_observed_species = min(50, nrow(taxonomy)),
  occupancy_mean = 0.20,
  occupancy_precision = 6,
  mean_abundance = 8,
  abundance_nb_size = 1.5,
  occupancy_abundance_slope = 0.6,
  observed_selection_bias = 0,
  environmental_gradient_strength = 0,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  assert_scalar(n_sites, "n_sites", lower = 2, integer = TRUE)
  assert_scalar(occupancy_mean, "occupancy_mean", lower = 1e-5, upper = 1 - 1e-5)
  assert_scalar(occupancy_precision, "occupancy_precision", lower = 1e-5)
  assert_scalar(mean_abundance, "mean_abundance", lower = 1e-5)
  assert_scalar(abundance_nb_size, "abundance_nb_size", lower = 1e-5)
  assert_scalar(environmental_gradient_strength, "environmental_gradient_strength", lower = 0)

  observed_species <- draw_observed_species(
    taxonomy = taxonomy,
    n_observed_species = n_observed_species,
    selection_bias = observed_selection_bias
  )

  alpha_beta <- occupancy_mean * occupancy_precision
  beta_beta <- (1 - occupancy_mean) * occupancy_precision

  species_occupancy <- stats::rbeta(
    length(observed_species),
    shape1 = alpha_beta,
    shape2 = beta_beta
  )
  species_occupancy <- clamp(species_occupancy, 1e-4, 1 - 1e-4)

  site_gradient <- seq(-1, 1, length.out = n_sites)
  species_optimum <- stats::runif(length(observed_species), min = -1, max = 1)

  presence_probability <- outer(
    site_gradient,
    seq_along(observed_species),
    FUN = function(site_value, species_index) {
      base_logit <- stats::qlogis(species_occupancy[species_index])
      penalty <- environmental_gradient_strength * (site_value - species_optimum[species_index])^2
      stats::plogis(base_logit - penalty)
    }
  )

  presence <- matrix(
    stats::rbinom(
      n = n_sites * length(observed_species),
      size = 1,
      prob = as.vector(presence_probability)
    ),
    nrow = n_sites,
    ncol = length(observed_species),
    dimnames = list(
      paste0("site_", sprintf("%03d", seq_len(n_sites))),
      observed_species
    )
  )

  abundance_mean <- mean_abundance *
    (species_occupancy / mean(species_occupancy))^occupancy_abundance_slope
  abundance_mean <- pmax(abundance_mean, 1e-4)

  counts_observed <- matrix(
    stats::rnbinom(
      n = n_sites * length(observed_species),
      mu = rep(abundance_mean, each = n_sites),
      size = abundance_nb_size
    ),
    nrow = n_sites,
    ncol = length(observed_species),
    dimnames = dimnames(presence)
  )

  counts_observed[presence == 0] <- 0
  counts_observed[presence == 1 & counts_observed == 0] <- 1

  # Expand to the full regional pool. Species absent from the observed
  # network start at zero but can appear as labels after regional-pool errors.
  comm_true <- matrix(
    0,
    nrow = n_sites,
    ncol = nrow(taxonomy),
    dimnames = list(rownames(counts_observed), taxonomy$species)
  )
  comm_true[, observed_species] <- counts_observed

  species_state <- taxonomy %>%
    filter(species %in% observed_species) %>%
    mutate(
      latent_occupancy = species_occupancy[match(species, observed_species)],
      niche_optimum = species_optimum[match(species, observed_species)],
      total_abundance = colSums(comm_true)[species],
      realised_occupancy = colMeans(comm_true[, species, drop = FALSE] > 0)
    )

  list(
    comm_true = comm_true,
    observed_species = observed_species,
    site_data = tibble(site = rownames(comm_true), environmental_gradient = site_gradient),
    species_state = species_state
  )
}
