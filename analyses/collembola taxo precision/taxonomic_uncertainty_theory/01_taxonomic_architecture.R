# ============================================================
# 01_taxonomic_architecture.R
# Regional family -> genus -> species architecture
# ============================================================

# Allocate `total` objects among `n_groups`, ensuring at least one
# object per group. Gamma-distributed weights control unevenness.
# concentration < 1: a few rich groups and many poor groups
# concentration = 1: geometric-like allocation
# concentration > 1: more even allocation
allocate_positive_counts <- function(total, n_groups, concentration = 1) {
  assert_scalar(total, "total", lower = 1, integer = TRUE)
  assert_scalar(n_groups, "n_groups", lower = 1, integer = TRUE)
  assert_scalar(concentration, "concentration", lower = 1e-8)

  if (n_groups > total) {
    stop("`n_groups` cannot exceed `total` when every group must contain one member.", call. = FALSE)
  }

  remaining <- as.integer(total - n_groups)
  if (remaining == 0L) return(rep.int(1L, n_groups))

  weights <- stats::rgamma(n_groups, shape = concentration, rate = 1)
  probs <- weights / sum(weights)
  as.integer(1L + stats::rmultinom(1L, size = remaining, prob = probs)[, 1])
}

generate_taxonomic_architecture <- function(
  n_regional_species,
  n_genera,
  n_families,
  species_per_genus_concentration = 1,
  genera_per_family_concentration = 1,
  difficult_genus_fraction = 0,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  assert_scalar(n_regional_species, "n_regional_species", lower = 1, integer = TRUE)
  assert_scalar(n_genera, "n_genera", lower = 1, integer = TRUE)
  assert_scalar(n_families, "n_families", lower = 1, integer = TRUE)
  assert_scalar(difficult_genus_fraction, "difficult_genus_fraction", lower = 0, upper = 1)

  if (n_genera > n_regional_species) {
    stop("`n_genera` cannot exceed `n_regional_species`.", call. = FALSE)
  }
  if (n_families > n_genera) {
    stop("`n_families` cannot exceed `n_genera`.", call. = FALSE)
  }

  n_genera_per_family <- allocate_positive_counts(
    total = n_genera,
    n_groups = n_families,
    concentration = genera_per_family_concentration
  )

  family <- paste0("fam_", sprintf("%03d", seq_len(n_families)))
  genus <- paste0("gen_", sprintf("%04d", seq_len(n_genera)))

  family_by_genus <- rep(family, times = n_genera_per_family)

  n_species_per_genus <- allocate_positive_counts(
    total = n_regional_species,
    n_groups = n_genera,
    concentration = species_per_genus_concentration
  )

  species <- paste0("sp_", sprintf("%05d", seq_len(n_regional_species)))
  genus_by_species <- rep(genus, times = n_species_per_genus)

  taxonomy <- tibble(
    species = species,
    genus = genus_by_species,
    family = family_by_genus[match(genus_by_species, genus)]
  )

  n_difficult <- round(difficult_genus_fraction * n_genera)
  difficult_genera <- if (n_difficult > 0) sample(genus, n_difficult, replace = FALSE) else character(0)

  taxonomy <- taxonomy %>%
    mutate(
      difficult_genus = genus %in% difficult_genera,
      species_per_genus = n_species_per_genus[match(genus, genus)],
      genera_per_family = n_genera_per_family[match(family, family)]
    )

  taxonomy
}

summarise_taxonomic_architecture <- function(taxonomy) {
  required <- c("species", "genus", "family")
  if (!all(required %in% names(taxonomy))) {
    stop("`taxonomy` must contain species, genus and family.", call. = FALSE)
  }

  sp_per_genus <- taxonomy %>%
    distinct(species, genus) %>%
    count(genus, name = "n_species")

  gen_per_family <- taxonomy %>%
    distinct(genus, family) %>%
    count(family, name = "n_genera")

  tibble(
    n_regional_species = n_distinct(taxonomy$species),
    n_genera = n_distinct(taxonomy$genus),
    n_families = n_distinct(taxonomy$family),
    mean_species_per_genus = mean(sp_per_genus$n_species),
    cv_species_per_genus = stats::sd(sp_per_genus$n_species) / mean(sp_per_genus$n_species),
    prop_monotypic_genera = mean(sp_per_genus$n_species == 1),
    mean_genera_per_family = mean(gen_per_family$n_genera),
    cv_genera_per_family = stats::sd(gen_per_family$n_genera) / mean(gen_per_family$n_genera),
    prop_difficult_genera = mean(taxonomy %>% distinct(genus, difficult_genus) %>% pull(difficult_genus))
  )
}
