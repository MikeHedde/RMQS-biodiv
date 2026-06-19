# ============================================================
# 03_observation_process.R
# Identification error K and reporting / workflow process C
# ============================================================

active_species <- function(comm) {
  colnames(comm)[colSums(comm) > 0]
}

build_error_probabilities <- function(
  comm,
  taxonomy,
  eligible_species,
  target_mean_error = 0,
  rarity_exponent = 0,
  difficulty_multiplier = 1,
  max_error = 0.95
) {
  assert_scalar(target_mean_error, "target_mean_error", lower = 0, upper = 1)
  assert_scalar(rarity_exponent, "rarity_exponent", lower = 0)
  assert_scalar(difficulty_multiplier, "difficulty_multiplier", lower = 0)
  assert_scalar(max_error, "max_error", lower = 0, upper = 1)

  totals <- colSums(comm)
  p <- stats::setNames(rep(0, ncol(comm)), colnames(comm))

  eligible_species <- intersect(eligible_species, names(totals))
  eligible_species <- eligible_species[totals[eligible_species] > 0]

  if (target_mean_error == 0 || length(eligible_species) == 0) {
    return(p)
  }

  tax <- taxonomy %>% filter(species %in% eligible_species)
  raw <- (totals[tax$species] + 1)^(-rarity_exponent)
  raw <- raw * ifelse(tax$difficult_genus, difficulty_multiplier, 1)

  weighted_mean_raw <- stats::weighted.mean(raw, w = totals[tax$species])
  if (!is.finite(weighted_mean_raw) || weighted_mean_raw <= 0) return(p)

  p_eligible <- pmin(target_mean_error * raw / weighted_mean_raw, max_error)
  p[tax$species] <- p_eligible
  p
}

candidate_species_for_error <- function(
  source_species,
  taxonomy,
  candidate_pool,
  observed_pool
) {
  source_row <- taxonomy %>% filter(species == source_species)
  if (nrow(source_row) != 1) return(character(0))

  candidates <- taxonomy %>%
    filter(genus == source_row$genus, species != source_species) %>%
    pull(species)

  if (candidate_pool == "observed") {
    candidates <- intersect(candidates, observed_pool)
  } else if (candidate_pool != "regional") {
    stop("`candidate_pool` must be 'observed' or 'regional'.", call. = FALSE)
  }

  candidates
}

apply_identification_error <- function(
  comm,
  taxonomy,
  target_mean_error = 0,
  candidate_pool = c("observed", "regional"),
  rarity_exponent = 0,
  difficulty_multiplier = 1,
  max_error = 0.95,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  candidate_pool <- match.arg(candidate_pool)

  if (!all(colnames(comm) %in% taxonomy$species)) {
    stop("All columns of `comm` must be named with taxonomy$species.", call. = FALSE)
  }

  observed_pool <- active_species(comm)

  candidate_lookup <- setNames(
    lapply(observed_pool, candidate_species_for_error,
           taxonomy = taxonomy,
           candidate_pool = candidate_pool,
           observed_pool = observed_pool),
    observed_pool
  )

  eligible <- names(candidate_lookup)[lengths(candidate_lookup) > 0]

  p_error <- build_error_probabilities(
    comm = comm,
    taxonomy = taxonomy,
    eligible_species = eligible,
    target_mean_error = target_mean_error,
    rarity_exponent = rarity_exponent,
    difficulty_multiplier = difficulty_multiplier,
    max_error = max_error
  )

  out <- comm
  n_moved_total <- 0

  # Source counts are always taken from `comm`, not from progressively
  # modified `out`, preventing artificial cascading errors.
  for (source in eligible) {
    p_i <- p_error[[source]]
    if (!is.finite(p_i) || p_i <= 0) next

    candidates <- candidate_lookup[[source]]
    source_counts <- comm[, source]

    for (site_index in seq_len(nrow(comm))) {
      n_source <- source_counts[[site_index]]
      if (n_source <= 0) next

      n_swap <- stats::rbinom(1, size = n_source, prob = p_i)
      if (n_swap == 0) next

      allocation <- stats::rmultinom(
        n = 1,
        size = n_swap,
        prob = rep(1 / length(candidates), length(candidates))
      )[, 1]

      out[site_index, source] <- out[site_index, source] - n_swap
      out[site_index, candidates] <- out[site_index, candidates] + allocation
      n_moved_total <- n_moved_total + n_swap
    }
  }

  out <- out[, colSums(out) > 0, drop = FALSE]

  error_summary <- tibble(
    requested_mean_error = target_mean_error,
    realised_mean_error = n_moved_total / sum(comm),
    candidate_pool = candidate_pool,
    rarity_exponent = rarity_exponent,
    difficulty_multiplier = difficulty_multiplier,
    n_eligible_source_species = length(eligible),
    n_moved_individuals = n_moved_total
  )

  list(comm = out, error_summary = error_summary)
}

aggregate_by_taxonomic_rank <- function(comm, taxonomy, rank = c("genus", "family")) {
  rank <- match.arg(rank)
  tax <- taxonomy %>% filter(species %in% colnames(comm))
  groups <- tax[[rank]][match(colnames(comm), tax$species)]

  out <- t(rowsum(t(comm), group = groups, reorder = TRUE))
  colnames(out) <- paste0(rank, "::", colnames(out))
  out[, colSums(out) > 0, drop = FALSE]
}

build_reporting_probabilities <- function(
  comm,
  taxonomy,
  unresolved_rate = 0,
  rarity_exponent = 0,
  difficulty_multiplier = 1,
  max_unresolved = 1
) {
  build_error_probabilities(
    comm = comm,
    taxonomy = taxonomy,
    eligible_species = active_species(comm),
    target_mean_error = unresolved_rate,
    rarity_exponent = rarity_exponent,
    difficulty_multiplier = difficulty_multiplier,
    max_error = max_unresolved
  )
}

report_mixed_rtu <- function(
  comm,
  taxonomy,
  unresolved_rate = 0,
  rarity_exponent = 0,
  difficulty_multiplier = 1,
  max_unresolved = 1,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  taxonomy_active <- taxonomy %>% filter(species %in% colnames(comm))
  p_unresolved <- build_reporting_probabilities(
    comm = comm,
    taxonomy = taxonomy_active,
    unresolved_rate = unresolved_rate,
    rarity_exponent = rarity_exponent,
    difficulty_multiplier = difficulty_multiplier,
    max_unresolved = max_unresolved
  )

  species_labels <- paste0("species::", colnames(comm))
  genus_labels <- paste0("genus::", unique(taxonomy_active$genus))

  out <- matrix(
    0,
    nrow = nrow(comm),
    ncol = length(species_labels) + length(genus_labels),
    dimnames = list(rownames(comm), c(species_labels, genus_labels))
  )

  for (species in colnames(comm)) {
    genus <- taxonomy_active$genus[match(species, taxonomy_active$species)]
    counts <- comm[, species]
    p_i <- p_unresolved[[species]]

    unresolved <- stats::rbinom(nrow(comm), size = counts, prob = p_i)
    out[, paste0("species::", species)] <- counts - unresolved
    out[, paste0("genus::", genus)] <- out[, paste0("genus::", genus)] + unresolved
  }

  out <- out[, colSums(out) > 0, drop = FALSE]

  list(
    comm = out,
    p_unresolved = p_unresolved,
    realised_unresolved_rate = 1 - sum(out[, grepl("^species::", colnames(out)), drop = FALSE]) / sum(comm)
  )
}

apply_reporting_workflow <- function(
  comm,
  taxonomy,
  workflow = c(
    "species", "mixed_rtu", "drop_unresolved",
    "rare_to_genus", "difficult_to_genus",
    "genus", "family"
  ),
  unresolved_rate = 0,
  rarity_exponent = 0,
  difficulty_multiplier = 1,
  rare_quantile = 0.10,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  workflow <- match.arg(workflow)

  if (workflow == "species") {
    out <- comm[, colSums(comm) > 0, drop = FALSE]
    colnames(out) <- paste0("species::", colnames(out))
    return(list(comm = out, workflow_summary = tibble(workflow = workflow, realised_unresolved_rate = 0)))
  }

  if (workflow == "genus") {
    return(list(
      comm = aggregate_by_taxonomic_rank(comm, taxonomy, rank = "genus"),
      workflow_summary = tibble(workflow = workflow, realised_unresolved_rate = NA_real_)
    ))
  }

  if (workflow == "family") {
    return(list(
      comm = aggregate_by_taxonomic_rank(comm, taxonomy, rank = "family"),
      workflow_summary = tibble(workflow = workflow, realised_unresolved_rate = NA_real_)
    ))
  }

  if (workflow == "mixed_rtu" || workflow == "drop_unresolved") {
    mixed <- report_mixed_rtu(
      comm = comm,
      taxonomy = taxonomy,
      unresolved_rate = unresolved_rate,
      rarity_exponent = rarity_exponent,
      difficulty_multiplier = difficulty_multiplier,
      seed = seed
    )

    out <- mixed$comm
    if (workflow == "drop_unresolved") {
      out <- out[, grepl("^species::", colnames(out)), drop = FALSE]
      out <- out[, colSums(out) > 0, drop = FALSE]
    }

    return(list(
      comm = out,
      workflow_summary = tibble(
        workflow = workflow,
        realised_unresolved_rate = mixed$realised_unresolved_rate
      )
    ))
  }

  active <- active_species(comm)
  totals <- colSums(comm)[active]
  tax_active <- taxonomy %>% filter(species %in% active)

  p_unresolved <- stats::setNames(rep(0, length(active)), active)

  if (workflow == "rare_to_genus") {
    assert_scalar(rare_quantile, "rare_quantile", lower = 0, upper = 1)
    threshold <- stats::quantile(totals, probs = rare_quantile, na.rm = TRUE, names = FALSE)
    p_unresolved[names(totals)[totals <= threshold]] <- 1
  }

  if (workflow == "difficult_to_genus") {
    difficult_species <- tax_active %>% filter(difficult_genus) %>% pull(species)
    p_unresolved[difficult_species] <- 1
  }

  # Build a deterministic mixed RTU matrix with p = 0 or 1.
  species_labels <- paste0("species::", active)
  genus_labels <- paste0("genus::", unique(tax_active$genus))
  out <- matrix(
    0,
    nrow = nrow(comm),
    ncol = length(species_labels) + length(genus_labels),
    dimnames = list(rownames(comm), c(species_labels, genus_labels))
  )

  for (species in active) {
    genus <- tax_active$genus[match(species, tax_active$species)]
    counts <- comm[, species]
    unresolved <- counts * p_unresolved[[species]]
    out[, paste0("species::", species)] <- counts - unresolved
    out[, paste0("genus::", genus)] <- out[, paste0("genus::", genus)] + unresolved
  }

  out <- out[, colSums(out) > 0, drop = FALSE]

  list(
    comm = out,
    workflow_summary = tibble(
      workflow = workflow,
      realised_unresolved_rate = sum(p_unresolved * totals) / sum(totals)
    )
  )
}

apply_observation_process <- function(
  comm_true,
  taxonomy,
  process = list(),
  seed = NULL
) {
  defaults <- list(
    id_error_rate = 0,
    candidate_pool = "observed",
    id_rarity_exponent = 0,
    id_difficulty_multiplier = 1,
    id_max_error = 0.95,
    workflow = "species",
    unresolved_rate = 0,
    unresolved_rarity_exponent = 0,
    unresolved_difficulty_multiplier = 1,
    rare_quantile = 0.10
  )
  process <- utils::modifyList(defaults, process)

  error_out <- apply_identification_error(
    comm = comm_true,
    taxonomy = taxonomy,
    target_mean_error = process$id_error_rate,
    candidate_pool = process$candidate_pool,
    rarity_exponent = process$id_rarity_exponent,
    difficulty_multiplier = process$id_difficulty_multiplier,
    max_error = process$id_max_error,
    seed = seed
  )

  workflow_out <- apply_reporting_workflow(
    comm = error_out$comm,
    taxonomy = taxonomy,
    workflow = process$workflow,
    unresolved_rate = process$unresolved_rate,
    rarity_exponent = process$unresolved_rarity_exponent,
    difficulty_multiplier = process$unresolved_difficulty_multiplier,
    rare_quantile = process$rare_quantile,
    seed = if (is.null(seed)) NULL else seed + 100000L
  )

  list(
    comm = workflow_out$comm,
    error_summary = error_out$error_summary,
    workflow_summary = workflow_out$workflow_summary,
    process = process
  )
}
