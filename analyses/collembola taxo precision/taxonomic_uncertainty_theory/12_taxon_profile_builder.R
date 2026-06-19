# ============================================================
# 12_taxon_profile_builder.R
# Standardised empirical taxon profiles for cross-taxon comparison
# ============================================================
#
# Purpose
# -------
# Build a descriptive profile for a biological group from two standardised files:
#   1) a site x species community table in long format;
#   2) a regional species -> genus -> family taxonomy.
#
# The profile is NOT a calibration target. It only describes the taxonomic
# architecture and sampling context used to position an empirical group
# relative to the theoretical simulations.

safe_cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || mean(x) == 0) return(NA_real_)
  stats::sd(x) / mean(x)
}

assert_standard_taxonomy <- function(taxonomy, name = "taxonomy") {
  required <- c("species", "genus", "family")
  if (!is.data.frame(taxonomy)) {
    stop(sprintf("`%s` must be a data frame.", name), call. = FALSE)
  }
  missing <- setdiff(required, names(taxonomy))
  if (length(missing) > 0) {
    stop(sprintf("`%s` lacks columns: %s", name, paste(missing, collapse = ", ")), call. = FALSE)
  }

  out <- taxonomy %>%
    transmute(
      species = as.character(species),
      genus = as.character(genus),
      family = as.character(family),
      difficult_genus = if ("difficult_genus" %in% names(taxonomy)) as.logical(difficult_genus) else FALSE
    ) %>%
    filter(!is.na(species), !is.na(genus), !is.na(family), species != "", genus != "", family != "") %>%
    distinct(species, .keep_all = TRUE)

  if (nrow(out) == 0) stop(sprintf("`%s` contains no valid species.", name), call. = FALSE)
  out
}

assert_standard_community_long <- function(community_long, name = "community_long") {
  required <- c("site", "species", "abundance")
  if (!is.data.frame(community_long)) {
    stop(sprintf("`%s` must be a data frame.", name), call. = FALSE)
  }
  missing <- setdiff(required, names(community_long))
  if (length(missing) > 0) {
    stop(sprintf("`%s` lacks columns: %s", name, paste(missing, collapse = ", ")), call. = FALSE)
  }

  out <- community_long %>%
    transmute(
      site = as.character(site),
      species = as.character(species),
      abundance = suppressWarnings(as.numeric(abundance))
    ) %>%
    filter(!is.na(site), !is.na(species), !is.na(abundance), abundance > 0, site != "", species != "")

  if (nrow(out) == 0) stop(sprintf("`%s` contains no positive records.", name), call. = FALSE)
  out
}

community_long_to_matrix <- function(community_long) {
  community_long <- assert_standard_community_long(community_long)

  community_long %>%
    group_by(site, species) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = species,
      values_from = abundance,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("site") %>%
    as.matrix()
}

read_standard_taxon_inputs <- function(community_file, taxonomy_file, community_delim = ",", taxonomy_delim = ",") {
  community_long <- readr::read_delim(
    community_file, delim = community_delim, show_col_types = FALSE, progress = FALSE
  )
  regional_taxonomy <- readr::read_delim(
    taxonomy_file, delim = taxonomy_delim, show_col_types = FALSE, progress = FALSE
  )

  list(
    community_long = assert_standard_community_long(community_long),
    regional_taxonomy = assert_standard_taxonomy(regional_taxonomy)
  )
}

summarise_confusion_capacity <- function(regional_taxonomy, observed_species) {
  regional_taxonomy <- assert_standard_taxonomy(regional_taxonomy)
  observed_species <- unique(as.character(observed_species))

  regional_genus_sizes <- regional_taxonomy %>%
    count(genus, name = "regional_species_per_genus")

  observed_taxonomy <- regional_taxonomy %>%
    filter(species %in% observed_species)

  observed_genus_sizes <- observed_taxonomy %>%
    count(genus, name = "observed_species_per_genus")

  by_observed_species <- observed_taxonomy %>%
    left_join(regional_genus_sizes, by = "genus") %>%
    left_join(observed_genus_sizes, by = "genus") %>%
    mutate(
      regional_congener_candidates = pmax(regional_species_per_genus - 1L, 0L),
      observed_congener_candidates = pmax(observed_species_per_genus - 1L, 0L)
    )

  tibble(
    n_observed_species_with_regional_congeners = sum(by_observed_species$regional_congener_candidates > 0),
    prop_observed_species_with_regional_congeners = mean(by_observed_species$regional_congener_candidates > 0),
    mean_regional_congener_candidates = mean(by_observed_species$regional_congener_candidates),
    median_regional_congener_candidates = stats::median(by_observed_species$regional_congener_candidates),
    n_observed_species_with_observed_congeners = sum(by_observed_species$observed_congener_candidates > 0),
    prop_observed_species_with_observed_congeners = mean(by_observed_species$observed_congener_candidates > 0),
    mean_observed_congener_candidates = mean(by_observed_species$observed_congener_candidates),
    median_observed_congener_candidates = stats::median(by_observed_species$observed_congener_candidates)
  )
}

build_empirical_taxon_profile <- function(taxon, community_long, regional_taxonomy) {
  if (length(taxon) != 1 || is.na(taxon) || taxon == "") {
    stop("`taxon` must be one non-empty string.", call. = FALSE)
  }

  community_long <- assert_standard_community_long(community_long)
  regional_taxonomy <- assert_standard_taxonomy(regional_taxonomy)

  unmatched <- setdiff(unique(community_long$species), regional_taxonomy$species)
  if (length(unmatched) > 0) {
    warning(
      sprintf(
        "%s community records contain %d species absent from regional taxonomy; these records are excluded from the profile.",
        taxon, length(unmatched)
      ),
      call. = FALSE
    )
    community_long <- community_long %>% filter(species %in% regional_taxonomy$species)
  }

  comm <- community_long_to_matrix(community_long)
  observed_species <- colnames(comm)[colSums(comm) > 0]
  observed_taxonomy <- regional_taxonomy %>% filter(species %in% observed_species)

  if (length(observed_species) == 0) stop("No observed species remain after taxonomy matching.", call. = FALSE)

  regional_architecture <- summarise_taxonomic_architecture(regional_taxonomy)
  observed_architecture <- summarise_taxonomic_architecture(observed_taxonomy)
  comm_metrics <- community_metrics(comm)
  confusion <- summarise_confusion_capacity(regional_taxonomy, observed_species)

  summary <- bind_cols(
    tibble(taxon = taxon),
    regional_architecture %>% rename_with(~ paste0("regional_", .x)),
    observed_architecture %>% rename_with(~ paste0("observed_", .x)),
    tibble(
      n_sites = nrow(comm),
      n_observed_species_from_community = length(observed_species),
      observed_alpha_q0 = comm_metrics$alpha_q0,
      observed_alpha_q1 = comm_metrics$alpha_q1,
      observed_alpha_q2 = comm_metrics$alpha_q2,
      observed_gamma = comm_metrics$gamma,
      observed_mean_occupancy = comm_metrics$mean_occupancy,
      observed_mean_bray_curtis = comm_metrics$mean_bray_curtis,
      regional_observed_species_ratio = regional_architecture$n_regional_species / length(observed_species)
    ),
    confusion
  )

  list(
    summary = summary,
    community_matrix = comm,
    community_long = community_long,
    regional_taxonomy = regional_taxonomy,
    observed_taxonomy = observed_taxonomy
  )
}

write_empirical_taxon_profile <- function(profile, output_dir = THEORY_OUT_DIR) {
  if (!is.list(profile) || is.null(profile$summary)) {
    stop("`profile` must be an object returned by build_empirical_taxon_profile().", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  taxon_key <- gsub("[^A-Za-z0-9]+", "_", profile$summary$taxon[[1]])

  readr::write_csv(profile$summary, file.path(output_dir, paste0("profile_", taxon_key, "_summary.csv")))
  readr::write_csv(profile$regional_taxonomy, file.path(output_dir, paste0("profile_", taxon_key, "_regional_taxonomy.csv")))
  readr::write_csv(profile$observed_taxonomy, file.path(output_dir, paste0("profile_", taxon_key, "_observed_taxonomy.csv")))
  readr::write_csv(profile$community_long, file.path(output_dir, paste0("profile_", taxon_key, "_community_long.csv")))
  invisible(profile)
}

profile_axes_long <- function(profile_summaries) {
  required <- c(
    "taxon", "regional_observed_species_ratio", "regional_mean_species_per_genus",
    "regional_mean_genera_per_family", "observed_mean_occupancy"
  )
  missing <- setdiff(required, names(profile_summaries))
  if (length(missing) > 0) {
    stop("Profile summary table lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  profile_summaries %>%
    select(all_of(required)) %>%
    pivot_longer(-taxon, names_to = "axis", values_to = "value") %>%
    mutate(
      axis = recode(
        axis,
        regional_observed_species_ratio = "Regional / observed species ratio",
        regional_mean_species_per_genus = "Mean species per genus",
        regional_mean_genera_per_family = "Mean genera per family",
        observed_mean_occupancy = "Observed mean occupancy"
      )
    )
}
