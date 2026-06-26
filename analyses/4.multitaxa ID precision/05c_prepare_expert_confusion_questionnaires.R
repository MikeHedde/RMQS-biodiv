# =============================================================================
# 05c_prepare_expert_confusion_questionnaires.R
# -----------------------------------------------------------------------------
# Prepare taxonomist-facing CSV questionnaires to curate expert confusion maps
# for the RMQS 2024 multi-taxon taxonomic-uncertainty analysis.
#
# PURPOSE
# -------
# This script does NOT simulate anything. It creates one review sheet per
# assemblage × protocol from:
#   1. observed RMQS species and their abundance / occurrence;
#   2. the accepted mainland-France TAXREF pool built by script 05;
#   3. all same-genus source -> target candidate pairs.
#
# A taxonomist reviews only biologically plausible directional confusions:
#   "Could a true individual of SOURCE be reported as TARGET under the
#    RMQS identification workflow?"
#
# The returned, completed sheets can then be compiled into the final files:
#   expert_confusion_maps/<assemblage_id>__expert_confusions.csv
# required by:
#   01_run_multitaxa_uncertainty_2024_v3_taxref_coherent_confusion.R
#
# IMPORTANT
# ---------
# - Same-genus candidates are intentionally over-inclusive. The taxonomist
#   should mark most impossible pairs as "no".
# - Direction matters. SOURCE -> TARGET need not imply TARGET -> SOURCE.
# - The aim is NOT to recreate a formal taxonomic key; it is to encode
#   realistic identification-confusion pathways under routine monitoring.
# - Cross-genus confusions can be added in the separate "additional_pairs" CSV.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# ---- 0. User settings --------------------------------------------------------

# Inputs created by 00_prepare_multitaxa_inputs_2024.R and
# 05_build_multitaxa_taxref_pools_v2.R.
INPUT_DIR <- "outputs_multitaxa_2024"
REGIONAL_POOL_DIR <- "regional_pools"

# Folder that you can send directly to taxonomists.
OUT_DIR <- "expert_review_to_send"

# Folder where returned, completed CSVs should be placed before switching
# MODE to "compile_reviewed".
RETURNED_DIR <- "expert_review_returned"

# "create_template": produce review sheets for taxonomists.
# "compile_reviewed": transform completed sheets into the expert confusion
#                     maps used by the main analysis.
MODE <- "create_template"

# Generate questionnaires only for the main assemblages by default.
# Add secondary assemblages if you have an available expert and a reason to
# review their identification uncertainty.
ASSEMBLAGES_TO_REVIEW <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "carabidae_pitfall",
  "formicidae_pitfall",
  "isopoda_pitfall"
)

# Candidate pairs are generated within genus. Do NOT change this to a broader
# taxonomic rank: this sheet is already intentionally inclusive.
INCLUDE_SAME_GENUS_CANDIDATES <- TRUE

# Keep every observed source species by default, including rare taxa. Set to 2
# or 3 only if a taxonomist explicitly asks for a shorter first-pass review.
MIN_SOURCE_TOTAL_ABUNDANCE <- 1L

# Optional hard cap for enormous genera. Inf means no cap. When a cap is used,
# the script keeps the most ecologically relevant targets first:
# (i) observed RMQS targets, then (ii) regional TAXREF targets.
MAX_CANDIDATES_PER_SOURCE <- Inf

# ---- 1. Labels and controlled vocabularies ----------------------------------

ASSEMBLAGE_LABELS <- c(
  collembola_soil_core = "Collembola — carottes de sol",
  araneae_pitfall = "Araneae — pièges Barber",
  carabidae_pitfall = "Carabidae — pièges Barber",
  formicidae_pitfall = "Formicidae — pièges Barber",
  isopoda_pitfall = "Isopoda — pièges Barber",
  isopoda_hand_sorting = "Isopoda — tri manuel",
  diplopoda_pitfall = "Diplopoda — pièges Barber",
  diplopoda_hand_sorting = "Diplopoda — tri manuel"
)

# These are written into each CSV header / README for human use.
REVIEW_STATUS_CHOICES <- c(
  "yes_likely",   # a likely or recurrent confusion
  "yes_possible", # possible, but uncommon / conditional
  "uncertain",    # taxonomist cannot decide from current information
  "no"            # not a plausible confusion
)

CONTEXT_CHOICES <- c(
  "adult",
  "juvenile_or_larva",
  "sex_specific",
  "damaged_or_partial_specimen",
  "routine_morphological_variation",
  "not_stage_dependent",
  "other"
)

REPORTING_CHOICES <- c(
  "species_reliable",
  "species_with_caution",
  "report_to_genus",
  "report_to_family",
  "not_applicable"
)

# ---- 2. Generic helpers ------------------------------------------------------

required_file <- function(path, what) {
  if (!file.exists(path)) {
    stop(
      "Missing ", what, ":\n", path, "\n\n",
      "For regional pools, run 05_build_multitaxa_taxref_pools_v2.R first."
    )
  }
  path
}

clean_chr <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_squish() %>%
    na_if("")
}

write_excel_csv_fr <- function(x, path) {
  # UTF-8 BOM + semicolon delimiter: convenient for French Excel installations.
  readr::write_excel_csv2(x, path, na = "")
}

source_stats_one <- function(assemblage_id) {
  species_path <- required_file(
    file.path(INPUT_DIR, "assemblages", paste0(assemblage_id, "__species_long.csv")),
    paste0("species-long input for ", assemblage_id)
  )
  
  readr::read_csv(species_path, show_col_types = FALSE) %>%
    mutate(
      station = as.character(station),
      taxon_unit = as.character(taxon_unit)
    ) %>%
    group_by(taxon_unit) %>%
    summarise(
      source_total_abundance = sum(abundance, na.rm = TRUE),
      source_n_stations = n_distinct(station[abundance > 0]),
      .groups = "drop"
    )
}

read_match_audit_one <- function(assemblage_id) {
  audit_path <- required_file(
    file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__observed_taxref_match_audit.csv")),
    paste0("TAXREF match audit for ", assemblage_id)
  )
  
  audit <- readr::read_csv(audit_path, show_col_types = FALSE) %>%
    janitor::clean_names()
  
  required <- c("observed_taxon_unit", "observed_label", "genus")
  missing <- setdiff(required, names(audit))
  if (length(missing)) {
    stop(
      "The TAXREF audit for ", assemblage_id,
      " does not contain: ", paste(missing, collapse = ", ")
    )
  }
  
  if (!"family" %in% names(audit)) audit$family <- NA_character_
  if (!"accepted_name" %in% names(audit)) audit$accepted_name <- NA_character_
  if (!"in_regional_pool" %in% names(audit)) audit$in_regional_pool <- NA
  
  audit %>%
    transmute(
      source_taxon_unit = as.character(observed_taxon_unit),
      source_label_rmq = clean_chr(observed_label),
      source_accepted_name = clean_chr(accepted_name),
      source_genus = clean_chr(genus),
      source_family = clean_chr(family),
      source_matched_to_taxref = as.logical(in_regional_pool)
    ) %>%
    distinct()
}

read_regional_pool_one <- function(assemblage_id) {
  pool_path <- required_file(
    file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__regional_pool.csv")),
    paste0("regional TAXREF pool for ", assemblage_id)
  )
  
  pool <- readr::read_csv(pool_path, show_col_types = FALSE) %>%
    janitor::clean_names()
  
  required <- c("taxon_unit", "accepted_name", "genus")
  missing <- setdiff(required, names(pool))
  if (length(missing)) {
    stop(
      "The regional pool for ", assemblage_id,
      " does not contain: ", paste(missing, collapse = ", ")
    )
  }
  
  if (!"family" %in% names(pool)) pool$family <- NA_character_
  if (!"candidate_origin" %in% names(pool)) pool$candidate_origin <- NA_character_
  
  pool %>%
    transmute(
      target_taxon_unit = as.character(taxon_unit),
      target_accepted_name = clean_chr(accepted_name),
      target_genus = clean_chr(genus),
      target_family = clean_chr(family),
      target_origin = clean_chr(candidate_origin)
    ) %>%
    distinct()
}

candidate_priority <- function(
    source_total_abundance, source_n_stations,
    target_origin, target_total_abundance, target_n_stations
) {
  # A transparent triage score only. It does NOT determine which pairs are
  # scientifically valid; the taxonomist does that.
  2 * log1p(source_total_abundance) +
    0.5 * source_n_stations +
    ifelse(target_origin == "observed_RMQS", 2, 0) +
    0.25 * log1p(replace_na(target_total_abundance, 0)) +
    0.10 * replace_na(target_n_stations, 0)
}

priority_class <- function(x) {
  case_when(
    is.na(x) ~ "normal",
    x >= quantile(x, 0.80, na.rm = TRUE) ~ "high",
    x >= quantile(x, 0.50, na.rm = TRUE) ~ "medium",
    TRUE ~ "normal"
  )
}

write_readme <- function(path) {
  lines <- c(
    "RMQS 2024 — revue experte des confusions taxonomiques",
    "====================================================",
    "",
    "Objectif",
    "--------",
    "Pour chaque ligne, indiquer si un individu réellement appartenant à",
    "SOURCE pourrait être rapporté comme TARGET dans le cadre d'une",
    "identification de routine comparable au protocole RMQS.",
    "",
    "Les paires proposées sont volontairement larges : même genre et pool",
    "TAXREF France continentale. Elles ne sont PAS toutes censées être",
    "plausibles. Merci de renseigner 'no' pour les impossibilités.",
    "",
    "Colonnes à renseigner",
    "---------------------",
    "review_status :",
    "  - yes_likely   : confusion probable / déjà rencontrée",
    "  - yes_possible : confusion possible mais contextuelle ou rare",
    "  - uncertain    : impossible à trancher avec les éléments disponibles",
    "  - no           : confusion non plausible",
    "",
    "relative_weight :",
    "  - 3 = voie de confusion principale ou fréquente",
    "  - 2 = plausible, contribution secondaire",
    "  - 1 = possible mais rare",
    "  - 0 = non plausible (équivalent pratique à review_status = no)",
    "",
    "confusion_context : choisir parmi adult, juvenile_or_larva, sex_specific,",
    "damaged_or_partial_specimen, routine_morphological_variation,",
    "not_stage_dependent, other.",
    "",
    "recommended_reporting_level :",
    "  - species_reliable",
    "  - species_with_caution",
    "  - report_to_genus",
    "  - report_to_family",
    "  - not_applicable",
    "",
    "Important",
    "---------",
    "1. La direction compte : SOURCE -> TARGET ne suppose pas nécessairement",
    "   TARGET -> SOURCE.",
    "2. 'source_total_abundance' et 'source_n_stations' servent uniquement à",
    "   situer le poids potentiel de la confusion dans RMQS.",
    "3. Les cibles 'observed_RMQS' sont déjà observées dans l'assemblage ; les",
    "   cibles 'TAXREF_mainland' viennent du pool national, mais ne sont pas",
    "   nécessairement présentes localement.",
    "4. Pour une confusion inter-générique ou une paire absente du tableau,",
    "   utiliser le fichier *_additional_pairs.csv.",
    "",
    "Retour attendu",
    "--------------",
    "Conserver toutes les colonnes techniques, compléter les colonnes de revue,",
    "puis renvoyer le fichier au format CSV encodé UTF-8."
  )
  
  writeLines(lines, con = path, useBytes = TRUE)
}

# ---- 3. Create questionnaires ------------------------------------------------

create_questionnaire_one <- function(assemblage_id) {
  # Keep scalar identifiers outside dplyr data masks. Within mutate(),
  # `assemblage_id` may otherwise resolve to a data-frame column.
  assemblage_id_value <- as.character(assemblage_id)[1]
  assemblage_label_value <- unname(ASSEMBLAGE_LABELS[[assemblage_id_value]])
  
  if (is.null(assemblage_label_value) || is.na(assemblage_label_value)) {
    stop("No display label defined for assemblage_id: ", assemblage_id_value)
  }
  
  source_stats <- source_stats_one(assemblage_id_value)
  audit <- read_match_audit_one(assemblage_id_value)
  pool <- read_regional_pool_one(assemblage_id_value)
  
  observed_target_stats <- source_stats %>%
    transmute(
      target_taxon_unit = taxon_unit,
      target_total_abundance = source_total_abundance,
      target_n_stations = source_n_stations
    )
  
  sources <- audit %>%
    left_join(
      source_stats %>% transmute(
        source_taxon_unit = taxon_unit,
        source_total_abundance,
        source_n_stations
      ),
      by = "source_taxon_unit"
    ) %>%
    filter(
      !is.na(source_total_abundance),
      source_total_abundance >= MIN_SOURCE_TOTAL_ABUNDANCE,
      !is.na(source_genus)
    ) %>%
    distinct()
  
  if (!nrow(sources)) {
    warning("No source species retained for ", assemblage_id_value)
    return(list(pairs = tibble(), source_overview = tibble()))
  }
  
  pairs <- if (INCLUDE_SAME_GENUS_CANDIDATES) {
    sources %>%
      inner_join(
        pool,
        by = c("source_genus" = "target_genus"),
        relationship = "many-to-many"
      ) %>%
      # dplyr retains the left-hand join key (`source_genus`) and drops the
      # right-hand key (`target_genus`). Restore it explicitly for the review
      # table: source and target are intentionally congeneric at this stage.
      mutate(target_genus = source_genus) %>%
      filter(source_taxon_unit != target_taxon_unit) %>%
      left_join(observed_target_stats, by = "target_taxon_unit") %>%
      mutate(
        target_observed_in_rmq = target_origin == "observed_RMQS",
        priority_score = candidate_priority(
          source_total_abundance, source_n_stations, target_origin,
          target_total_abundance, target_n_stations
        )
      ) %>%
      group_by(source_taxon_unit) %>%
      arrange(desc(target_observed_in_rmq), desc(priority_score), target_accepted_name, .by_group = TRUE) %>%
      mutate(
        candidate_rank_for_source = row_number(),
        n_candidates_for_source = n()
      ) %>%
      ungroup()
  } else {
    tibble()
  }
  
  if (is.finite(MAX_CANDIDATES_PER_SOURCE)) {
    pairs <- pairs %>%
      filter(candidate_rank_for_source <= MAX_CANDIDATES_PER_SOURCE)
  }
  
  pairs <- pairs %>%
    mutate(
      review_priority = priority_class(priority_score),
      assemblage_id = assemblage_id_value,
      assemblage_label = assemblage_label_value,
      review_status = NA_character_,
      relative_weight = NA_real_,
      confusion_context = NA_character_,
      recommended_reporting_level = NA_character_,
      taxonomist_comment = NA_character_,
      reviewer_name = NA_character_,
      review_date = NA_character_
    ) %>%
    select(
      assemblage_id,
      assemblage_label,
      review_priority,
      candidate_rank_for_source,
      n_candidates_for_source,
      source_taxon_unit,
      source_label_rmq,
      source_accepted_name,
      source_genus,
      source_family,
      source_total_abundance,
      source_n_stations,
      source_matched_to_taxref,
      target_taxon_unit,
      target_accepted_name,
      target_genus,
      target_family,
      target_origin,
      target_observed_in_rmq,
      target_total_abundance,
      target_n_stations,
      review_status,
      relative_weight,
      confusion_context,
      recommended_reporting_level,
      taxonomist_comment,
      reviewer_name,
      review_date
    )
  
  source_overview <- sources %>%
    left_join(
      pairs %>%
        count(source_taxon_unit, name = "n_candidate_pairs_presented"),
      by = "source_taxon_unit"
    ) %>%
    mutate(
      n_candidate_pairs_presented = replace_na(n_candidate_pairs_presented, 0L),
      assemblage_id = assemblage_id_value,
      assemblage_label = assemblage_label_value,
      taxonomist_general_note = NA_character_
    ) %>%
    select(
      assemblage_id, assemblage_label,
      source_taxon_unit, source_label_rmq, source_accepted_name,
      source_genus, source_family,
      source_total_abundance, source_n_stations,
      source_matched_to_taxref,
      n_candidate_pairs_presented,
      taxonomist_general_note
    ) %>%
    arrange(desc(source_total_abundance), source_genus, source_accepted_name)
  
  list(pairs = pairs, source_overview = source_overview)
}

create_additional_pairs_template <- function(assemblage_id) {
  tibble(
    assemblage_id = assemblage_id,
    source_taxon_unit = NA_character_,
    source_label = NA_character_,
    target_taxon_unit = NA_character_,
    target_label = NA_character_,
    review_status = NA_character_,
    relative_weight = NA_real_,
    confusion_context = NA_character_,
    recommended_reporting_level = NA_character_,
    taxonomist_comment = NA_character_,
    reviewer_name = NA_character_,
    review_date = NA_character_
  )
}

# ---- 4. Compile returned reviews --------------------------------------------

read_completed_review <- function(path) {
  # Completed files are semicolon-delimited because templates use
  # write_excel_csv2(). We allow commas as a fallback.
  first_line <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (stringr::str_count(first_line, ";") >= stringr::str_count(first_line, ",")) {
    readr::read_csv2(path, show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
}

compile_one <- function(assemblage_id) {
  review_path <- file.path(
    RETURNED_DIR,
    paste0(assemblage_id, "__expert_confusion_review_COMPLETED.csv")
  )
  
  if (!file.exists(review_path)) {
    warning("No completed review found for ", assemblage_id, ": ", review_path)
    return(tibble())
  }
  
  reviewed <- read_completed_review(review_path) %>%
    janitor::clean_names()
  
  required <- c(
    "source_taxon_unit",
    "target_taxon_unit",
    "review_status",
    "relative_weight"
  )
  missing <- setdiff(required, names(reviewed))
  if (length(missing)) {
    stop(
      "Completed review for ", assemblage_id,
      " is missing: ", paste(missing, collapse = ", ")
    )
  }
  
  if (!"taxonomist_comment" %in% names(reviewed)) reviewed$taxonomist_comment <- NA_character_
  
  reviewed %>%
    mutate(
      review_status = tolower(clean_chr(review_status)),
      relative_weight = suppressWarnings(as.numeric(relative_weight)),
      relative_weight = replace_na(relative_weight, 0),
      enabled = review_status %in% c("yes_likely", "yes_possible") & relative_weight > 0
    ) %>%
    filter(enabled, source_taxon_unit != target_taxon_unit) %>%
    transmute(
      source_taxon_unit,
      target_taxon_unit,
      weight = relative_weight,
      enabled = TRUE,
      comment = taxonomist_comment
    ) %>%
    distinct()
}

# ---- 5. Execute --------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RETURNED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("expert_confusion_maps", recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(file.path(INPUT_DIR, "assemblages"))) {
  stop(
    "Prepared assemblages are missing. Run 00_prepare_multitaxa_inputs_2024.R first."
  )
}

if (MODE == "create_template") {
  readme_path <- file.path(OUT_DIR, "README_A_LIRE_POUR_LES_TAXONOMISTES.txt")
  write_readme(readme_path)
  
  summary_rows <- list()
  
  for (assemblage_id in ASSEMBLAGES_TO_REVIEW) {
    message("Preparing expert questionnaire: ", assemblage_id)
    
    result <- create_questionnaire_one(assemblage_id)
    
    pair_path <- file.path(
      OUT_DIR,
      paste0(assemblage_id, "__expert_confusion_review.csv")
    )
    overview_path <- file.path(
      OUT_DIR,
      paste0(assemblage_id, "__source_taxa_overview.csv")
    )
    additional_path <- file.path(
      OUT_DIR,
      paste0(assemblage_id, "__additional_pairs.csv")
    )
    
    write_excel_csv_fr(result$pairs, pair_path)
    write_excel_csv_fr(result$source_overview, overview_path)
    write_excel_csv_fr(create_additional_pairs_template(assemblage_id), additional_path)
    
    summary_rows[[assemblage_id]] <- tibble(
      assemblage_id = assemblage_id,
      assemblage_label = unname(ASSEMBLAGE_LABELS[[as.character(assemblage_id)[1]]]),
      n_source_species = nrow(result$source_overview),
      n_source_target_pairs = nrow(result$pairs),
      n_high_priority_pairs = sum(result$pairs$review_priority == "high", na.rm = TRUE),
      questionnaire_file = basename(pair_path),
      source_overview_file = basename(overview_path),
      additional_pairs_file = basename(additional_path)
    )
  }
  
  summary_tbl <- bind_rows(summary_rows)
  write_excel_csv_fr(
    summary_tbl,
    file.path(OUT_DIR, "questionnaire_manifest.csv")
  )
  
  message("\nQuestionnaires created in: ", normalizePath(OUT_DIR))
  message("Send each specialist:")
  message("  - README_A_LIRE_POUR_LES_TAXONOMISTES.txt")
  message("  - <assemblage>__expert_confusion_review.csv")
  message("  - <assemblage>__source_taxa_overview.csv")
  message("  - <assemblage>__additional_pairs.csv")
  message("\nThe taxonomist should rename the completed main sheet to:")
  message("  <assemblage>__expert_confusion_review_COMPLETED.csv")
  message("and return it in: ", normalizePath(RETURNED_DIR))
}

if (MODE == "compile_reviewed") {
  compiled <- purrr::map_dfr(ASSEMBLAGES_TO_REVIEW, compile_one)
  
  if (!nrow(compiled)) {
    stop(
      "No enabled expert-confusion pair was compiled. Check returned CSVs and ",
      "review_status / relative_weight entries."
    )
  }
  
  split_maps <- split(compiled, compiled$assemblage_id)
  
  for (assemblage_id in names(split_maps)) {
    final_map <- split_maps[[assemblage_id]] %>%
      select(source_taxon_unit, target_taxon_unit, weight, enabled, comment)
    
    final_path <- file.path(
      "expert_confusion_maps",
      paste0(assemblage_id, "__expert_confusions.csv")
    )
    
    readr::write_csv(final_map, final_path)
    message("Wrote final expert map: ", final_path)
  }
  
  readr::write_csv(
    compiled,
    file.path("expert_confusion_maps", "expert_confusion_maps_compilation_audit.csv")
  )
}
