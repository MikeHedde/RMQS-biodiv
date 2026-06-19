# =============================================================================
# 05_build_multitaxa_taxref_pools.R
# Build France-mainland TAXREF candidate pools and expert-confusion templates
# for the RMQS 2024 multi-taxon taxonomic-uncertainty analysis.
#
# INPUTS
#   1. outputs_multitaxa_2024/assemblages/*__taxon_lookup.csv
#   2. A full TAXREF export (semicolon- or tab-delimited) supplied locally.
#
# OUTPUTS
#   regional_pools/<assemblage_id>__regional_pool.csv
#   regional_pools/<assemblage_id>__observed_taxref_match_audit.csv
#   expert_confusion_maps/<assemblage_id>__expert_confusions_TEMPLATE.csv
#   regional_pools/regional_pool_build_summary.csv
#
# Design:
#   - TAXREF is used to identify accepted, mainland-France candidate species.
#   - Observed species remain represented by their original RMQS taxon_unit.
#   - Unobserved TAXREF candidates use "species:cdref:<CD_REF>".
#   - This avoids creating duplicate IDs when an observed CD_NOM is a synonym
#     of a TAXREF accepted CD_REF.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(janitor)
  library(purrr)
})

# ---- 0. User settings --------------------------------------------------------
INPUT_DIR <- "outputs_multitaxa_2024"

# Point this to the complete TAXREF export downloaded from the INPN/MNHN.
# The script supports common TAXREF column names and semicolon/tab delimiters.
TAXREF_FILE <- "TAXREF_full.csv"

REGIONAL_POOL_DIR <- "regional_pools"
EXPERT_MAP_DIR <- "expert_confusion_maps"

# Mainland France statuses retained in the published Collembola workflow.
FR_STATUS_KEEP <- c("P", "E", "C")

# Species ranks retained. TAXREF exports generally use ES; S is accepted too.
SPECIES_RANKS <- c("ES", "S")

# Assemblage-specific TAXREF filters.
# Each rule is applied case-insensitively after cleaning taxonomic names.
ASSEMBLAGE_FILTERS <- tibble::tribble(
  ~assemblage_id,              ~rank_column, ~rank_value,
  "collembola_soil_core",      "ordre",      "Collembola",
  "araneae_pitfall",           "ordre",      "Araneae",
  "carabidae_pitfall",         "famille",    "Carabidae",
  "formicidae_pitfall",        "famille",    "Formicidae",
  "isopoda_pitfall",           "ordre",      "Isopoda",
  "isopoda_hand_sorting",      "ordre",      "Isopoda",
  "diplopoda_pitfall",         "classe",     "Diplopoda",
  "diplopoda_hand_sorting",    "classe",     "Diplopoda"
)

# ---- 1. Helpers --------------------------------------------------------------
clean_chr <- function(x) {
  x %>%
    as.character() %>%
    str_squish() %>%
    na_if("")
}

norm_taxon <- function(x) {
  x %>%
    clean_chr() %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    str_to_lower()
}

extract_binomial <- function(x) {
  x <- clean_chr(x)
  out <- str_extract(x, "^[A-Z][a-zA-Z-]+\\s+[a-z][a-zA-Z-]+")
  clean_chr(out)
}

extract_genus <- function(x) {
  x <- clean_chr(x)
  clean_chr(str_extract(x, "^[A-Z][a-zA-Z-]+"))
}

read_taxref_flexible <- function(path) {
  if (!file.exists(path)) {
    stop(
      "TAXREF file not found: ", path, "\n",
      "Set TAXREF_FILE at the top of this script."
    )
  }
  
  first_line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  delim <- if (str_count(first_line, ";") >= str_count(first_line, "\t")) ";" else "\t"
  
  readr::read_delim(
    path,
    delim = delim,
    locale = readr::locale(encoding = "Latin1"),
    show_col_types = FALSE,
    name_repair = "unique"
  ) %>%
    janitor::clean_names()
}

pick_col <- function(dat, candidates, required = TRUE) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0L) {
    if (required) stop("Missing required TAXREF column. Tried: ", paste(candidates, collapse = ", "))
    return(rep(NA_character_, nrow(dat)))
  }
  dat[[hit[[1]]]]
}

make_taxref_standard <- function(taxref_raw) {
  out <- tibble(
    cd_nom = clean_chr(pick_col(taxref_raw, c("cd_nom"))),
    cd_ref = clean_chr(pick_col(taxref_raw, c("cd_ref"))),
    rang = toupper(clean_chr(pick_col(taxref_raw, c("rang", "rank")))),
    lb_nom = clean_chr(pick_col(taxref_raw, c("lb_nom", "nom_complet"), required = FALSE)),
    nom_valide = clean_chr(pick_col(taxref_raw, c("nom_valide"), required = FALSE)),
    regne = clean_chr(pick_col(taxref_raw, c("regne"), required = FALSE)),
    phylum = clean_chr(pick_col(taxref_raw, c("phylum"), required = FALSE)),
    classe = clean_chr(pick_col(taxref_raw, c("classe", "class"), required = FALSE)),
    ordre = clean_chr(pick_col(taxref_raw, c("ordre", "order"), required = FALSE)),
    famille = clean_chr(pick_col(taxref_raw, c("famille", "family"), required = FALSE)),
    fr = clean_chr(pick_col(taxref_raw, c("fr"), required = FALSE))
  ) %>%
    mutate(
      accepted_name = coalesce(extract_binomial(nom_valide), extract_binomial(lb_nom)),
      genus = extract_genus(accepted_name),
      cd_ref = coalesce(cd_ref, cd_nom),
      is_accepted_row = is.na(cd_nom) | is.na(cd_ref) | cd_nom == cd_ref
    ) %>%
    filter(
      rang %in% SPECIES_RANKS,
      fr %in% FR_STATUS_KEEP,
      !is.na(cd_ref),
      !is.na(accepted_name),
      !is.na(genus)
    )
  
  # Keep accepted CD_REF rows where possible. In case an export includes only a
  # synonym row for a valid CD_REF, retain one representative deterministically.
  out %>%
    arrange(desc(is_accepted_row)) %>%
    distinct(cd_ref, .keep_all = TRUE)
}

read_lookup_species <- function(assemblage_id) {
  path <- file.path(INPUT_DIR, "assemblages", paste0(assemblage_id, "__taxon_lookup.csv"))
  if (!file.exists(path)) {
    stop("Missing taxon lookup for ", assemblage_id, ": ", path)
  }
  
  readr::read_csv(path, show_col_types = FALSE) %>%
    transmute(
      observed_taxon_unit = taxon_unit_species,
      observed_genus = genus,
      observed_family = family,
      observed_label = coalesce(binomial, valid_name, name)
    ) %>%
    filter(!is.na(observed_taxon_unit)) %>%
    mutate(
      observed_cd_nom = str_match(observed_taxon_unit, "^species:cdnom:(.+)$")[, 2],
      observed_label = clean_chr(observed_label),
      observed_genus = clean_chr(observed_genus),
      observed_family = clean_chr(observed_family)
    ) %>%
    distinct()
}

filter_taxref_for_assemblage <- function(taxref_std, assemblage_id) {
  rule <- ASSEMBLAGE_FILTERS %>% filter(.data$assemblage_id == !!assemblage_id)
  if (nrow(rule) != 1L) stop("No unique TAXREF filter rule for: ", assemblage_id)
  
  col <- rule$rank_column[[1]]
  val <- norm_taxon(rule$rank_value[[1]])
  
  if (!col %in% names(taxref_std)) {
    stop("TAXREF standardisation does not contain filter column: ", col)
  }
  
  taxref_std %>%
    filter(norm_taxon(.data[[col]]) == val)
}

build_regional_pool_one <- function(assemblage_id, taxref_std) {
  observed <- read_lookup_species(assemblage_id)
  taxref_group <- filter_taxref_for_assemblage(taxref_std, assemblage_id)
  
  # First map observed CD_NOM to accepted CD_REF. Fallback on exact binomial
  # matching only when a CD_NOM is absent or unmatched.
  observed_match <- observed %>%
    left_join(
      taxref_group %>%
        select(
          cd_nom, cd_ref, genus_taxref = genus, family_taxref = famille,
          accepted_name
        ),
      by = c("observed_cd_nom" = "cd_nom")
    ) %>%
    mutate(
      exact_name_key = norm_taxon(observed_label)
    ) %>%
    left_join(
      taxref_group %>%
        transmute(
          exact_name_key = norm_taxon(accepted_name),
          cd_ref_name = cd_ref,
          genus_name = genus,
          family_name = famille,
          accepted_name_name = accepted_name
        ),
      by = "exact_name_key"
    ) %>%
    mutate(
      cd_ref = coalesce(cd_ref, cd_ref_name),
      genus = coalesce(genus_taxref, genus_name, observed_genus),
      family = coalesce(family_taxref, family_name, observed_family),
      accepted_name = coalesce(accepted_name, accepted_name_name, observed_label),
      match_method = case_when(
        !is.na(observed_cd_nom) & !is.na(cd_ref) ~ "cd_nom",
        is.na(observed_cd_nom) & !is.na(cd_ref) ~ "accepted_name",
        is.na(cd_ref) ~ "unmatched"
      )
    ) %>%
    select(
      observed_taxon_unit, observed_cd_nom, observed_label,
      observed_genus, observed_family,
      cd_ref, genus, family, accepted_name, match_method
    )
  
  # One regional candidate per accepted CD_REF. If it is already observed under
  # any synonymous CD_NOM, preserve the existing RMQS taxon unit instead of
  # generating a duplicate cdref-based unit.
  observed_by_ref <- observed_match %>%
    filter(!is.na(cd_ref)) %>%
    group_by(cd_ref) %>%
    summarise(
      observed_taxon_unit = first(observed_taxon_unit),
      observed_label = first(observed_label),
      .groups = "drop"
    )
  
  pool <- taxref_group %>%
    left_join(observed_by_ref, by = "cd_ref") %>%
    transmute(
      taxon_unit = coalesce(observed_taxon_unit, paste0("species:cdref:", cd_ref)),
      cd_ref,
      genus,
      family,
      accepted_name,
      observed_in_assemblage = !is.na(observed_taxon_unit),
      candidate_origin = if_else(observed_in_assemblage, "observed_RMQS", "TAXREF_mainland")
    ) %>%
    distinct(taxon_unit, .keep_all = TRUE) %>%
    arrange(genus, accepted_name)
  
  # Candidate counts are useful for interpreting the size of the regional
  # confusion pool per observed species.
  source_audit <- observed_match %>%
    left_join(
      pool %>%
        count(genus, name = "n_regional_candidates"),
      by = "genus"
    ) %>%
    mutate(
      n_alternative_candidates = pmax(coalesce(n_regional_candidates, 0L) - 1L, 0L),
      in_regional_pool = !is.na(cd_ref),
      assemblage_id = assemblage_id
    )
  
  list(pool = pool, audit = source_audit)
}

make_expert_template <- function(assemblage_id, regional_pool, source_audit) {
  source <- source_audit %>%
    filter(!is.na(observed_taxon_unit), !is.na(genus)) %>%
    transmute(
      source_taxon_unit = observed_taxon_unit,
      source_label = observed_label,
      source_genus = genus
    ) %>%
    distinct()
  
  source %>%
    inner_join(
      regional_pool %>%
        transmute(
          target_taxon_unit = taxon_unit,
          target_label = accepted_name,
          target_genus = genus,
          target_origin = candidate_origin
        ),
      by = c("source_genus" = "target_genus")
    ) %>%
    filter(source_taxon_unit != target_taxon_unit) %>%
    transmute(
      source_taxon_unit,
      source_label,
      source_genus,
      target_taxon_unit,
      target_label,
      target_genus,
      target_origin,
      weight = 1,
      enabled = FALSE,
      comment = NA_character_
    ) %>%
    arrange(source_genus, source_label, target_label)
}

# ---- 2. Run ------------------------------------------------------------------
dir.create(REGIONAL_POOL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EXPERT_MAP_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(file.path(INPUT_DIR, "assemblages"))) {
  stop("Missing prepared assemblies. Run 00_prepare_multitaxa_inputs_2024.R first.")
}

taxref_raw <- read_taxref_flexible(TAXREF_FILE)
taxref_std <- make_taxref_standard(taxref_raw)

message("TAXREF candidates retained after species/mainland filters: ", nrow(taxref_std))

summary_rows <- list()

for (assemblage_id in ASSEMBLAGE_FILTERS$assemblage_id) {
  message("Building regional pool: ", assemblage_id)
  
  out <- build_regional_pool_one(assemblage_id, taxref_std)
  
  pool_path <- file.path(
    REGIONAL_POOL_DIR,
    paste0(assemblage_id, "__regional_pool.csv")
  )
  audit_path <- file.path(
    REGIONAL_POOL_DIR,
    paste0(assemblage_id, "__observed_taxref_match_audit.csv")
  )
  template_path <- file.path(
    EXPERT_MAP_DIR,
    paste0(assemblage_id, "__expert_confusions_TEMPLATE.csv")
  )
  
  template <- make_expert_template(assemblage_id, out$pool, out$audit)
  
  readr::write_csv(out$pool, pool_path)
  readr::write_csv(out$audit, audit_path)
  readr::write_csv(template, template_path)
  
  summary_rows[[assemblage_id]] <- tibble(
    assemblage_id = assemblage_id,
    n_regional_species = nrow(out$pool),
    n_observed_species = n_distinct(out$audit$observed_taxon_unit),
    n_observed_matched_to_taxref = sum(out$audit$in_regional_pool, na.rm = TRUE),
    match_rate = mean(out$audit$in_regional_pool, na.rm = TRUE),
    n_observed_species_with_alternative_candidate =
      sum(out$audit$n_alternative_candidates > 0, na.rm = TRUE),
    median_n_alternative_candidates =
      median(out$audit$n_alternative_candidates, na.rm = TRUE),
    n_expert_template_rows = nrow(template)
  )
}

summary_tbl <- bind_rows(summary_rows)
readr::write_csv(
  summary_tbl,
  file.path(REGIONAL_POOL_DIR, "regional_pool_build_summary.csv")
)

message("\nCompleted.")
message("Regional pools: ", normalizePath(REGIONAL_POOL_DIR))
message("Expert templates: ", normalizePath(EXPERT_MAP_DIR))
message("\nBefore the main v2 analysis:")
message("  1. Inspect regional_pool_build_summary.csv and each match audit.")
message("  2. Copy any expert template to a non-TEMPLATE filename, e.g.")
message("     expert_confusion_maps/collembola_soil_core__expert_confusions.csv")
message("  3. Set enabled = TRUE only for biologically justified source-target pairs.")
