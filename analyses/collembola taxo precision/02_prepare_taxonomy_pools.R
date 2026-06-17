# ============================================================
# 02_prepare_taxonomy_pools — TAXREF, espèces observées et pools de confusion
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances
load_required_step("01_import_clean")

STEP_ID <- "02_prepare_taxonomy_pools"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("02_prepare_taxonomy_pools — TAXREF, espèces observées et pools de confusion : cache chargé")
} else {
  # -----------------------------
  # 3. TAXREF MAINLAND + POOLS DE CONFUSION
  # -----------------------------

  message_header("Preparation TAXREF mainland et pools de confusion")

  taxref_mainland <- taxref_raw %>%
    mutate(
      cd_nom = clean_chr(cd_nom),
      cd_ref = clean_chr(cd_ref),
      rang = toupper(clean_chr(rang)),
      fr = clean_chr(fr),
      lb_nom = clean_chr(lb_nom),
      nom_valide = clean_chr(nom_valide),
      famille = clean_chr(famille),
      lb_binomial = extract_binomial(lb_nom),
      valid_binomial = extract_binomial(nom_valide),
      # Pour les candidats TAXREF, on travaille sur le taxon valide accepté.
      species_cd_ref = cd_ref,
      species_key = if_else(!is.na(species_cd_ref), paste0("cdref:", species_cd_ref), NA_character_),
      species_label = coalesce(valid_binomial, lb_binomial),
      genus = stringr::word(species_label, 1)
    ) %>%
    filter(
      rang %in% c("ES", "S"),
      fr %in% MAINLAND_FR_STATUS_KEEP,
      !is.na(species_key),
      !is.na(species_label),
      !is.na(genus)
    ) %>%
    distinct(species_key, species_cd_ref, species_label, genus, famille, fr)

  observed_species <- col_used %>%
    filter(taxo_level == "species") %>%
    group_by(species_key, species_cd_ref, species_label, binomial, genus, famille) %>%
    summarise(total_abundance = sum(abundance), n_rows = n(), .groups = "drop") %>%
    arrange(genus, species_label)

  rare_species_keys <- observed_species %>%
    filter(total_abundance <= RARE_MAX_TOTAL_ABUNDANCE) %>%
    pull(species_key)

  observed_same_genus_pool <- observed_species %>%
    transmute(
      candidate_key = species_key,
      candidate_label = species_label,
      genus,
      source_pool = "observed"
    ) %>%
    distinct()

  taxref_same_genus_pool <- taxref_mainland %>%
    transmute(
      candidate_key = species_key,
      candidate_label = species_label,
      genus,
      source_pool = "taxref_mainland"
    ) %>%
    distinct()

  confusion_pool_summary <- bind_rows(observed_same_genus_pool, taxref_same_genus_pool) %>%
    group_by(source_pool, genus) %>%
    summarise(n_candidate_species = n_distinct(candidate_key), .groups = "drop") %>%
    arrange(source_pool, genus)

  readr::write_csv(observed_species, file.path(OUT_DIR, "observed_species_summary.csv"))
  readr::write_csv(confusion_pool_summary, file.path(OUT_DIR, "confusion_pool_summary.csv"))

  # Audit TAXREF des espèces observées : appariement prioritaire par CD_NOM -> CD_REF.
  taxref_audit <- observed_species %>%
    mutate(
      in_taxref_by_cdref = !is.na(species_cd_ref),
      in_mainland_kept_status = species_key %in% taxref_mainland$species_key
    ) %>%
    left_join(
      taxref_mainland %>%
        transmute(
          species_key,
          taxref_mainland_label = species_label,
          taxref_direct_fr_status = fr,
          taxref_nom_valide = species_label
        ) %>%
        distinct(species_key, .keep_all = TRUE),
      by = "species_key"
    )

  readr::write_csv(taxref_audit, file.path(OUT_DIR, "taxref_cdref_audit_observed_species.csv"))

  message("Espèces observées au rang espèce : ", nrow(observed_species))
  message("Espèces rares <= ", RARE_MAX_TOTAL_ABUNDANCE, " individus : ", length(rare_species_keys))
  message("Genres observés avec au moins deux espèces observées : ",
          confusion_pool_summary %>% filter(source_pool == "observed", n_candidate_species >= 2) %>% nrow())

  save_step(STEP_ID, c(
    "taxref_mainland",
    "observed_species",
    "rare_species_keys",
    "observed_same_genus_pool",
    "taxref_same_genus_pool",
    "confusion_pool_summary",
    "taxref_audit"
  ))
}
