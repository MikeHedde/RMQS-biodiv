# ============================================================
# 01_import_clean — import, nettoyage et audit initial
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

STEP_ID <- "01_import_clean"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("01_import_clean — import, nettoyage et audit initial : cache chargé")
} else {
  # -----------------------------
  # 2. IMPORT ET NETTOYAGE
  # -----------------------------

  message_header("Import des donnees")

  col_raw <- safe_read_semicolon(PATH_COLLEMBOLA, encoding = "Latin1") %>%
    janitor::clean_names(ascii = TRUE)

  env_raw <- readr::read_csv(PATH_ENV, show_col_types = FALSE, name_repair = "unique") %>%
    janitor::clean_names(ascii = TRUE)

  taxref_raw <- safe_read_semicolon(PATH_TAXREF, encoding = "UTF-8") %>%
    janitor::clean_names(ascii = TRUE)

  # Alias robustes pour quelques variables dont le nettoyage des noms peut varier.

  required_col <- c("station", "repetition", "cd_nom", "lb_nom", "nom_valide", "rank", "famille", "ordre", "abondance_totale")
  missing_col <- setdiff(required_col, names(col_raw))
  if (length(missing_col) > 0) {
    stop("Colonnes manquantes dans collembola.csv : ", paste(missing_col, collapse = ", "))
  }

  # TAXREF est utilisé comme référentiel d'identifiants :
  # la base d'observation contient CD_NOM ; on récupère le CD_REF accepté via TAXREF.
  # Les noms restent seulement des labels lisibles.
  taxref_key_by_cdnom <- taxref_raw %>%
    mutate(
      cd_nom = clean_chr(cd_nom),
      cd_ref = clean_chr(cd_ref),
      rang = toupper(clean_chr(rang)),
      fr = clean_chr(fr),
      lb_nom = clean_chr(lb_nom),
      nom_valide = clean_chr(nom_valide),
      famille = clean_chr(famille),
      lb_binomial = extract_binomial(lb_nom),
      valid_binomial = extract_binomial(nom_valide)
    ) %>%
    transmute(
      cd_nom,
      taxref_cd_ref = cd_ref,
      taxref_rang = rang,
      taxref_fr = fr,
      taxref_lb_nom = lb_nom,
      taxref_nom_valide = nom_valide,
      taxref_lb_binomial = lb_binomial,
      taxref_valid_binomial = valid_binomial,
      taxref_famille = famille
    ) %>%
    distinct(cd_nom, .keep_all = TRUE)

  col_dat <- col_raw %>%
    mutate(
      station = as.character(station),
      repetition = as.character(repetition),
      sample_id = paste(station, repetition, sep = "_"),
      abundance = as.numeric(abondance_totale),
      cd_nom = clean_chr(cd_nom),
      lb_nom = clean_chr(lb_nom),
      nom_valide = clean_chr(nom_valide),
      rank = toupper(clean_chr(rank)),
      famille = clean_chr(famille),
      ordre = clean_chr(ordre),
      binomial_lb = extract_binomial(lb_nom),
      binomial_valid = extract_binomial(nom_valide)
    ) %>%
    left_join(taxref_key_by_cdnom, by = "cd_nom") %>%
    mutate(
      # Binôme lisible : priorité au nom valide TAXREF quand le CD_NOM matche.
      binomial = if_else(
        rank %in% c("S", "ES"),
        coalesce(taxref_valid_binomial, binomial_lb, binomial_valid),
        NA_character_
      ),
    
      # Clé espèce robuste : CD_REF TAXREF si disponible, sinon fallback explicite.
      species_cd_ref = if_else(rank %in% c("S", "ES") & !is.na(taxref_cd_ref), taxref_cd_ref, NA_character_),
      species_key = case_when(
        rank %in% c("S", "ES") & !is.na(species_cd_ref) ~ paste0("cdref:", species_cd_ref),
        rank %in% c("S", "ES") & !is.na(binomial) ~ paste0("name_fallback:", binomial),
        TRUE ~ NA_character_
      ),
      species_label = if_else(rank %in% c("S", "ES"), coalesce(taxref_valid_binomial, binomial), NA_character_),
    
      genus = case_when(
        !is.na(species_label) ~ stringr::word(species_label, 1),
        rank %in% c("G", "GN") & !is.na(lb_nom) ~ stringr::word(lb_nom, 1),
        TRUE ~ NA_character_
      ),
      taxo_level = case_when(
        rank %in% c("S", "ES") & !is.na(species_key) ~ "species",
        rank %in% c("G", "GN") & !is.na(genus) ~ "genus",
        rank %in% c("F", "FM") & !is.na(famille) ~ "family",
        rank %in% c("O", "OR") & !is.na(ordre) ~ "order",
        TRUE ~ "unusable"
      ),
    
      # RTU = reported taxonomic unit. Les espèces sont identifiées par CD_REF quand possible.
      rtu_best = case_when(
        taxo_level == "species" ~ paste0("species:", species_key),
        taxo_level == "genus" ~ paste0("genus:", genus),
        taxo_level == "family" ~ paste0("family:", famille),
        TRUE ~ NA_character_
      ),
      genus_unit = if_else(taxo_level %in% c("species", "genus") & !is.na(genus), paste0("genus:", genus), NA_character_),
      family_unit = if_else(taxo_level %in% c("species", "genus", "family") & !is.na(famille), paste0("family:", famille), NA_character_),
      species_unit = if_else(taxo_level == "species" & !is.na(species_key), paste0("species:", species_key), NA_character_)
    ) %>%
    filter(!is.na(abundance), abundance > 0)

  taxref_cdref_match_audit <- col_dat %>%
    mutate(
      taxref_match_status = case_when(
        taxo_level == "species" & !is.na(species_cd_ref) ~ "species_matched_to_cd_ref",
        taxo_level == "species" & is.na(species_cd_ref) ~ "species_unmatched_name_fallback",
        taxo_level %in% c("genus", "family", "order") & !is.na(taxref_cd_ref) ~ paste0(taxo_level, "_matched_to_cd_ref"),
        taxo_level %in% c("genus", "family", "order") & is.na(taxref_cd_ref) ~ paste0(taxo_level, "_unmatched"),
        TRUE ~ "unusable"
      )
    ) %>%
    count(taxo_level, rank, taxref_match_status, name = "n_rows") %>%
    arrange(taxo_level, rank, taxref_match_status)

  readr::write_csv(taxref_cdref_match_audit, file.path(OUT_DIR, "taxref_cdref_match_audit.csv"))

  taxref_unmatched_records <- col_dat %>%
    filter(is.na(taxref_cd_ref)) %>%
    distinct(cd_nom, lb_nom, nom_valide, rank, taxo_level, binomial, species_key, species_label) %>%
    arrange(rank, lb_nom)

  readr::write_csv(taxref_unmatched_records, file.path(OUT_DIR, "taxref_unmatched_records.csv"))

  # Audit : les individus seulement à l'ordre ne sont pas utilisés dans les matrices communautés.
  taxo_audit <- col_dat %>%
    group_by(taxo_level) %>%
    summarise(
      n_rows = n(),
      total_abundance = sum(abundance),
      n_lb_nom = n_distinct(lb_nom, na.rm = TRUE),
      n_species_binomial = n_distinct(binomial, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(taxo_level, c("species", "genus", "family", "order", "unusable")))

  readr::write_csv(taxo_audit, file.path(OUT_DIR, "taxonomic_resolution_audit.csv"))

  excluded_too_coarse <- col_dat %>%
    filter(taxo_level %in% c("order", "unusable"))
  readr::write_csv(excluded_too_coarse, file.path(OUT_DIR, "excluded_order_or_unusable_records.csv"))

  col_used <- col_dat %>%
    filter(taxo_level %in% c("species", "genus", "family"))

  message("Lignes totales collemboles : ", nrow(col_dat))
  message("Lignes utilisées dans les analyses communauté : ", nrow(col_used))
  message("Lignes exclues car ordre ou inutilisables : ", nrow(excluded_too_coarse))
  message("Abondance exclue car ordre ou inutilisable : ", sum(excluded_too_coarse$abundance, na.rm = TRUE))
  message("Stations : ", n_distinct(col_dat$station))
  message("Echantillons station x repetition : ", n_distinct(col_dat$sample_id))

  save_step(STEP_ID, c(
    "col_raw",
    "env_raw",
    "taxref_raw",
    "taxref_key_by_cdnom",
    "col_dat",
    "taxref_cdref_match_audit",
    "taxref_unmatched_records",
    "taxo_audit",
    "excluded_too_coarse",
    "col_used"
  ))
}
