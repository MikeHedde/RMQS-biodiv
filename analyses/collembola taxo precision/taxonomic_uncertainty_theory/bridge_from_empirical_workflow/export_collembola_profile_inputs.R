# ============================================================
# export_collembola_profile_inputs.R
# One-way bridge: empirical Collembola workflow -> theory module
# ============================================================
#
# Run this ONLY after the empirical workflow has created the objects
# `col_used` and `taxref_mainland` (steps 01 and 02).
# This script does not alter empirical analyses or scenario outputs.
# It only exports two standardised files used by the theoretical module.

if (!exists("col_used") || !exists("taxref_mainland")) {
  stop(
    "Objects `col_used` and/or `taxref_mainland` are missing. Source 01_import_clean.R and 02_prepare_taxonomy_pools.R first.",
    call. = FALSE
  )
}

PROFILE_EXPORT_DIR <- getOption(
  "taxo.profile_export_dir",
  file.path(getwd(), "theory_profile_inputs")
)
dir.create(PROFILE_EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)

# Strict species-level community: this is the latent taxonomic matrix used to
# describe the empirical group, not a replacement for mixed-RTU analyses.
collembola_community_long <- col_used %>%
  filter(taxo_level == "species", !is.na(species_key)) %>%
  transmute(
    site = as.character(station),
    species = as.character(species_key),
    abundance = as.numeric(abundance)
  ) %>%
  group_by(site, species) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# Regional French species pool. Difficult genera can be edited below or joined
# from an expert table with two columns: genus, difficult_genus.
DIFFICULT_GENERA <- c(
  "Folsomia", "Mesaphorura", "Isotoma", "Lepidocyrtus",
  "Entomobrya", "Ceratophysella", "Sminthurides"
)

collembola_regional_taxonomy <- taxref_mainland %>%
  transmute(
    species = as.character(species_key),
    genus = as.character(genus),
    family = as.character(famille),
    difficult_genus = genus %in% DIFFICULT_GENERA
  ) %>%
  filter(!is.na(species), !is.na(genus), !is.na(family)) %>%
  distinct(species, .keep_all = TRUE)

readr::write_csv(
  collembola_community_long,
  file.path(PROFILE_EXPORT_DIR, "collembola_community_long.csv")
)
readr::write_csv(
  collembola_regional_taxonomy,
  file.path(PROFILE_EXPORT_DIR, "collembola_regional_taxonomy.csv")
)

message("Profile inputs written to: ", normalizePath(PROFILE_EXPORT_DIR))
