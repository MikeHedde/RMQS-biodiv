# ============================================================
# taxon_profile_inputs.R
# Configuration template for the cross-taxon profile stage
# ============================================================
#
# Standard community file: site, species, abundance
# Standard taxonomy file: species, genus, family[, difficult_genus]
#
# Use absolute paths or paths relative to this file's parent directory.
# Keep only taxa for which both files exist; comment out the others.

TAXON_PROFILE_INPUTS <- list(
  Collembola = list(
    community_file = "inputs/collembola_community_long.csv",
    taxonomy_file = "inputs/collembola_regional_taxonomy.csv",
    community_delim = ",",
    taxonomy_delim = ","
  )

  #,
  # Earthworms = list(
  #   community_file = "inputs/earthworms_community_long.csv",
  #   taxonomy_file = "inputs/earthworms_regional_taxonomy.csv",
  #   community_delim = ",",
  #   taxonomy_delim = ","
  # ),
  # Ants = list(
  #   community_file = "inputs/ants_community_long.csv",
  #   taxonomy_file = "inputs/ants_regional_taxonomy.csv",
  #   community_delim = ",",
  #   taxonomy_delim = ","
  # )
)
