# Orchestration du workflow

# 00_run_all.R
source("analyses/Occupancy/01_config.R")
source("analyses/Occupancy/02_packages.R")
source("analyses/Occupancy/03_read_parse.R")   # crée dat0_all, hab_present, out_dir

# Conserver le out_dir racine
out_dir_root <- out_dir

# Si tu veux garder aussi une analyse "tous habitats", laisse ce bloc avant la boucle
source("analyses/Occupancy/04_detection_effort.R")
source("analyses/Occupancy/05_occupancy_fit.R")
source("analyses/Occupancy/06_occupancy_figs.R")
source("analyses/Occupancy/07_traits_detectability.R")
#source("analyses/Occupancy/08_diversity_indices.R")
#source("analyses/Occupancy/09_composition_beta.R")
#source("analyses/Occupancy//11_accumulation_curves.R")
#message("=== FIN WORKFLOW ===  Résultats : ", normalizePath(out_dir))

for (hab in hab_present) {
  message("=== RUN HABITAT_OC: ", hab, " ===")
  # 1) Filtrer les données pour cet habitat
  dat0 <- dat0_all %>% dplyr::filter(HABITAT_OC == hab)
  
  # 2) Dossier de sortie spécifique
  out_dir <- file.path(out_dir_root, "by_HABITAT_OC", hab)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 3) Lancer le pipeline inchangé (il écrira dans out_dir)
  source("analyses/Occupancy/04_detection_effort.R")
  source("analyses/Occupancy/05_occupancy_fit.R")
  source("analyses/Occupancy/06_occupancy_figs.R")
  source("analyses/Occupancy/07_traits_detectability.R")
  source("R/08_diversity_indices.R")
  source("R/09_composition_beta.R")
  source("R/11_accumulation_curves.R")
}

message("=== FIN par habitat — résultats dans ", normalizePath(file.path(out_dir_root, "by_HABITAT_OC")), " ===")
