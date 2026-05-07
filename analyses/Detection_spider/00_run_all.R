# Orchestration du workflow

# preparation
source("analyses/Ecol Indic - Detection/01_packages.R")
source("analyses/Ecol Indic - Detection/02_config.R")
source("analyses/Ecol Indic - Detection/03_read_parse.R")   

# calcul detection
source("analyses/Ecol Indic - Detection/04_detection_effort.R")
source("analyses/Ecol Indic - Detection/05_occupancy_fit.R")

# stat detection
source("analyses/Ecol Indic - Detection/06_occupancy_figs.R")
source("analyses/Ecol Indic - Detection/07_occupancy_stat.R")
source("analyses/Ecol Indic - Detection/08_DARMHa checks.R")
source("analyses/Ecol Indic - Detection/09_traits_detectability.R")

# diversité
source("analyses/Ecol Indic - Detection/10_diversity_indices.R")
source("analyses/Ecol Indic - Detection/11_diversity_stat.R")
source("analyses/Occupancy/12_composition_beta.R")

source("analyses/Ecol Indic - Detection//13_accumulation_curves.R")
#message("=== FIN WORKFLOW ===  Résultats : ", normalizePath(out_dir))
