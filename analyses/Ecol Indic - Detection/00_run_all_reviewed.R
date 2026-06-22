# Orchestration du workflow — reviewer-ready version

source("analyses/Ecol Indic - Detection/01_packages.R")
source("analyses/Ecol Indic - Detection/02_config.R")
source("analyses/Ecol Indic - Detection/03_read_parse.R")

# Detection / occupancy
source("analyses/Ecol Indic - Detection/04_detection_effort_reviewed.R")
source("analyses/Ecol Indic - Detection/05_occupancy_fit_reviewed.R")
source("analyses/Ecol Indic - Detection/15_gpd_IL_sensitivity_refit.R")

# Occupancy outputs and stats
source("analyses/Ecol Indic - Detection/06_occupancy_figs.R")
source("analyses/Ecol Indic - Detection/07_occupancy_stat.R")
source("analyses/Ecol Indic - Detection/09_traits_detectability_reviewed.R")

# Diversity
source("analyses/Ecol Indic - Detection/10_diversity_indices.R")
source("analyses/Ecol Indic - Detection/11_diversity_stat.R")

# Reviewer-requested tables and sensitivity analyses
source("analyses/Ecol Indic - Detection/12_supplementary_reviewer_response_analyses.R")

message("=== FIN WORKFLOW REVIEWED === Résultats : ", normalizePath(out_dir))
