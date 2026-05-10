# Orchestration du workflow — reviewer-ready version

source("analyses/Detection_spider/01_packages.R")
source("analyses/Detection_spider/02_config.R")
source("analyses/Detection_spider/03_read_parse.R")

# Detection / occupancy
source("analyses/Detection_spider/04_detection_effort_reviewed.R")
source("analyses/Detection_spider/05_occupancy_fit_reviewed.R")

# Occupancy outputs and stats
source("analyses/Detection_spider/06_occupancy_figs.R")
source("analyses/Detection_spider/07_occupancy_stat.R")
source("analyses/Detection_spider/08_DARMHa checks.R")
source("analyses/Detection_spider/09_traits_detectability.R")

# Diversity
source("analyses/Detection_spider/10_diversity_indices.R")
source("analyses/Detection_spider/11_diversity_stat.R")
source("analyses/Detection_spider/12_composition_beta.R")
source("analyses/Detection_spider/13_accumulation_curves.R")

# Reviewer-requested tables and sensitivity analyses
source("analyses/Detection_spider/14_reviewer_analyses.R")

message("=== FIN WORKFLOW REVIEWED === Résultats : ", normalizePath(out_dir))
