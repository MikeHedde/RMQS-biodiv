# ============================================================
# PARAMÈTRES GLOBAUX
# ============================================================
in_path        <- "data/raw-data/1.faune/carabidae.csv"
out_dir        <- "figures/detectability/carabidae"
trap_coord_path<- "data/derived-data/pitfall_trap_coordinates.csv"
target_project <- "RMQS_2025"
taxo_level     <- c("ES", "S")
min_sites      <- 7
require_contrast <- TRUE

# Méthodes analysées
methods_use <- c("Pitfall", "GPD")   # parmi "Pitfall","GPD","DVAC","TRI"
require_complete_pitfall_subsets <- FALSE

# Réplicats Pitfall
pit_reps <- c(10, 8, 6, 4, 2)

# Constantes effort
pit_diam_m         <- 0.05
pit_circ_m         <- pi * pit_diam_m
gpd_fence_length_m <- 1.0
gpd_open_circ_m    <- pi * pit_diam_m
gpd_cross_m        <- gpd_fence_length_m
gpd_unit_intensity <- gpd_open_circ_m + gpd_cross_m
tm_unit            <- 0.0625

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
