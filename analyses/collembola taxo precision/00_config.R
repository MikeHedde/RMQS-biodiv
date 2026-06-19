# ============================================================
# 00_config.R
# Configuration commune du workflow taxonomic uncertainty
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(vegan)
  library(divent)
  library(broom)
  library(patchwork)
  library(forcats)
})

# Dossier projet : par défaut, le working directory R.
# Si les scripts sont dans un sous-dossier, tu peux faire avant source() :
# options(taxo.project_dir = "C:/chemin/vers/le/projet")
PROJECT_DIR <- getOption("taxo.project_dir", getwd())

resolve_project_path <- function(path) {
  if (grepl("^([A-Za-z]:)?[\\/]", path)) {
    path
  } else {
    file.path(PROJECT_DIR, path)
  }
}

PATH_COLLEMBOLA <- resolve_project_path("data/raw-data/1.faune/collembola.csv")
PATH_ENV        <- resolve_project_path("data/derived-data/all_env_variables.csv")
PATH_TAXREF     <- resolve_project_path("data/raw-data/taxref_collembola.csv")

OUT_DIR <- resolve_project_path("outputs_taxonomic_uncertainty_v9")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CACHE_DIR <- file.path(OUT_DIR, "_rds_cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

PUB_DIR <- file.path(OUT_DIR, "figures_publication_ready")
dir.create(PUB_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

ANALYSIS_UNIT <- "station"
N_SIM <- 10

RUN_DRIVER_BLOCK <- TRUE
PERMANOVA_N_PERM <- 999
PERMANOVA_MAX_STOCHASTIC_ITER <- 100
PERMANOVA_MAX_NUMERIC_DRIVERS <- 3

GENUS_CONFUSION_RATES <- c(0.03, 0.05, 0.10)
RARE_WEIGHTED_MEAN_ERROR <- c(0.03, 0.05, 0.10)
RARE_WEIGHTED_MAX_ERROR <- 0.25
RARE_MAX_TOTAL_ABUNDANCE <- 3
MAINLAND_FR_STATUS_KEEP <- c("P", "E", "C")

FOCAL_ENV_DRIVERS <- c("t360_mean", "mos", "p_h")
DRIVER_CANDIDATES_NUMERIC <- FOCAL_ENV_DRIVERS
DRIVER_CANDIDATES_FACTOR <- character(0)
MAX_MISSING_PROP_DRIVER <- 0.25
MAX_ABS_COR_DRIVER <- 0.95

DIVENT_HILL_ESTIMATOR <- "naive"
DIVENT_COVERAGE_ESTIMATORS <- c("ZhangHuang", "Chao", "Turing")
COVERAGE_STANDARDIZATION_LEVEL <- NULL


ID_ERROR_APPENDIX_RATES <- seq(0.01, 0.20, by = 0.01)
ID_ERROR_APPENDIX_N_SIM <- getOption("taxo.appendix_n_sim", N_SIM)
ID_ERROR_APPENDIX_POOL <- getOption("taxo.appendix_pool", "observed")  # or "taxref_mainland"
ID_ERROR_APPENDIX_CRITICAL_REL_CHANGE <- getOption("taxo.appendix_critical_rel_change", 0.05)
ID_ERROR_APPENDIX_CRITICAL_STABILITY <- getOption("taxo.appendix_critical_stability", 0.95)

cache_path <- function(step_id) {
  file.path(CACHE_DIR, paste0(step_id, ".rds"))
}

should_rebuild_step <- function(step_id) {
  isTRUE(getOption("taxo.force_rebuild", FALSE)) ||
    step_id %in% getOption("taxo.force_steps", character(0))
}

save_step <- function(step_id, objects, envir = parent.frame()) {
  found <- objects[objects %in% ls(envir = envir)]
  missing <- setdiff(objects, found)
  if (length(missing) > 0) {
    warning("Objets non trouvés et non sauvegardés pour ", step_id, " : ", paste(missing, collapse = ", "))
  }
  saveRDS(mget(found, envir = envir), cache_path(step_id))
  message("Cache RDS écrit : ", cache_path(step_id))
  invisible(cache_path(step_id))
}

load_step <- function(step_id, envir = parent.frame()) {
  path <- cache_path(step_id)
  if (!file.exists(path)) return(FALSE)
  objs <- readRDS(path)
  list2env(objs, envir = envir)
  message("Cache RDS chargé : ", path)
  TRUE
}

load_required_step <- function(step_id, envir = parent.frame()) {
  if (!load_step(step_id, envir = envir)) {
    stop("Cache manquant pour l'étape ", step_id, ". Lance d'abord le script correspondant.")
  }
  invisible(TRUE)
}

options(
  taxo.gbif_mode = "download", # ou "auto" / "search"
  taxo.gbif_buffer_km = 250,
  taxo.gbif_n_sim = 20,
  taxo.gbif_pool_strategy = "buffer_then_region_fallback",
  taxo.gbif_weighted_shuffling = TRUE
)