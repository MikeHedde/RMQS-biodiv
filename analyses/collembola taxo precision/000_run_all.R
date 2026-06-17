# ============================================================
# run_all.R
# Lance tout le workflow dans l'ordre.
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- if (!is.na(.workflow_file)) dirname(.workflow_file) else getwd()

# Pour forcer le recalcul complet :
# options(taxo.force_rebuild = TRUE)
# Pour forcer seulement certaines étapes :
# options(taxo.force_steps = c("04_build_scenarios", "05_compute_alpha_metrics"))

source(file.path(WORKFLOW_DIR, "01_import_clean.R"))
source(file.path(WORKFLOW_DIR, "02_prepare_taxonomy_pools.R"))
source(file.path(WORKFLOW_DIR, "03_prepare_environment_drivers.R"))
source(file.path(WORKFLOW_DIR, "04_build_scenarios.R"))
source(file.path(WORKFLOW_DIR, "05_compute_alpha_metrics.R"))
source(file.path(WORKFLOW_DIR, "06_compute_stability.R"))
source(file.path(WORKFLOW_DIR, "07_compute_drivers.R"))
source(file.path(WORKFLOW_DIR, "08_rare_and_inventory_audits.R"))
source(file.path(WORKFLOW_DIR, "09_appendix_error_gradient_analysis.R"))
source(file.path(WORKFLOW_DIR, "10_gbif_regional_error_gradient.R"))
source(file.path(WORKFLOW_DIR, "figures.R"))
source(file.path(WORKFLOW_DIR, "appendix_figures_error_gradient.R"))


message("Workflow terminé.")
