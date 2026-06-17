# ============================================================
# 08_rare_and_inventory_audits — espèces rares et inventaires par scénario
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances
load_required_step("01_import_clean")
load_required_step("02_prepare_taxonomy_pools")
load_required_step("04_build_scenarios")

STEP_ID <- "08_rare_and_inventory_audits"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("08_rare_and_inventory_audits — espèces rares et inventaires par scénario : cache chargé")
} else {
  # -----------------------------
  # 9. ESPECES RARES : AUDIT PROPRE
  # -----------------------------

  message_header("Audit des espèces rares")

  rare_contribution <- col_used %>%
    mutate(
      rare_status = case_when(
        taxo_level == "species" & species_key %in% rare_species_keys ~ "rare_species",
        taxo_level == "species" ~ "non_rare_species",
        TRUE ~ "non_species_RTU"
      )
    ) %>%
    group_by(rare_status) %>%
    summarise(
      abundance = sum(abundance),
      n_rows = n(),
      n_lb_nom = n_distinct(lb_nom, na.rm = TRUE),
      n_species_binomial = n_distinct(binomial, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(rare_contribution, file.path(OUT_DIR, "rare_species_contribution.csv"))

  # -----------------------------


  # -----------------------------
  # 5b. INVENTAIRES TAXONOMIQUES PAR SCENARIO
  # -----------------------------

  scenario_inventory_by_iter <- scenario_data %>%
    distinct(
      scenario_family,
      scenario,
      baseline_scenario,
      unit_type,
      iter,
      taxon_unit
    ) %>%
    arrange(scenario_family, scenario, iter, taxon_unit)

  readr::write_csv(
    scenario_inventory_by_iter,
    file.path(OUT_DIR, "scenario_inventory_by_iter.csv")
  )

  # Changements d'inventaire par rapport à la baseline propre à chaque famille.
  # Pour chaque scénario, on mesure :
  # - unités retenues depuis la baseline
  # - unités perdues par rapport à la baseline
  # - unités ajoutées par rapport à la baseline

  baseline_inventory <- scenario_inventory_by_iter %>%
    filter(scenario == baseline_scenario) %>%
    group_by(scenario_family, baseline_scenario) %>%
    summarise(
      baseline_units = list(unique(taxon_unit)),
      gamma_baseline = n_distinct(taxon_unit),
      .groups = "drop"
    )

  scenario_inventory_change <- scenario_inventory_by_iter %>%
    filter(scenario != baseline_scenario) %>%
    group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
    summarise(
      scenario_units = list(unique(taxon_unit)),
      gamma_scenario = n_distinct(taxon_unit),
      .groups = "drop"
    ) %>%
    left_join(
      baseline_inventory,
      by = c("scenario_family", "baseline_scenario")
    ) %>%
    mutate(
      n_retained = purrr::map2_int(
        scenario_units,
        baseline_units,
        ~ length(intersect(.x, .y))
      ),
      n_lost = purrr::map2_int(
        baseline_units,
        scenario_units,
        ~ length(setdiff(.x, .y))
      ),
      n_gained = purrr::map2_int(
        scenario_units,
        baseline_units,
        ~ length(setdiff(.x, .y))
      ),
      gamma_change = gamma_scenario - gamma_baseline,
      gamma_change_pct = 100 * gamma_change / gamma_baseline,
      lost_pct = 100 * n_lost / gamma_baseline,
      gained_pct = 100 * n_gained / gamma_baseline
    ) %>%
    select(
      scenario_family,
      scenario,
      baseline_scenario,
      unit_type,
      iter,
      gamma_baseline,
      gamma_scenario,
      gamma_change,
      gamma_change_pct,
      n_retained,
      n_lost,
      n_gained,
      lost_pct,
      gained_pct
    )

  readr::write_csv(
    scenario_inventory_change,
    file.path(OUT_DIR, "scenario_inventory_change.csv")
  )

  save_step(STEP_ID, c(
    "rare_contribution",
    "scenario_inventory_by_iter",
    "baseline_inventory",
    "scenario_inventory_change"
  ))
}
