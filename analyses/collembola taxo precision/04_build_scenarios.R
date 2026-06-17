# ============================================================
# 04_build_scenarios — scénarios taxonomiques déterministes et stochastiques
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances
load_required_step("01_import_clean")
load_required_step("02_prepare_taxonomy_pools")

STEP_ID <- "04_build_scenarios"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("04_build_scenarios — scénarios taxonomiques déterministes et stochastiques : cache chargé")
} else {
  # -----------------------------
  # 5. SCENARIOS PROPRES D'INCERTITUDE
  # -----------------------------

  message_header("Construction des scénarios propres")

  base_cols <- c("station", "repetition", "sample_id", "abundance", "taxo_level", "binomial", "species_key", "species_label", "genus", "famille", "ordre")

  mk_scenario <- function(dat, scenario_name, scenario_family, baseline_scenario, taxon_col, iter = 1, unit_type = "RTU") {
    dat %>%
      mutate(
        taxon_unit = .data[[taxon_col]],
        scenario = scenario_name,
        scenario_family = scenario_family,
        baseline_scenario = baseline_scenario,
        iter = iter,
        unit_type = unit_type
      ) %>%
      filter(!is.na(taxon_unit)) %>%
      select(all_of(base_cols), taxon_unit, scenario, scenario_family, baseline_scenario, iter, unit_type)
  }

  # Baselines et scénarios déterministes.
  scenario_det <- bind_rows(
    # Workflow sur RTU mixtes : espèces + genres + familles, ordre exclu.
    mk_scenario(
      col_used,
      "rtu_mixed_best_available",
      "workflow_filtering",
      "rtu_mixed_best_available",
      "rtu_best",
      unit_type = "mixed_RTU"
    ),
    # Sous-ensemble strictement spécifique : ce n'est PAS une référence RTU mixte.
    mk_scenario(
      col_used %>% filter(taxo_level == "species"),
      "rtu_species_only_drop_unresolved",
      "workflow_filtering",
      "rtu_mixed_best_available",
      "species_unit",
      unit_type = "species_only"
    ),
    # Espèces rares rabattues au genre, en conservant genres/familles déjà non spécifiques.
    col_used %>%
      mutate(
        rare_to_genus_unit = case_when(
          taxo_level == "species" & species_key %in% rare_species_keys ~ genus_unit,
          TRUE ~ rtu_best
        )
      ) %>%
      mk_scenario(
        "rtu_rare_species_to_genus",
        "workflow_filtering",
        "rtu_mixed_best_available",
        "rare_to_genus_unit",
        unit_type = "mixed_RTU"
      ),
  
    # Passage au genre sur le sous-ensemble réellement résoluble au genre.
    mk_scenario(
      col_used %>% filter(taxo_level %in% c("species", "genus")),
      "rtu_mixed_genus_resolvable",
      "resolution_genus",
      "rtu_mixed_genus_resolvable",
      "rtu_best",
      unit_type = "mixed_RTU_genus_resolvable"
    ),
    mk_scenario(
      col_used %>% filter(taxo_level %in% c("species", "genus")),
      "rtu_genus_level",
      "resolution_genus",
      "rtu_mixed_genus_resolvable",
      "genus_unit",
      unit_type = "genus_units"
    ),
  
    # Passage à la famille sur espèces + genres + familles. Ordre exclu.
    mk_scenario(
      col_used,
      "rtu_mixed_family_resolvable",
      "resolution_family",
      "rtu_mixed_family_resolvable",
      "rtu_best",
      unit_type = "mixed_RTU_family_resolvable"
    ),
    mk_scenario(
      col_used,
      "rtu_family_level",
      "resolution_family",
      "rtu_mixed_family_resolvable",
      "family_unit",
      unit_type = "family_units"
    ),
  
    # Baseline strictement spécifique pour les scénarios d'erreur entre espèces.
    mk_scenario(
      col_used %>% filter(taxo_level == "species"),
      "species_baseline",
      "species_error_strict",
      "species_baseline",
      "species_unit",
      unit_type = "species_units"
    )
  )

  # Probabilités d'erreur pondérées par la rareté.
  species_abundance <- col_used %>%
    filter(taxo_level == "species") %>%
    group_by(species_key) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop")

  make_rare_weight_error <- function(mean_error) {
    tmp <- species_abundance %>%
      mutate(raw_weight = 1 / sqrt(total_abundance))
    scale_factor <- mean_error / weighted.mean(tmp$raw_weight, w = tmp$total_abundance)
    tmp %>%
      transmute(
        species_key,
        p_error = pmin(raw_weight * scale_factor, RARE_WEIGHTED_MAX_ERROR)
      )
  }

  choose_same_genus_candidate <- function(genus, current_species_key, pool_tbl) {
    cand <- pool_tbl %>%
      filter(.data$genus == !!genus, .data$candidate_key != !!current_species_key) %>%
      distinct(candidate_key, candidate_label, genus)
  
    if (nrow(cand) == 0) {
      return(tibble(candidate_key = NA_character_, candidate_label = NA_character_, genus = NA_character_))
    }
  
    cand %>% slice_sample(n = 1)
  }

  simulate_same_genus_error <- function(dat_species, pool_tbl, scenario_name, iter_id, fixed_rate = NULL, p_tbl = NULL) {
    dat0 <- dat_species %>%
      filter(taxo_level == "species", !is.na(species_key), !is.na(genus)) %>%
      select(all_of(base_cols)) %>%
      mutate(
        p_error = if (!is.null(fixed_rate)) fixed_rate else NA_real_
      )
  
    if (!is.null(p_tbl)) {
      dat0 <- dat0 %>%
        select(-p_error) %>%
        left_join(p_tbl, by = "species_key") %>%
        mutate(p_error = replace_na(p_error, 0))
    }
  
    purrr::pmap_dfr(dat0, function(station, repetition, sample_id, abundance, taxo_level,
                                   binomial, species_key, species_label, genus, famille, ordre, p_error) {
      n <- as.integer(round(abundance))
      if (is.na(p_error) || p_error <= 0 || n <= 0) {
        return(tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = abundance, taxo_level = taxo_level, binomial = binomial,
          species_key = species_key, species_label = species_label,
          genus = genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", species_key)
        ))
      }
    
      target <- choose_same_genus_candidate(genus, species_key, pool_tbl)
      if (nrow(target) == 0 || is.na(target$candidate_key[1])) {
        return(tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = abundance, taxo_level = taxo_level, binomial = binomial,
          species_key = species_key, species_label = species_label,
          genus = genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", species_key)
        ))
      }
    
      target_key <- target$candidate_key[1]
      target_label <- target$candidate_label[1]
      target_genus <- target$genus[1]
    
      n_swap <- rbinom(1, size = n, prob = p_error)
      if (n_swap == 0) {
        return(tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = abundance, taxo_level = taxo_level, binomial = binomial,
          species_key = species_key, species_label = species_label,
          genus = genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", species_key)
        ))
      }
    
      bind_rows(
        if (n - n_swap > 0) {
          tibble(
            station = station, repetition = repetition, sample_id = sample_id,
            abundance = n - n_swap, taxo_level = taxo_level, binomial = binomial,
            species_key = species_key, species_label = species_label,
            genus = genus, famille = famille, ordre = ordre,
            taxon_unit = paste0("species:", species_key)
          )
        },
        tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = n_swap, taxo_level = taxo_level, binomial = target_label,
          species_key = target_key, species_label = target_label,
          genus = target_genus, famille = famille, ordre = ordre,
          taxon_unit = paste0("species:", target_key)
        )
      )
    }) %>%
      mutate(
        scenario = scenario_name,
        scenario_family = "species_error_strict",
        baseline_scenario = "species_baseline",
        iter = iter_id,
        unit_type = "species_units"
      )
  }

  species_records <- col_used %>% filter(taxo_level == "species")

  scenario_stoch <- purrr::map_dfr(seq_len(N_SIM), function(iter_i) {
    set.seed(1000 + iter_i)
  
    fixed_observed <- purrr::map_dfr(GENUS_CONFUSION_RATES, function(rate_i) {
      simulate_same_genus_error(
        species_records,
        pool_tbl = observed_same_genus_pool,
        scenario_name = paste0("species_same_genus_observed_", round(rate_i * 100), "pct"),
        iter_id = iter_i,
        fixed_rate = rate_i
      )
    })
  
    rare_weighted <- purrr::map_dfr(RARE_WEIGHTED_MEAN_ERROR, function(mean_i) {
      simulate_same_genus_error(
        species_records,
        pool_tbl = observed_same_genus_pool,
        scenario_name = paste0("species_same_genus_observed_rare_weighted_mean_", round(mean_i * 100), "pct"),
        iter_id = iter_i,
        p_tbl = make_rare_weight_error(mean_i)
      )
    })
  
    taxref_fixed <- purrr::map_dfr(GENUS_CONFUSION_RATES, function(rate_i) {
      simulate_same_genus_error(
        species_records,
        pool_tbl = taxref_same_genus_pool,
        scenario_name = paste0("species_same_genus_taxref_mainland_", round(rate_i * 100), "pct"),
        iter_id = iter_i,
        fixed_rate = rate_i
      )
    })
  
    bind_rows(fixed_observed, rare_weighted, taxref_fixed)
  })

  scenario_data <- bind_rows(scenario_det, scenario_stoch) %>%
    mutate(
      scenario = as.character(scenario),
      scenario_family = as.character(scenario_family),
      baseline_scenario = as.character(baseline_scenario)
    )

  scenario_definitions <- scenario_data %>%
    distinct(scenario, scenario_family, baseline_scenario, unit_type) %>%
    arrange(scenario_family, scenario)
  readr::write_csv(scenario_definitions, file.path(OUT_DIR, "scenario_definitions.csv"))

  readr::write_csv(scenario_definitions, file.path(OUT_DIR, "scenario_definitions.csv"))
  message("Scénarios construits : ", n_distinct(scenario_data$scenario))
  print(scenario_definitions)


  # Ordre d'affichage par logique analytique, sauvegardé pour figures.R.
  scenario_order <- scenario_definitions %>%
    mutate(order_key = case_when(
      scenario == baseline_scenario ~ 0,
      scenario_family == "workflow_filtering" ~ 1,
      scenario_family == "resolution_genus" ~ 2,
      scenario_family == "resolution_family" ~ 3,
      scenario_family == "species_error_strict" ~ 4,
      TRUE ~ 5
    )) %>%
    arrange(order_key, scenario_family, scenario) %>%
    pull(scenario) %>%
    unique()

  save_step(STEP_ID, c(
    "scenario_det",
    "species_abundance",
    "species_records",
    "scenario_stoch",
    "scenario_data",
    "scenario_definitions",
    "scenario_order"
  ))
}
