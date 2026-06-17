# ============================================================
# 03_prepare_environment_drivers — variables environnementales et drivers
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- "analyses/collembola taxo precision"
source(file.path(WORKFLOW_DIR,  "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

# Dépendances
load_required_step("01_import_clean")

STEP_ID <- "03_prepare_environment_drivers"

if (!should_rebuild_step(STEP_ID) && load_step(STEP_ID)) {
  message_header("03_prepare_environment_drivers — variables environnementales et drivers : cache chargé")
} else {
  # -----------------------------
  # 4. VARIABLES ENVIRONNEMENTALES ET DRIVERS
  # -----------------------------

  message_header("Preparation variables environnementales")

  # all_env_variables peut contenir deux colonnes altitude après clean_names : altitude / altitude_2.
  env_dat <- env_raw %>%
    mutate(station = as.character(station)) %>%
    {
      if ("altitude_2" %in% names(.)) {
        mutate(., altitude = coalesce(as.numeric(.data$altitude), as.numeric(.data$altitude_2)))
      } else {
        .
      }
    } %>%
    group_by(station) %>%
    summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      across(where(is.character), ~ first_non_missing(.x)),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

  if ("hab" %in% names(env_dat)) {
    env_dat <- env_dat %>% mutate(hab = as.factor(hab))
  } else {
    env_dat$hab <- factor(NA)
  }

  # habitat_broad issu du fichier collemboles si absent du fichier env.
  station_habitat_from_col <- col_dat %>%
    group_by(station) %>%
    summarise(habitat_broad = first_non_missing(habitat), .groups = "drop") %>%
    mutate(habitat_broad = as.factor(habitat_broad))

  env_dat <- env_dat %>%
    left_join(station_habitat_from_col, by = "station")

  station_report <- tibble(
    station = sort(unique(c(col_dat$station, env_dat$station))),
    in_collembola = station %in% col_dat$station,
    in_env = station %in% env_dat$station
  )
  readr::write_csv(station_report, file.path(OUT_DIR, "station_join_report.csv"))

  # Sélection a priori des drivers environnementaux : T360_mean, MOS et pH.
  # Ces variables sont conservées si elles existent, ont peu de NA et présentent assez de variance.
  available_numeric <- intersect(FOCAL_ENV_DRIVERS, names(env_dat))
  missing_focal_drivers <- setdiff(FOCAL_ENV_DRIVERS, names(env_dat))

  driver_screening_numeric <- purrr::map_dfr(FOCAL_ENV_DRIVERS, function(v) {
    if (!v %in% names(env_dat)) {
      return(tibble(
        variable = v,
        missing_prop = NA_real_,
        sd = NA_real_,
        n_unique = NA_integer_,
        keep_missing = FALSE,
        keep_variance = FALSE,
        keep_prelim = FALSE,
        reason = "missing_from_env_file"
      ))
    }
    x <- suppressWarnings(as.numeric(env_dat[[v]]))
    missing_prop <- mean(is.na(x))
    sd_x <- sd(x, na.rm = TRUE)
    n_unique <- n_distinct(x, na.rm = TRUE)
    keep_missing <- missing_prop <= MAX_MISSING_PROP_DRIVER
    keep_variance <- !is.na(sd_x) & sd_x > 0 & n_unique >= 3
    tibble(
      variable = v,
      missing_prop = missing_prop,
      sd = sd_x,
      n_unique = n_unique,
      keep_missing = keep_missing,
      keep_variance = keep_variance,
      keep_prelim = keep_missing & keep_variance,
      reason = case_when(
        !keep_missing ~ "too_many_missing_values",
        !keep_variance ~ "insufficient_variance",
        TRUE ~ "kept_a_priori"
      )
    )
  })

  selected_numeric_drivers_driverblock <- driver_screening_numeric %>%
    filter(keep_prelim) %>%
    pull(variable)

  selected_factor_drivers <- character(0)
  selected_drivers <- selected_numeric_drivers_driverblock

  # Audit de corrélation entre les trois drivers a priori. On ne les exclut pas automatiquement.
  if (length(selected_numeric_drivers_driverblock) >= 2) {
    driver_correlation_matrix <- env_dat %>%
      select(all_of(selected_numeric_drivers_driverblock)) %>%
      mutate(across(everything(), as.numeric)) %>%
      cor(use = "pairwise.complete.obs") %>%
      as.data.frame() %>%
      rownames_to_column("driver")
  } else {
    driver_correlation_matrix <- tibble()
  }

  # IMPORTANT : objet explicitement disponible en session après source().
  driver_selection_used <- tibble(
    variable = selected_drivers,
    type = rep("numeric_a_priori", length(selected_drivers)),
    ecological_role = recode(
      variable,
      t360_mean = "previous-year mean temperature",
      mos = "soil organic matter",
      ph = "soil pH",
      .default = "a priori driver"
    )
  )

  readr::write_csv(driver_screening_numeric, file.path(OUT_DIR, "driver_screening_numeric.csv"))
  readr::write_csv(driver_selection_used, file.path(OUT_DIR, "driver_selection_used.csv"))
  readr::write_csv(driver_correlation_matrix, file.path(OUT_DIR, "driver_correlation_matrix.csv"))

  message("Drivers focaux retenus : ", paste(selected_drivers, collapse = ", "))
  if (length(missing_focal_drivers) > 0) {
    warning("Drivers focaux absents du fichier environnemental : ", paste(missing_focal_drivers, collapse = ", "))
  }
  message("Objet R disponible : driver_selection_used")

  meta_station <- col_dat %>%
    group_by(station) %>%
    summarise(
      longitude = first_non_missing(x_longitude),
      latitude = first_non_missing(y_latitude),
      altitude_collembola = first_non_missing(altitude),
      habitat = first_non_missing(habitat),
      n_repetitions = n_distinct(repetition),
      .groups = "drop"
    ) %>%
    left_join(env_dat, by = "station")

  save_step(STEP_ID, c(
    "env_dat",
    "station_habitat_from_col",
    "station_report",
    "available_numeric",
    "missing_focal_drivers",
    "driver_screening_numeric",
    "selected_numeric_drivers_driverblock",
    "selected_factor_drivers",
    "selected_drivers",
    "driver_correlation_matrix",
    "driver_selection_used",
    "meta_station"
  ))
}
