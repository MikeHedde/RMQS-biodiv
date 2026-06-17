# ============================================================
# 10_gbif_regional_error_gradient.R
# ------------------------------------------------------------
# Objectif : construire un pool de confusion taxonomique informe par GBIF
#            et relancer un gradient d'erreur d'identification de 1% a 20%.
#
# Principe : une espece candidate au shuffling doit etre :
#   1) congenerique ;
#   2) connue de GBIF dans une zone regionale large ;
#   3) idealement observee dans un buffer autour de la station ;
#   4) ou observee dans ton propre jeu de donnees a cette station.
#
# Dependances du workflow :
#   - R/00_config.R
#   - R/functions_taxonomic_uncertainty.R
#   - cache 01_import_clean.rds
#   - cache 02_prepare_taxonomy_pools.rds
#
# Packages additionnels :
#   install.packages(c("rgbif", "sf"))
#   install.packages("CoordinateCleaner") # optionnel mais recommande
#
# GBIF : pour l'analyse finale, creer un compte GBIF puis definir dans ~/.Renviron :
#   GBIF_USER=ton_username_gbif
#   GBIF_PWD=ton_password_gbif
#   GBIF_EMAIL=ton.email@domaine.fr
#
# Puis redemarrer R.
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- if (!is.na(.workflow_file)) dirname(.workflow_file) else getwd()
source(file.path(WORKFLOW_DIR, "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(rgbif)
  library(sf)
  library(forcats)
})

# ------------------------------------------------------------
# 0. PARAMETRES MODIFIABLES
# ------------------------------------------------------------

# Pays retenus pour definir le contexte regional large autour de la France.
# A adapter si ton jeu de donnees couvre seulement une partie de la France.
GBIF_COUNTRIES <- getOption(
  "taxo.gbif_countries",
  c("FR", "ES", "IT", "CH", "DE", "BE", "LU", "AD", "MC")
)

# Projection metrique europeenne pour buffers et AOO.
GBIF_CRS_METRIC <- getOption("taxo.gbif_crs_metric", 3035)

# Buffer autour de chaque station pour le filtre regional local.
GBIF_BUFFER_KM <- getOption("taxo.gbif_buffer_km", 250)

# Grille AOO 2 x 2 km, standard classiquement utilise pour l'AOO IUCN.
GBIF_AOO_CELLSIZE_M <- getOption("taxo.gbif_aoo_cellsize_m", 2000)

# Filtre de precision des coordonnees.
GBIF_MAX_COORD_UNCERTAINTY_M <- getOption("taxo.gbif_max_coord_uncertainty_m", 10000)

# Gradient d'erreur teste.
GBIF_ERROR_RATES <- getOption("taxo.gbif_error_rates", seq(0.01, 0.20, by = 0.01))
GBIF_N_SIM <- getOption("taxo.gbif_n_sim", 20)

# Mode GBIF :
#   "auto"     : occ_download si identifiants GBIF disponibles, sinon occ_search.
#   "download" : occ_download, recommande pour l'analyse finale et pour citation GBIF.
#   "search"   : occ_search, pratique pour un test rapide mais moins propre pour publication.
GBIF_MODE <- getOption("taxo.gbif_mode", "auto")

# Mode search seulement : limite par espece. Ne pas utiliser pour une analyse finale.
GBIF_SEARCH_LIMIT_PER_SPECIES <- getOption("taxo.gbif_search_limit_per_species", 500)

# Strategie du pool station-specifique :
#   "buffer_or_observed"        : candidats observes dans le buffer GBIF ou dans ton dataset a la station.
#   "buffer_then_region_fallback": ajoute le pool regional large si station x genre a < 2 candidats.
#   "country_region"           : tous les candidats connus dans les pays GBIF_COUNTRIES pour toutes les stations.
GBIF_POOL_STRATEGY <- getOption("taxo.gbif_pool_strategy", "buffer_then_region_fallback")

# Tirage pondere des candidats.
GBIF_WEIGHTED_SHUFFLING <- getOption("taxo.gbif_weighted_shuffling", TRUE)

# Seuils d'interpretation graphique.
GBIF_CRITICAL_REL_CHANGE <- getOption("taxo.gbif_critical_rel_change", 0.05)
GBIF_CRITICAL_STABILITY <- getOption("taxo.gbif_critical_stability", 0.95)

# Si TRUE, recharge les fichiers GBIF/intermediaires existants si presents.
GBIF_USE_CACHE <- getOption("taxo.gbif_use_cache", TRUE)

# ------------------------------------------------------------
# 1. DEPENDANCES DU WORKFLOW
# ------------------------------------------------------------

load_required_step("01_import_clean")
load_required_step("02_prepare_taxonomy_pools")

message_header("GBIF regional shuffling — initialisation")

GBIF_DIR <- file.path(OUT_DIR, "gbif_regional_shuffling")
GBIF_FIG_DIR <- file.path(PUB_DIR, "gbif_regional_shuffling")
dir.create(GBIF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GBIF_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

safe_scalar_chr <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  as.character(x[1])
}

safe_scalar_int <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) return(NA_integer_)
  suppressWarnings(as.integer(x[1]))
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

# ------------------------------------------------------------
# 2. POOL TAXONOMIQUE CANDIDAT
# ------------------------------------------------------------

message_header("Construction du pool taxonomique candidat")

candidate_taxonomic_pool <- bind_rows(
  observed_same_genus_pool,
  taxref_same_genus_pool
) %>%
  distinct(candidate_key, candidate_label, genus) %>%
  filter(!is.na(candidate_key), !is.na(candidate_label), !is.na(genus)) %>%
  arrange(genus, candidate_label)

readr::write_csv(candidate_taxonomic_pool, file.path(GBIF_DIR, "candidate_taxonomic_pool.csv"))

message("Nombre de candidats taxonomiques : ", nrow(candidate_taxonomic_pool))
message("Nombre de genres candidats : ", n_distinct(candidate_taxonomic_pool$genus))

# ------------------------------------------------------------
# 3. APPARIEMENT DES NOMS AVEC GBIF BACKBONE
# ------------------------------------------------------------

gbif_name_match_path <- file.path(GBIF_DIR, "gbif_name_matching.csv")

if (GBIF_USE_CACHE && file.exists(gbif_name_match_path)) {
  gbif_name_matching <- readr::read_csv(gbif_name_match_path, show_col_types = FALSE)
  message("Name matching GBIF recharge depuis : ", gbif_name_match_path)
} else {
  message_header("Appariement des noms avec GBIF Backbone")
  
  gbif_match_one <- function(candidate_key, candidate_label, genus_expected) {
    res <- tryCatch(
      rgbif::name_backbone(name = candidate_label, rank = "species"),
      error = function(e) NULL
    )
    
    if (is.null(res) || is.null(res$usageKey)) {
      return(tibble(
        candidate_key = candidate_key,
        candidate_label = candidate_label,
        genus = genus_expected,
        gbif_usage_key = NA_integer_,
        gbif_raw_usage_key = NA_integer_,
        gbif_accepted_usage_key = NA_integer_,
        gbif_scientific_name = NA_character_,
        gbif_canonical_name = NA_character_,
        gbif_rank = NA_character_,
        gbif_status = NA_character_,
        gbif_match_type = NA_character_,
        gbif_confidence = NA_real_,
        gbif_genus = NA_character_,
        gbif_match_ok = FALSE
      ))
    }
    
    raw_key <- safe_scalar_int(res$usageKey)
    accepted_key <- safe_scalar_int(res$acceptedUsageKey)
    usage_key <- accepted_key %||% raw_key
    gbif_genus <- safe_scalar_chr(res$genus)
    match_type <- safe_scalar_chr(res$matchType)
    confidence <- suppressWarnings(as.numeric(res$confidence %||% NA_real_))
    
    tibble(
      candidate_key = candidate_key,
      candidate_label = candidate_label,
      genus = genus_expected,
      gbif_usage_key = usage_key,
      gbif_raw_usage_key = raw_key,
      gbif_accepted_usage_key = accepted_key,
      gbif_scientific_name = safe_scalar_chr(res$scientificName),
      gbif_canonical_name = safe_scalar_chr(res$canonicalName),
      gbif_rank = safe_scalar_chr(res$rank),
      gbif_status = safe_scalar_chr(res$status),
      gbif_match_type = match_type,
      gbif_confidence = confidence,
      gbif_genus = gbif_genus,
      gbif_match_ok = !is.na(usage_key) &
        !is.na(gbif_genus) &
        gbif_genus == genus_expected &
        match_type %in% c("EXACT", "FUZZY")
    )
  }
  
  gbif_name_matching <- purrr::pmap_dfr(
    candidate_taxonomic_pool %>% select(candidate_key, candidate_label, genus),
    gbif_match_one
  )
  
  readr::write_csv(gbif_name_matching, gbif_name_match_path)
}

matched_keys <- gbif_name_matching %>%
  filter(gbif_match_ok, !is.na(gbif_usage_key)) %>%
  distinct(gbif_usage_key) %>%
  pull(gbif_usage_key)

message("Noms apparies GBIF retenus : ", length(matched_keys), " / ", nrow(candidate_taxonomic_pool))

if (length(matched_keys) == 0) {
  stop("Aucun taxon GBIF apparié correctement. Verifie les noms candidats et gbif_name_matching.csv.")
}

# ------------------------------------------------------------
# 4. RECUPERATION DES OCCURRENCES GBIF
# ------------------------------------------------------------

gbif_occ_raw_path <- file.path(GBIF_DIR, "gbif_occurrences_raw.rds")

gbif_user <- Sys.getenv("GBIF_USER", unset = NA_character_)
gbif_pwd <- Sys.getenv("GBIF_PWD", unset = NA_character_)
gbif_email <- Sys.getenv("GBIF_EMAIL", unset = NA_character_)
gbif_credentials_ok <- all(!is.na(c(gbif_user, gbif_pwd, gbif_email))) &&
  all(nchar(c(gbif_user, gbif_pwd, gbif_email)) > 0)

if (GBIF_MODE == "auto") {
  GBIF_MODE_USED <- if (gbif_credentials_ok) "download" else "search"
} else {
  GBIF_MODE_USED <- match.arg(GBIF_MODE, choices = c("download", "search"))
}

if (GBIF_USE_CACHE && file.exists(gbif_occ_raw_path)) {
  gbif_occ_raw <- readRDS(gbif_occ_raw_path)
  message("Occurrences GBIF brutes rechargees depuis : ", gbif_occ_raw_path)
} else {
  message_header(paste0("Recuperation des occurrences GBIF — mode ", GBIF_MODE_USED))
  
  if (GBIF_MODE_USED == "download") {
    if (!gbif_credentials_ok) {
      stop(
        "Mode download demande mais identifiants GBIF absents. ",
        "Definis GBIF_USER, GBIF_PWD et GBIF_EMAIL dans ~/.Renviron, ",
        "ou utilise options(taxo.gbif_mode = 'search') pour un test."
      )
    }
    
    dl_key <- rgbif::occ_download(
      rgbif::pred_in("taxonKey", matched_keys),
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      rgbif::pred("occurrenceStatus", "PRESENT"),
      rgbif::pred_in("country", GBIF_COUNTRIES),
      format = "SIMPLE_CSV",
      user = gbif_user,
      pwd = gbif_pwd,
      email = gbif_email
    )
    
    message("GBIF download key : ", dl_key)
    message("Attente de la preparation du download GBIF...")
    rgbif::occ_download_wait(dl_key)
    
    dl_file <- rgbif::occ_download_get(dl_key, path = GBIF_DIR, overwrite = FALSE)
    gbif_occ_raw <- rgbif::occ_download_import(dl_file)
    
    writeLines(as.character(dl_key), file.path(GBIF_DIR, "gbif_download_key.txt"))
  }
  
  if (GBIF_MODE_USED == "search") {
    warning(
      "Mode occ_search utilise. C'est acceptable pour tester le pipeline, ",
      "mais occ_download est preferable pour l'analyse finale et la citation GBIF."
    )
    
    gbif_occ_raw <- purrr::map_dfr(seq_along(matched_keys), function(i) {
      k <- matched_keys[i]
      message("occ_search ", i, "/", length(matched_keys), " — taxonKey ", k)
      
      res <- tryCatch(
        rgbif::occ_search(
          taxonKey = k,
          hasCoordinate = TRUE,
          hasGeospatialIssue = FALSE,
          occurrenceStatus = "PRESENT",
          limit = GBIF_SEARCH_LIMIT_PER_SPECIES
        ),
        error = function(e) NULL
      )
      
      if (is.null(res) || is.null(res$data) || nrow(res$data) == 0) {
        return(tibble())
      }
      
      as_tibble(res$data) %>% mutate(query_taxon_key = k)
    })
  }
  
  saveRDS(gbif_occ_raw, gbif_occ_raw_path)
}

message("Occurrences GBIF brutes : ", nrow(gbif_occ_raw))

if (nrow(gbif_occ_raw) == 0) {
  stop("Aucune occurrence GBIF récupérée.")
}

# ------------------------------------------------------------
# 5. NETTOYAGE DES OCCURRENCES
# ------------------------------------------------------------

message_header("Nettoyage des occurrences GBIF")

gbif_occ <- gbif_occ_raw %>%
  as_tibble() %>%
  janitor::clean_names()

# Colonnes robustes selon occ_download ou occ_search.
if (!"decimal_latitude" %in% names(gbif_occ) && "decimallatitude" %in% names(gbif_occ)) {
  gbif_occ <- gbif_occ %>% rename(decimal_latitude = decimallatitude)
}
if (!"decimal_longitude" %in% names(gbif_occ) && "decimallongitude" %in% names(gbif_occ)) {
  gbif_occ <- gbif_occ %>% rename(decimal_longitude = decimallongitude)
}
if (!"coordinate_uncertainty_in_meters" %in% names(gbif_occ)) {
  gbif_occ$coordinate_uncertainty_in_meters <- NA_real_
}
if (!"country_code" %in% names(gbif_occ)) {
  if ("country" %in% names(gbif_occ)) {
    gbif_occ$country_code <- gbif_occ$country
  } else {
    gbif_occ$country_code <- NA_character_
  }
}
if (!"basis_of_record" %in% names(gbif_occ)) {
  gbif_occ$basis_of_record <- NA_character_
}
if (!"scientific_name" %in% names(gbif_occ)) {
  gbif_occ$scientific_name <- NA_character_
}
if (!"species" %in% names(gbif_occ)) {
  gbif_occ$species <- gbif_occ$scientific_name
}

# Cle taxonomique GBIF robuste.
gbif_occ$gbif_usage_key <- NA_integer_
for (cc in c("accepted_taxon_key", "taxon_key", "species_key", "query_taxon_key")) {
  if (cc %in% names(gbif_occ)) {
    gbif_occ$gbif_usage_key <- dplyr::coalesce(
      gbif_occ$gbif_usage_key,
      suppressWarnings(as.integer(gbif_occ[[cc]]))
    )
  }
}

gbif_clean <- gbif_occ %>%
  mutate(
    decimal_latitude = safe_numeric(decimal_latitude),
    decimal_longitude = safe_numeric(decimal_longitude),
    coordinate_uncertainty_in_meters = safe_numeric(coordinate_uncertainty_in_meters),
    country_code = toupper(as.character(country_code)),
    basis_of_record = toupper(as.character(basis_of_record)),
    species_for_cleaning = coalesce(as.character(species), as.character(scientific_name), as.character(gbif_usage_key))
  ) %>%
  filter(
    !is.na(gbif_usage_key),
    gbif_usage_key %in% matched_keys,
    is.finite(decimal_latitude),
    is.finite(decimal_longitude),
    decimal_latitude >= -90,
    decimal_latitude <= 90,
    decimal_longitude >= -180,
    decimal_longitude <= 180,
    is.na(coordinate_uncertainty_in_meters) |
      coordinate_uncertainty_in_meters <= GBIF_MAX_COORD_UNCERTAINTY_M,
    is.na(country_code) | country_code %in% GBIF_COUNTRIES,
    !basis_of_record %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")
  )

message("Occurrences apres filtres simples : ", nrow(gbif_clean))

# Nettoyage spatial optionnel avec CoordinateCleaner.
if (requireNamespace("CoordinateCleaner", quietly = TRUE) && nrow(gbif_clean) > 0) {
  cc_flags <- tryCatch(
    CoordinateCleaner::clean_coordinates(
      x = gbif_clean,
      lon = "decimal_longitude",
      lat = "decimal_latitude",
      species = "species_for_cleaning",
      tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
      value = "spatialvalid"
    ),
    error = function(e) {
      warning("CoordinateCleaner a echoue : ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(cc_flags)) {
    if (is.data.frame(cc_flags) && ".summary" %in% names(cc_flags)) {
      gbif_clean <- gbif_clean[cc_flags$.summary, , drop = FALSE]
    } else if (is.logical(cc_flags) && length(cc_flags) == nrow(gbif_clean)) {
      gbif_clean <- gbif_clean[cc_flags, , drop = FALSE]
    }
  }
} else {
  warning("Package CoordinateCleaner absent : nettoyage spatial avance non applique.")
}

readr::write_csv(gbif_clean, file.path(GBIF_DIR, "gbif_occurrences_clean.csv"))
message("Occurrences apres nettoyage spatial : ", nrow(gbif_clean))

if (nrow(gbif_clean) == 0) {
  stop("Toutes les occurrences GBIF ont ete filtrees. Relache les filtres ou verifie le matching.")
}

# ------------------------------------------------------------
# 6. AOO OBSERVEE GBIF ET OCCURRENCES REGIONALES
# ------------------------------------------------------------

message_header("Calcul AOO observee et filtre regional")

occ_sf_3035 <- gbif_clean %>%
  st_as_sf(coords = c("decimal_longitude", "decimal_latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(GBIF_CRS_METRIC)

xy <- st_coordinates(occ_sf_3035)
occ_sf_3035$grid_x_2km <- floor(xy[, "X"] / GBIF_AOO_CELLSIZE_M)
occ_sf_3035$grid_y_2km <- floor(xy[, "Y"] / GBIF_AOO_CELLSIZE_M)
occ_sf_3035$aoo_cell_id_2km <- paste(occ_sf_3035$grid_x_2km, occ_sf_3035$grid_y_2km, sep = ":")

gbif_aoo_by_species <- occ_sf_3035 %>%
  st_drop_geometry() %>%
  group_by(gbif_usage_key) %>%
  summarise(
    n_occ_clean = n(),
    n_countries = n_distinct(country_code, na.rm = TRUE),
    countries = paste(sort(unique(na.omit(country_code))), collapse = ";"),
    n_aoo_cells_2km = n_distinct(aoo_cell_id_2km),
    aoo_km2 = n_aoo_cells_2km * (GBIF_AOO_CELLSIZE_M / 1000)^2,
    .groups = "drop"
  )

readr::write_csv(gbif_aoo_by_species, file.path(GBIF_DIR, "gbif_aoo_by_species.csv"))

# Coordonnees des stations.
if (exists("meta_station") && all(c("station", "longitude", "latitude") %in% names(meta_station))) {
  station_coords <- meta_station %>%
    transmute(
      station = as.character(station),
      longitude = safe_numeric(longitude),
      latitude = safe_numeric(latitude)
    )
} else if (all(c("station", "x_longitude", "y_latitude") %in% names(col_used))) {
  station_coords <- col_used %>%
    group_by(station) %>%
    summarise(
      longitude = safe_numeric(first_non_missing(x_longitude)),
      latitude = safe_numeric(first_non_missing(y_latitude)),
      .groups = "drop"
    )
} else {
  stop(
    "Impossible de trouver les coordonnees station. ",
    "Il faut meta_station$longitude/$latitude ou col_used$x_longitude/$y_latitude."
  )
}

station_coords <- station_coords %>%
  filter(is.finite(longitude), is.finite(latitude)) %>%
  distinct(station, .keep_all = TRUE)

if (nrow(station_coords) == 0) {
  stop("Aucune station avec coordonnees valides.")
}

stations_sf_3035 <- station_coords %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(GBIF_CRS_METRIC)

station_buffers <- stations_sf_3035 %>%
  select(station) %>%
  st_buffer(dist = GBIF_BUFFER_KM * 1000)

regional_occ_by_station <- st_join(
  station_buffers,
  occ_sf_3035 %>% select(gbif_usage_key, country_code, aoo_cell_id_2km),
  join = st_intersects,
  left = FALSE
) %>%
  st_drop_geometry() %>%
  group_by(station, gbif_usage_key) %>%
  summarise(
    n_occ_buffer = n(),
    n_aoo_cells_buffer_2km = n_distinct(aoo_cell_id_2km),
    buffer_aoo_km2 = n_aoo_cells_buffer_2km * (GBIF_AOO_CELLSIZE_M / 1000)^2,
    .groups = "drop"
  )

readr::write_csv(regional_occ_by_station, file.path(GBIF_DIR, "gbif_occurrences_by_station_buffer.csv"))

# ------------------------------------------------------------
# 7. POOL DE CONFUSION STATION-SPECIFIQUE
# ------------------------------------------------------------

message_header("Construction du pool de confusion biogeographique")

candidate_gbif <- candidate_taxonomic_pool %>%
  left_join(
    gbif_name_matching %>%
      filter(gbif_match_ok, !is.na(gbif_usage_key)) %>%
      select(candidate_key, candidate_label, genus, gbif_usage_key, gbif_match_type, gbif_confidence),
    by = c("candidate_key", "candidate_label", "genus")
  ) %>%
  left_join(gbif_aoo_by_species, by = "gbif_usage_key")

# Pool GBIF buffer stationnel.
station_pool_buffer <- regional_occ_by_station %>%
  inner_join(candidate_gbif, by = "gbif_usage_key") %>%
  mutate(source_buffer = TRUE)

# Pool observe dans ton propre dataset a la station : garde-fou contre les lacunes GBIF.
observed_station_species <- col_used %>%
  filter(taxo_level == "species", !is.na(species_key)) %>%
  transmute(
    station = as.character(station),
    candidate_key = species_key
  ) %>%
  distinct()

station_pool_observed_dataset <- observed_station_species %>%
  inner_join(candidate_gbif, by = "candidate_key") %>%
  mutate(
    n_occ_buffer = 0L,
    n_aoo_cells_buffer_2km = 0L,
    buffer_aoo_km2 = 0,
    source_observed_dataset_station = TRUE
  )

# Pool regional large : toutes les especes candidates connues dans GBIF_COUNTRIES.
candidate_region_all <- candidate_gbif %>%
  filter(!is.na(gbif_usage_key), !is.na(n_occ_clean))

station_pool_country_region <- tidyr::expand_grid(
  station = station_coords$station,
  candidate_row_id = seq_len(nrow(candidate_region_all))
) %>%
  left_join(
    candidate_region_all %>% mutate(candidate_row_id = row_number()),
    by = "candidate_row_id"
  ) %>%
  select(-candidate_row_id) %>%
  mutate(
    n_occ_buffer = 0L,
    n_aoo_cells_buffer_2km = 0L,
    buffer_aoo_km2 = 0,
    source_country_region = TRUE
  )

pool_buffer_or_observed <- bind_rows(
  station_pool_buffer %>% mutate(source_observed_dataset_station = FALSE, source_country_region = FALSE),
  station_pool_observed_dataset %>% mutate(source_buffer = FALSE, source_country_region = FALSE)
) %>%
  group_by(station, candidate_key, candidate_label, genus, gbif_usage_key) %>%
  summarise(
    n_occ_buffer = max(n_occ_buffer, na.rm = TRUE),
    n_aoo_cells_buffer_2km = max(n_aoo_cells_buffer_2km, na.rm = TRUE),
    buffer_aoo_km2 = max(buffer_aoo_km2, na.rm = TRUE),
    n_occ_clean = max(n_occ_clean, na.rm = TRUE),
    n_aoo_cells_2km = max(n_aoo_cells_2km, na.rm = TRUE),
    aoo_km2 = max(aoo_km2, na.rm = TRUE),
    countries = first_non_missing(countries),
    source_buffer = any(source_buffer, na.rm = TRUE),
    source_observed_dataset_station = any(source_observed_dataset_station, na.rm = TRUE),
    source_country_region = any(source_country_region, na.rm = TRUE),
    .groups = "drop"
  )

pool_country_region <- station_pool_country_region %>%
  group_by(station, candidate_key, candidate_label, genus, gbif_usage_key) %>%
  summarise(
    n_occ_buffer = max(n_occ_buffer, na.rm = TRUE),
    n_aoo_cells_buffer_2km = max(n_aoo_cells_buffer_2km, na.rm = TRUE),
    buffer_aoo_km2 = max(buffer_aoo_km2, na.rm = TRUE),
    n_occ_clean = max(n_occ_clean, na.rm = TRUE),
    n_aoo_cells_2km = max(n_aoo_cells_2km, na.rm = TRUE),
    aoo_km2 = max(aoo_km2, na.rm = TRUE),
    countries = first_non_missing(countries),
    source_buffer = FALSE,
    source_observed_dataset_station = FALSE,
    source_country_region = TRUE,
    .groups = "drop"
  )

if (GBIF_POOL_STRATEGY == "buffer_or_observed") {
  gbif_station_candidate_pool <- pool_buffer_or_observed
}

if (GBIF_POOL_STRATEGY == "country_region") {
  gbif_station_candidate_pool <- pool_country_region
}

if (GBIF_POOL_STRATEGY == "buffer_then_region_fallback") {
  station_genus_all <- tidyr::crossing(
    station = station_coords$station,
    genus = candidate_region_all %>% pull(genus) %>% unique()
  )
  
  buffer_counts <- pool_buffer_or_observed %>%
    count(station, genus, name = "n_candidates_buffer")
  
  fallback_station_genus <- station_genus_all %>%
    left_join(buffer_counts, by = c("station", "genus")) %>%
    mutate(n_candidates_buffer = replace_na(n_candidates_buffer, 0L)) %>%
    filter(n_candidates_buffer < 2)
  
  gbif_station_candidate_pool <- bind_rows(
    pool_buffer_or_observed,
    pool_country_region %>%
      semi_join(fallback_station_genus, by = c("station", "genus"))
  ) %>%
    group_by(station, candidate_key, candidate_label, genus, gbif_usage_key) %>%
    summarise(
      n_occ_buffer = max(n_occ_buffer, na.rm = TRUE),
      n_aoo_cells_buffer_2km = max(n_aoo_cells_buffer_2km, na.rm = TRUE),
      buffer_aoo_km2 = max(buffer_aoo_km2, na.rm = TRUE),
      n_occ_clean = max(n_occ_clean, na.rm = TRUE),
      n_aoo_cells_2km = max(n_aoo_cells_2km, na.rm = TRUE),
      aoo_km2 = max(aoo_km2, na.rm = TRUE),
      countries = first_non_missing(countries),
      source_buffer = any(source_buffer, na.rm = TRUE),
      source_observed_dataset_station = any(source_observed_dataset_station, na.rm = TRUE),
      source_country_region = any(source_country_region, na.rm = TRUE),
      .groups = "drop"
    )
}

# Poids de tirage : donne plus de poids aux especes observees dans le dataset/station,
# puis aux especes avec occurrences GBIF dans le buffer, puis au fallback regional.
gbif_station_candidate_pool <- gbif_station_candidate_pool %>%
  mutate(
    n_occ_buffer = if_else(is.infinite(n_occ_buffer), 0, n_occ_buffer),
    n_occ_clean = if_else(is.infinite(n_occ_clean), 0, n_occ_clean),
    candidate_weight_raw = case_when(
      source_observed_dataset_station ~ 5 + log1p(n_occ_clean),
      source_buffer & n_occ_buffer > 0 ~ 2 + log1p(n_occ_buffer),
      source_country_region ~ 0.5 + log1p(n_occ_clean),
      TRUE ~ 1
    ),
    candidate_weight = pmax(candidate_weight_raw, 0.001)
  )

# Evite qu'une espece hyper-documentee domine tous les tirages.
w_cap <- quantile(gbif_station_candidate_pool$candidate_weight, 0.95, na.rm = TRUE, names = FALSE)
gbif_station_candidate_pool <- gbif_station_candidate_pool %>%
  mutate(candidate_weight = pmin(candidate_weight, w_cap))

readr::write_csv(gbif_station_candidate_pool, file.path(GBIF_DIR, "gbif_station_candidate_pool.csv"))

pool_diagnostic <- gbif_station_candidate_pool %>%
  group_by(station, genus) %>%
  summarise(
    n_candidates = n_distinct(candidate_key),
    n_candidates_buffer = n_distinct(candidate_key[source_buffer]),
    n_candidates_observed_dataset = n_distinct(candidate_key[source_observed_dataset_station]),
    n_candidates_country_region = n_distinct(candidate_key[source_country_region]),
    .groups = "drop"
  )

readr::write_csv(pool_diagnostic, file.path(GBIF_DIR, "gbif_station_genus_pool_diagnostic.csv"))

message("Candidats station-specifiques : ", nrow(gbif_station_candidate_pool))
message("Combinaisons station x genre avec >=2 candidats : ", sum(pool_diagnostic$n_candidates >= 2))
message("Strategie de pool : ", GBIF_POOL_STRATEGY)

# ------------------------------------------------------------
# 8. SHUFFLING GBIF REGIONAL, GRADIENT 1-20%
# ------------------------------------------------------------

message_header("Simulation du gradient d'erreur GBIF regional")

base_cols <- c("station", "repetition", "sample_id", "abundance", "taxo_level", "binomial", "species_key", "species_label", "genus", "famille", "ordre")

species_records <- col_used %>%
  filter(taxo_level == "species", !is.na(species_key), !is.na(genus)) %>%
  select(all_of(base_cols))

choose_gbif_station_candidate <- function(station_i, genus_i, current_species_key, pool_tbl) {
  cand <- pool_tbl %>%
    filter(
      .data$station == !!as.character(station_i),
      .data$genus == !!genus_i,
      .data$candidate_key != !!current_species_key
    ) %>%
    distinct(candidate_key, candidate_label, genus, candidate_weight, .keep_all = TRUE)
  
  if (nrow(cand) == 0) {
    return(tibble(candidate_key = NA_character_, candidate_label = NA_character_, genus = NA_character_))
  }
  
  if (isTRUE(GBIF_WEIGHTED_SHUFFLING)) {
    cand %>% slice_sample(n = 1, weight_by = candidate_weight)
  } else {
    cand %>% slice_sample(n = 1)
  }
}

simulate_gbif_regional_error <- function(dat_species, pool_tbl, scenario_name, iter_id, fixed_rate) {
  dat0 <- dat_species %>% mutate(p_error = fixed_rate)
  
  purrr::pmap_dfr(dat0, function(station, repetition, sample_id, abundance, taxo_level,
                                 binomial, species_key, species_label, genus, famille, ordre, p_error) {
    n <- as.integer(round(abundance))
    
    no_swap_row <- tibble(
      station = station, repetition = repetition, sample_id = sample_id,
      abundance = abundance, taxo_level = taxo_level, binomial = binomial,
      species_key = species_key, species_label = species_label,
      genus = genus, famille = famille, ordre = ordre,
      source_species_key = species_key,
      taxon_unit = paste0("species:", species_key),
      is_swapped = FALSE
    )
    
    if (is.na(p_error) || p_error <= 0 || n <= 0) return(no_swap_row)
    
    target <- choose_gbif_station_candidate(station, genus, species_key, pool_tbl)
    if (nrow(target) == 0 || is.na(target$candidate_key[1])) return(no_swap_row)
    
    n_swap <- rbinom(1, size = n, prob = p_error)
    if (n_swap == 0) return(no_swap_row)
    
    target_key <- target$candidate_key[1]
    target_label <- target$candidate_label[1]
    target_genus <- target$genus[1]
    
    bind_rows(
      if (n - n_swap > 0) {
        tibble(
          station = station, repetition = repetition, sample_id = sample_id,
          abundance = n - n_swap, taxo_level = taxo_level, binomial = binomial,
          species_key = species_key, species_label = species_label,
          genus = genus, famille = famille, ordre = ordre,
          source_species_key = species_key,
          taxon_unit = paste0("species:", species_key),
          is_swapped = FALSE
        )
      },
      tibble(
        station = station, repetition = repetition, sample_id = sample_id,
        abundance = n_swap, taxo_level = taxo_level, binomial = target_label,
        species_key = target_key, species_label = target_label,
        genus = target_genus, famille = famille, ordre = ordre,
        source_species_key = species_key,
        taxon_unit = paste0("species:", target_key),
        is_swapped = TRUE
      )
    )
  }) %>%
    mutate(
      scenario = scenario_name,
      scenario_family = "gbif_regional_error_gradient",
      baseline_scenario = "gbif_species_baseline",
      iter = iter_id,
      unit_type = "species_units"
    )
}

gbif_error_grid <- tibble(
  error_rate = GBIF_ERROR_RATES,
  error_pct = round(100 * error_rate),
  scenario = paste0("gbif_regional_same_genus_", sprintf("%02d", error_pct), "pct")
)

scenario_baseline_gbif <- species_records %>%
  mutate(
    source_species_key = species_key,
    taxon_unit = paste0("species:", species_key),
    is_swapped = FALSE,
    scenario = "gbif_species_baseline",
    scenario_family = "gbif_regional_error_gradient",
    baseline_scenario = "gbif_species_baseline",
    iter = 1,
    unit_type = "species_units",
    error_rate = 0,
    error_pct = 0
  )

scenario_gradient_gbif <- purrr::map_dfr(seq_len(GBIF_N_SIM), function(iter_i) {
  set.seed(12000 + iter_i)
  purrr::pmap_dfr(gbif_error_grid, function(error_rate, error_pct, scenario) {
    simulate_gbif_regional_error(
      dat_species = species_records,
      pool_tbl = gbif_station_candidate_pool,
      scenario_name = scenario,
      iter_id = iter_i,
      fixed_rate = error_rate
    ) %>%
      mutate(error_rate = error_rate, error_pct = error_pct)
  })
})

gbif_scenario_data <- bind_rows(scenario_baseline_gbif, scenario_gradient_gbif) %>%
  mutate(
    pool_strategy = GBIF_POOL_STRATEGY,
    buffer_km = GBIF_BUFFER_KM,
    weighted_shuffling = GBIF_WEIGHTED_SHUFFLING
  )

saveRDS(gbif_scenario_data, file.path(GBIF_DIR, "gbif_scenario_data.rds"))

readr::write_csv(
  gbif_scenario_data %>%
    distinct(scenario, scenario_family, baseline_scenario, unit_type, error_rate, error_pct, pool_strategy, buffer_km, weighted_shuffling),
  file.path(GBIF_DIR, "gbif_scenario_definitions.csv")
)

gbif_effective_error_by_iter <- gbif_scenario_data %>%
  group_by(scenario, iter, error_rate, error_pct) %>%
  summarise(
    total_individuals = sum(abundance, na.rm = TRUE),
    swapped_individuals = sum(abundance[is_swapped], na.rm = TRUE),
    effective_error_rate = swapped_individuals / total_individuals,
    .groups = "drop"
  )

readr::write_csv(gbif_effective_error_by_iter, file.path(GBIF_DIR, "gbif_effective_error_by_iter.csv"))

# ------------------------------------------------------------
# 9. METRIQUES BIODIVERSITE ET STABILITE
# ------------------------------------------------------------

message_header("Calcul des metriques biodiversite et stabilite")

gbif_alpha_diversity_by_unit <- gbif_scenario_data %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct) %>%
  group_split() %>%
  purrr::map_dfr(function(df) {
    mat <- make_comm_matrix(df, unit_col = ANALYSIS_UNIT)
    hill_metrics_from_matrix(mat) %>%
      mutate(
        scenario_family = df$scenario_family[1],
        scenario = df$scenario[1],
        baseline_scenario = df$baseline_scenario[1],
        unit_type = df$unit_type[1],
        iter = df$iter[1],
        error_rate = df$error_rate[1],
        error_pct = df$error_pct[1]
      )
  }) %>%
  relocate(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct, unit)

readr::write_csv(gbif_alpha_diversity_by_unit, file.path(GBIF_DIR, "gbif_alpha_diversity_by_unit.csv"))

gbif_gamma_by_scenario_iter <- gbif_scenario_data %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct) %>%
  summarise(gamma_taxon_units = n_distinct(taxon_unit), .groups = "drop")

gbif_scenario_summary_by_iter <- gbif_alpha_diversity_by_unit %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct) %>%
  summarise(
    n_units = n(),
    total_abundance = sum(total_abundance, na.rm = TRUE),
    mean_local_q0 = mean(q0, na.rm = TRUE),
    mean_local_q1 = mean(q1, na.rm = TRUE),
    mean_local_q2 = mean(q2, na.rm = TRUE),
    mean_coverage_chao = mean(coverage_chao, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    gbif_gamma_by_scenario_iter,
    by = c("scenario_family", "scenario", "baseline_scenario", "unit_type", "iter", "error_rate", "error_pct")
  ) %>%
  left_join(
    gbif_effective_error_by_iter %>% select(scenario, iter, effective_error_rate),
    by = c("scenario", "iter")
  )

readr::write_csv(gbif_scenario_summary_by_iter, file.path(GBIF_DIR, "gbif_scenario_summary_by_iter.csv"))

gbif_baseline_metrics <- gbif_scenario_summary_by_iter %>%
  filter(scenario == "gbif_species_baseline") %>%
  summarise(
    mean_local_q0 = mean(mean_local_q0, na.rm = TRUE),
    mean_local_q1 = mean(mean_local_q1, na.rm = TRUE),
    mean_local_q2 = mean(mean_local_q2, na.rm = TRUE),
    mean_coverage_chao = mean(mean_coverage_chao, na.rm = TRUE),
    gamma_taxon_units = mean(gamma_taxon_units, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "baseline_value")

gbif_metric_key <- tribble(
  ~metric, ~metric_label,
  "mean_local_q0", "Alpha q0 richness",
  "mean_local_q1", "Alpha Hill q1",
  "mean_local_q2", "Alpha Hill q2",
  "mean_coverage_chao", "Mean sample coverage",
  "gamma_taxon_units", "Gamma richness"
)

gbif_biodiversity_curves <- gbif_scenario_summary_by_iter %>%
  filter(scenario != "gbif_species_baseline") %>%
  select(error_rate, error_pct, effective_error_rate, iter,
         mean_local_q0, mean_local_q1, mean_local_q2, mean_coverage_chao, gamma_taxon_units) %>%
  pivot_longer(
    cols = c(mean_local_q0, mean_local_q1, mean_local_q2, mean_coverage_chao, gamma_taxon_units),
    names_to = "metric",
    values_to = "value"
  ) %>%
  left_join(gbif_baseline_metrics, by = "metric") %>%
  left_join(gbif_metric_key, by = "metric") %>%
  mutate(
    response_ratio = value / baseline_value,
    percent_change = 100 * (response_ratio - 1),
    effective_error_pct = 100 * effective_error_rate
  )

gbif_biodiversity_curve_summary <- gbif_biodiversity_curves %>%
  group_by(metric, metric_label, error_rate, error_pct) %>%
  summarise(
    mean_effective_error_pct = mean(effective_error_pct, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE),
    p10_value = quantile(value, 0.10, na.rm = TRUE, names = FALSE),
    p90_value = quantile(value, 0.90, na.rm = TRUE, names = FALSE),
    mean_percent_change = mean(percent_change, na.rm = TRUE),
    p10_percent_change = quantile(percent_change, 0.10, na.rm = TRUE, names = FALSE),
    p90_percent_change = quantile(percent_change, 0.90, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

readr::write_csv(gbif_biodiversity_curve_summary, file.path(GBIF_DIR, "gbif_biodiversity_curve_summary.csv"))

# Stabilite vs baseline.
gbif_matrix_index <- gbif_scenario_data %>%
  distinct(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct) %>%
  arrange(error_rate, iter)

make_gbif_matrix_for_key <- function(scen, iter_i) {
  gbif_scenario_data %>%
    filter(scenario == scen, iter == iter_i) %>%
    make_comm_matrix(unit_col = ANALYSIS_UNIT)
}

gbif_stability_by_iter <- gbif_matrix_index %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct) %>%
  group_split() %>%
  purrr::map_dfr(function(key) {
    scen <- key$scenario[1]
    iter_i <- key$iter[1]
    mat_base <- make_gbif_matrix_for_key("gbif_species_baseline", 1)
    mat_scen <- make_gbif_matrix_for_key(scen, iter_i)
    
    stability_one(mat_base, mat_scen) %>%
      mutate(
        scenario_family = key$scenario_family[1],
        scenario = scen,
        baseline_scenario = key$baseline_scenario[1],
        unit_type = key$unit_type[1],
        iter = iter_i,
        error_rate = key$error_rate[1],
        error_pct = key$error_pct[1]
      )
  }) %>%
  relocate(scenario_family, scenario, baseline_scenario, unit_type, iter, error_rate, error_pct)

readr::write_csv(gbif_stability_by_iter, file.path(GBIF_DIR, "gbif_stability_by_iter.csv"))

gbif_stability_long <- gbif_stability_by_iter %>%
  filter(scenario != "gbif_species_baseline") %>%
  select(error_rate, error_pct, iter,
         q0_cor, q1_cor, q2_cor, coverage_chao_cor,
         bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2) %>%
  pivot_longer(
    cols = c(q0_cor, q1_cor, q2_cor, coverage_chao_cor,
             bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2),
    names_to = "metric",
    values_to = "stability"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      q0_cor = "Alpha q0 richness",
      q1_cor = "Alpha Hill q1",
      q2_cor = "Alpha Hill q2",
      coverage_chao_cor = "Mean sample coverage",
      bray_curtis_cor = "Bray-Curtis",
      jaccard_pa_cor = "Jaccard",
      sorensen_pa_cor = "Sørensen",
      procrustes_r2 = "Ordination (Procrustes R2)"
    )
  ) %>%
  left_join(
    gbif_effective_error_by_iter %>% mutate(effective_error_pct = 100 * effective_error_rate) %>% select(scenario, iter, effective_error_pct),
    by = c("scenario", "iter")
  )

gbif_stability_curve_summary <- gbif_stability_long %>%
  group_by(metric, metric_label, error_rate, error_pct) %>%
  summarise(
    mean_effective_error_pct = mean(effective_error_pct, na.rm = TRUE),
    mean_stability = mean(stability, na.rm = TRUE),
    p10_stability = quantile(stability, 0.10, na.rm = TRUE, names = FALSE),
    p90_stability = quantile(stability, 0.90, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

readr::write_csv(gbif_stability_curve_summary, file.path(GBIF_DIR, "gbif_stability_curve_summary.csv"))

# Seuils critiques.
gbif_threshold_biodiv <- gbif_biodiversity_curve_summary %>%
  group_by(metric, metric_label) %>%
  summarise(
    threshold_type = "relative_change",
    critical_definition = paste0("|relative change| >= ", round(100 * GBIF_CRITICAL_REL_CHANGE), "%"),
    critical_error_pct = suppressWarnings(min(error_pct[abs(mean_percent_change) >= 100 * GBIF_CRITICAL_REL_CHANGE], na.rm = TRUE)),
    critical_effective_error_pct = suppressWarnings(min(mean_effective_error_pct[abs(mean_percent_change) >= 100 * GBIF_CRITICAL_REL_CHANGE], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    critical_error_pct = ifelse(is.infinite(critical_error_pct), NA_real_, critical_error_pct),
    critical_effective_error_pct = ifelse(is.infinite(critical_effective_error_pct), NA_real_, critical_effective_error_pct)
  )

gbif_threshold_stability <- gbif_stability_curve_summary %>%
  group_by(metric, metric_label) %>%
  summarise(
    threshold_type = "stability",
    critical_definition = paste0("mean stability <= ", GBIF_CRITICAL_STABILITY),
    critical_error_pct = suppressWarnings(min(error_pct[mean_stability <= GBIF_CRITICAL_STABILITY], na.rm = TRUE)),
    critical_effective_error_pct = suppressWarnings(min(mean_effective_error_pct[mean_stability <= GBIF_CRITICAL_STABILITY], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    critical_error_pct = ifelse(is.infinite(critical_error_pct), NA_real_, critical_error_pct),
    critical_effective_error_pct = ifelse(is.infinite(critical_effective_error_pct), NA_real_, critical_effective_error_pct)
  )

gbif_critical_thresholds <- bind_rows(gbif_threshold_biodiv, gbif_threshold_stability) %>%
  relocate(threshold_type, metric, metric_label, critical_definition, critical_error_pct, critical_effective_error_pct)

readr::write_csv(gbif_critical_thresholds, file.path(GBIF_DIR, "gbif_critical_thresholds.csv"))

# ------------------------------------------------------------
# 10. FIGURES
# ------------------------------------------------------------

message_header("Production des figures")

metric_levels_biodiv <- c(
  "Alpha q0 richness",
  "Alpha Hill q1",
  "Alpha Hill q2",
  "Mean sample coverage",
  "Gamma richness"
)

p_gbif_biodiv <- gbif_biodiversity_curve_summary %>%
  mutate(metric_label = factor(metric_label, levels = metric_levels_biodiv)) %>%
  ggplot(aes(x = error_pct, y = mean_percent_change)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = c(-100 * GBIF_CRITICAL_REL_CHANGE, 100 * GBIF_CRITICAL_REL_CHANGE),
             linewidth = 0.25, linetype = "dotted") +
  geom_ribbon(aes(ymin = p10_percent_change, ymax = p90_percent_change), alpha = 0.18) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  labs(
    title = "GBIF-informed regional shuffling: biodiversity response curves",
    subtitle = paste0(
      "Pool strategy: ", GBIF_POOL_STRATEGY,
      " | buffer: ", GBIF_BUFFER_KM, " km",
      " | dotted lines: ±", round(100 * GBIF_CRITICAL_REL_CHANGE), "%"
    ),
    x = "Requested identification error rate (%)",
    y = "Relative change vs species baseline (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A1_biodiversity_response_curves.png"), p_gbif_biodiv,
       width = 9.5, height = 7.2, dpi = 600)
ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A1_biodiversity_response_curves.pdf"), p_gbif_biodiv,
       width = 9.5, height = 7.2)

metric_levels_stability <- c(
  "Alpha q0 richness", "Alpha Hill q1", "Alpha Hill q2", "Mean sample coverage",
  "Bray-Curtis", "Jaccard", "Sørensen", "Ordination (Procrustes R2)"
)

p_gbif_stability <- gbif_stability_curve_summary %>%
  mutate(metric_label = factor(metric_label, levels = metric_levels_stability)) %>%
  ggplot(aes(x = error_pct, y = mean_stability)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = GBIF_CRITICAL_STABILITY, linewidth = 0.25, linetype = "dotted") +
  geom_ribbon(aes(ymin = p10_stability, ymax = p90_stability), alpha = 0.18) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  labs(
    title = "GBIF-informed regional shuffling: inference stability",
    subtitle = paste0("Dotted line: stability threshold ", GBIF_CRITICAL_STABILITY),
    x = "Requested identification error rate (%)",
    y = "Mean stability relative to species baseline"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A2_stability_response_curves.png"), p_gbif_stability,
       width = 9.5, height = 9.0, dpi = 600)
ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A2_stability_response_curves.pdf"), p_gbif_stability,
       width = 9.5, height = 9.0)

p_gbif_effective <- gbif_effective_error_by_iter %>%
  filter(error_pct > 0) %>%
  group_by(error_pct) %>%
  summarise(
    mean_effective_error_pct = 100 * mean(effective_error_rate, na.rm = TRUE),
    p10 = 100 * quantile(effective_error_rate, 0.10, na.rm = TRUE, names = FALSE),
    p90 = 100 * quantile(effective_error_rate, 0.90, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = error_pct, y = mean_effective_error_pct)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.18) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  labs(
    title = "Requested vs realised identification error",
    subtitle = "Realised error can be lower when no biogeographically plausible congeneric candidate is available",
    x = "Requested identification error rate (%)",
    y = "Realised shuffled abundance (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A3_requested_vs_realised_error.png"), p_gbif_effective,
       width = 7.2, height = 5.0, dpi = 600)
ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A3_requested_vs_realised_error.pdf"), p_gbif_effective,
       width = 7.2, height = 5.0)

p_gbif_thresholds <- gbif_critical_thresholds %>%
  mutate(
    threshold_label = recode(
      threshold_type,
      relative_change = "Relative change criterion",
      stability = "Stability criterion"
    ),
    metric_label = forcats::fct_reorder(metric_label, critical_error_pct, .na_rm = FALSE, .desc = TRUE)
  ) %>%
  ggplot(aes(x = critical_error_pct, y = metric_label)) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ threshold_label, scales = "free_x") +
  labs(
    title = "Critical error rate under GBIF-informed regional shuffling",
    subtitle = "First requested error level reaching the chosen criticality criterion",
    x = "Critical requested error rate (%)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A4_critical_error_thresholds.png"), p_gbif_thresholds,
       width = 9.0, height = 5.0, dpi = 600)
ggsave(file.path(GBIF_FIG_DIR, "Fig_GBIF_A4_critical_error_thresholds.pdf"), p_gbif_thresholds,
       width = 9.0, height = 5.0)

# ------------------------------------------------------------
# 11. RAPPORT COURT CONSOLE
# ------------------------------------------------------------

message_header("Termine")
message("Dossier de sortie : ", GBIF_DIR)
message("Dossier figures : ", GBIF_FIG_DIR)
message("Fichiers importants :")
message(" - gbif_name_matching.csv")
message(" - gbif_occurrences_clean.csv")
message(" - gbif_aoo_by_species.csv")
message(" - gbif_station_candidate_pool.csv")
message(" - gbif_station_genus_pool_diagnostic.csv")
message(" - gbif_biodiversity_curve_summary.csv")
message(" - gbif_stability_curve_summary.csv")
message(" - gbif_critical_thresholds.csv")
message(" - Fig_GBIF_A1_biodiversity_response_curves.*")
message(" - Fig_GBIF_A2_stability_response_curves.*")
message(" - Fig_GBIF_A3_requested_vs_realised_error.*")
message(" - Fig_GBIF_A4_critical_error_thresholds.*")
