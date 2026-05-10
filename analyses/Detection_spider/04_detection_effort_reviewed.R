# ============================================================
# DÉTECTIONS PAR RÉPLICAT & MATRICES D’EFFORT — REVIEWED v2
# Focus recommandé pour le manuscrit: Pitfall + GPD
#
# Sorties: det_wide, det_cols, Y, eff_mat, eff_z, site_cov, sites, spp, Rnames
# + tables reviewer-ready dans reviewer_response/
#
# Principes clés:
# 1) Pitfall et GPD sont déployés sur tous les sites, mais l'effort effectif
#    varie selon les pièges effectivement récupérés/exploitables.
# 2) Une non-détection dans un piège/protocole récupéré = 0.
# 3) Une visite non exploitable / aucun piège récupéré = NA, pas 0.
# 4) GPD = un device, avec effort = longueur de barrière + ouvertures des pots
#    effectivement récupérés.
# 5) Pitfall2/4/6/8/10 = sous-ensembles nominaux imbriqués; l'effort utilise
#    le nombre effectif de pots récupérés dans le sous-ensemble.
# ============================================================

# ---------- 0. Paramètres locaux avec valeurs par défaut ----------
if (!exists("pit_reps")) pit_reps <- c(10, 8, 6, 4, 2)
if (!exists("pit_circ_m")) pit_circ_m <- pi * 0.05
if (!exists("gpd_fence_length_m")) gpd_fence_length_m <- 1.0
if (!exists("require_complete_pitfall_subsets")) require_complete_pitfall_subsets <- FALSE
if (!exists("reviewer_out_dir")) reviewer_out_dir <- file.path(out_dir, "reviewer_response")
dir.create(reviewer_out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. Utilitaires ----------
med_days <- function(x) {
  z <- suppressWarnings(stats::median(x, na.rm = TRUE))
  if (!is.finite(z) || is.na(z) || z < 1) 1 else z
}

# Sécurise les colonnes attendues
stopifnot(all(c("site_id", "method", "species", "det", "days") %in% names(dat0)))

# Focus paper: on garde uniquement les méthodes demandées dans methods_use.
dat_eff <- dat0 %>%
  dplyr::filter(method %in% methods_use)

# Univers complet sites x espèces: indispensable pour coder les absences en 0
# lorsque la visite existe, et éviter les faux NA espèce-site.
sites <- sort(unique(dat_eff$site_id))
spp   <- sort(unique(dat_eff$species))
site_species_grid <- tidyr::expand_grid(site_id = sites, species = spp)

# ---------- 2. Détection par sous-ensemble Pitfall ----------
det_tabs <- list()
avail_tabs <- list()
effort_tabs <- list()

if ("Pitfall" %in% methods_use) {
  pit_dat <- dat_eff %>%
    dplyr::filter(method == "Pitfall", !is.na(trap), trap %in% 1:max(pit_reps))

  # Jours par site pour Pitfall
  pit_days_site <- pit_dat %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(days_pit = med_days(days), .groups = "drop")

  # Nombre effectif de pots récupérés dans chaque sous-ensemble nominal 1:k
  pit_avail_long <- purrr::map_dfr(pit_reps, function(k) {
    pit_dat %>%
      dplyr::filter(trap %in% 1:k) %>%
      dplyr::group_by(site_id) %>%
      dplyr::summarise(n_pit_eff = dplyr::n_distinct(trap), .groups = "drop") %>%
      dplyr::right_join(tibble::tibble(site_id = sites), by = "site_id") %>%
      dplyr::mutate(
        n_pit_eff = tidyr::replace_na(n_pit_eff, 0L),
        nominal_k = k,
        det_col = paste0("Pitfall", k),
        visit_available = if (isTRUE(require_complete_pitfall_subsets)) {
          n_pit_eff >= nominal_k
        } else {
          n_pit_eff > 0
        }
      ) %>%
      dplyr::left_join(pit_days_site, by = "site_id") %>%
      dplyr::mutate(
        days_pit = tidyr::replace_na(days_pit, 1),
        effort = dplyr::if_else(visit_available, n_pit_eff * pit_circ_m * days_pit, NA_real_)
      )
  })

  # Table large des efforts Pitfall: eff_p10, eff_p8, ...
  eff_pit_site <- pit_avail_long %>%
    dplyr::transmute(site_id, eff_name = paste0("eff_p", nominal_k), effort) %>%
    tidyr::pivot_wider(names_from = eff_name, values_from = effort)

  # Table large des disponibilités par colonne Pitfall
  avail_pit_site <- pit_avail_long %>%
    dplyr::transmute(site_id, det_col, visit_available = as.integer(visit_available)) %>%
    tidyr::pivot_wider(names_from = det_col, values_from = visit_available, values_fill = 0L)

  effort_tabs$pit <- eff_pit_site
  avail_tabs$pit <- avail_pit_site

  # Détection espèce x site pour chaque sous-ensemble nominal
  for (k in pit_reps) {
    nm <- paste0("Pitfall", k)
    det_k <- pit_dat %>%
      dplyr::filter(trap %in% 1:k) %>%
      dplyr::group_by(site_id, species) %>%
      dplyr::summarise(!!nm := as.integer(any(det == 1, na.rm = TRUE)), .groups = "drop")
    det_tabs[[nm]] <- det_k
  }

  readr::write_csv(
    pit_avail_long %>% dplyr::arrange(site_id, dplyr::desc(nominal_k)),
    file.path(reviewer_out_dir, "pitfall_effective_traps_by_site.csv")
  )
}

# ---------- 3. Détection GPD et effort effectif ----------
if ("GPD" %in% methods_use) {
  gpd_dat <- dat_eff %>% dplyr::filter(method == "GPD")

  gpd_days_site <- gpd_dat %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(days_gpd = med_days(days), .groups = "drop")

  # IMPORTANT: n_gpd_eff dépend ici des sous-pièges présents dans dat0.
  # Si les pièges vides ne sont pas représentés dans le fichier taxonomique,
  # fournir idéalement une table d'effort externe. Sinon, cette estimation est conservatrice.
  gpd_avail_site <- gpd_dat %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(
      n_gpd_eff = dplyr::n_distinct(gpd_sub[!is.na(gpd_sub)]),
      .groups = "drop"
    ) %>%
    dplyr::right_join(tibble::tibble(site_id = sites), by = "site_id") %>%
    dplyr::mutate(n_gpd_eff = tidyr::replace_na(n_gpd_eff, 0L)) %>%
    dplyr::left_join(gpd_days_site, by = "site_id") %>%
    dplyr::mutate(
      days_gpd = tidyr::replace_na(days_gpd, 1),
      visit_available = n_gpd_eff > 0,
      IL_effective_m = dplyr::if_else(visit_available,
                                      gpd_fence_length_m + n_gpd_eff * pit_circ_m,
                                      NA_real_),
      eff_gpd = IL_effective_m * days_gpd,
      GPD = as.integer(visit_available)
    )

  effort_tabs$gpd <- gpd_avail_site %>% dplyr::select(site_id, eff_gpd, n_gpd_eff, IL_effective_m)
  avail_tabs$gpd <- gpd_avail_site %>% dplyr::select(site_id, GPD)

  det_tabs$GPD <- gpd_dat %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(GPD = as.integer(any(det == 1, na.rm = TRUE)), .groups = "drop")

  readr::write_csv(
    gpd_avail_site %>% dplyr::arrange(site_id),
    file.path(reviewer_out_dir, "gpd_effective_pots_by_site.csv")
  )
}

# ---------- 4. Optionnel: DVAC / TRI conservés si présents ----------
if ("DVAC" %in% methods_use) {
  det_tabs$DVAC <- dat_eff %>%
    dplyr::filter(method == "DVAC") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(DVAC = as.integer(any(det == 1, na.rm = TRUE)), .groups = "drop")
  effort_tabs$dvac <- dat_eff %>% dplyr::filter(method == "DVAC") %>%
    dplyr::distinct(site_id) %>% dplyr::mutate(eff_dvac = 1)
  avail_tabs$dvac <- dat_eff %>% dplyr::filter(method == "DVAC") %>%
    dplyr::distinct(site_id) %>% dplyr::mutate(DVAC = 1L)
}

if ("TRI" %in% methods_use) {
  det_tabs$TM <- dat_eff %>%
    dplyr::filter(method == "TRI") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(TM = as.integer(any(det == 1, na.rm = TRUE)), .groups = "drop")
  effort_tabs$tm <- dat_eff %>% dplyr::filter(method == "TRI") %>%
    dplyr::distinct(site_id) %>% dplyr::mutate(eff_tm = 6 * tm_unit)
  avail_tabs$tm <- dat_eff %>% dplyr::filter(method == "TRI") %>%
    dplyr::distinct(site_id) %>% dplyr::mutate(TM = 1L)
}

# ---------- 5. Tables larges: détection, disponibilité, effort ----------
# Détection: full grid + left joins; NA -> 0 seulement ensuite si la visite existe.
det_wide0 <- site_species_grid
if (length(det_tabs)) {
  for (nm in names(det_tabs)) {
    det_wide0 <- dplyr::left_join(det_wide0, det_tabs[[nm]], by = c("site_id", "species"))
  }
}

det_cols <- setdiff(names(det_wide0), c("site_id", "species"))

# Disponibilité site x colonne de détection
sampling_map <- tibble::tibble(site_id = sites)
if (length(avail_tabs)) {
  for (nm in names(avail_tabs)) {
    sampling_map <- dplyr::left_join(sampling_map, avail_tabs[[nm]], by = "site_id")
  }
}
for (cl in setdiff(names(sampling_map), "site_id")) {
  sampling_map[[cl]] <- tidyr::replace_na(as.integer(sampling_map[[cl]]), 0L)
}

# Effort site x variable effort
effort_wide <- tibble::tibble(site_id = sites)
if (length(effort_tabs)) {
  for (nm in names(effort_tabs)) {
    effort_wide <- dplyr::left_join(effort_wide, effort_tabs[[nm]], by = "site_id")
  }
}

# Application de la règle: si visite disponible -> NA detection devient 0; sinon NA.
det_wide <- det_wide0
for (cl in det_cols) {
  if (!cl %in% names(sampling_map)) next
  idx <- match(det_wide$site_id, sampling_map$site_id)
  available <- sampling_map[[cl]][idx] == 1L
  x <- det_wide[[cl]]
  x[available & is.na(x)] <- 0L
  x[!available] <- NA_integer_
  det_wide[[cl]] <- as.integer(x)
}

det_wide <- det_wide %>% dplyr::select(site_id, species, dplyr::all_of(det_cols))

# ---------- 6. Matrice Y ----------
Rnames <- det_cols
# ordre plus lisible si Pitfall + GPD
preferred_order <- c(paste0("Pitfall", pit_reps), "GPD", "DVAC", "TM")
Rnames <- intersect(preferred_order, Rnames)

Y <- array(NA_integer_,
           dim = c(length(sites), length(Rnames), length(spp)),
           dimnames = list(site = sites, rep = Rnames, species = spp))

det_long <- det_wide %>%
  dplyr::mutate(site_id = factor(site_id, levels = sites),
                species = factor(species, levels = spp)) %>%
  dplyr::arrange(species, site_id)

for (k in seq_along(spp)) {
  dk <- det_long %>%
    dplyr::filter(species == spp[k]) %>%
    dplyr::arrange(site_id)
  # avec full grid, nrow(dk) doit être length(sites)
  if (nrow(dk) != length(sites)) {
    dk <- tibble::tibble(site_id = factor(sites, levels = sites)) %>%
      dplyr::left_join(dk, by = "site_id")
  }
  for (r in seq_along(Rnames)) {
    Y[, r, k] <- as.integer(dk[[Rnames[r]]])
  }
}

# ---------- 7. Matrices effort brut et standardisé ----------
eff_mat <- matrix(NA_real_, nrow = length(sites), ncol = length(Rnames),
                  dimnames = list(sites, Rnames))

map_eff_name <- c(
  stats::setNames(paste0("eff_p", pit_reps), paste0("Pitfall", pit_reps)),
  "GPD" = "eff_gpd",
  "DVAC" = "eff_dvac",
  "TM" = "eff_tm"
)

for (rn in Rnames) {
  nm <- map_eff_name[[rn]]
  if (!is.null(nm) && nm %in% names(effort_wide)) {
    eff_mat[, rn] <- effort_wide[[nm]][match(sites, effort_wide$site_id)]
  }
}

# Masque des visites disponibles indépendant de l'espèce: au moins une valeur non-NA
# dans Y pour la première espèce suffit car disponibilité est site x méthode.
visit_mask <- matrix(FALSE, nrow = length(sites), ncol = length(Rnames),
                     dimnames = list(sites, Rnames))
for (rn in Rnames) {
  if (rn %in% names(sampling_map)) {
    visit_mask[, rn] <- sampling_map[[rn]][match(sites, sampling_map$site_id)] == 1L
  }
}
eff_mat[!visit_mask] <- NA_real_

eff_log <- log1p(as.vector(eff_mat))
eff_mean <- mean(eff_log, na.rm = TRUE)
eff_sd <- stats::sd(eff_log, na.rm = TRUE)
if (!is.finite(eff_sd) || is.na(eff_sd) || eff_sd == 0) eff_sd <- 1

eff_z <- matrix((eff_log - eff_mean) / eff_sd,
                nrow = length(sites), ncol = length(Rnames), byrow = FALSE,
                dimnames = dimnames(eff_mat))
eff_z[!visit_mask] <- NA_real_

# ---------- 8. Tables reviewer-ready ----------
effort_definition_table <- tibble::tibble(
  method_level = Rnames,
  nominal_unit = dplyr::case_when(
    stringr::str_detect(method_level, "^Pitfall") ~ "nested nominal pitfall subset",
    method_level == "GPD" ~ "one X-shaped directional pitfall device",
    TRUE ~ "other"
  ),
  effort_formula = dplyr::case_when(
    stringr::str_detect(method_level, "^Pitfall") ~ "n_effective_pots_in_subset × pitfall opening circumference × trapping duration",
    method_level == "GPD" ~ "(drift-fence length + n_effective_pots × pitfall opening circumference) × trapping duration",
    method_level == "DVAC" ~ "constant nominal effort",
    method_level == "TM" ~ "6 monoliths × 0.0625 m²",
    TRUE ~ NA_character_
  ),
  pitfall_opening_circumference_m = pit_circ_m,
  gpd_fence_length_m = ifelse(method_level == "GPD", gpd_fence_length_m, NA_real_)
)
readr::write_csv(effort_definition_table, file.path(reviewer_out_dir, "effort_definition_table.csv"))

visit_summary <- sampling_map %>%
  tidyr::pivot_longer(-site_id, names_to = "method_level", values_to = "available") %>%
  dplyr::group_by(method_level) %>%
  dplyr::summarise(
    n_sites_available = sum(available == 1L, na.rm = TRUE),
    n_sites_total = dplyr::n(),
    prop_sites_available = n_sites_available / n_sites_total,
    .groups = "drop"
  )
readr::write_csv(visit_summary, file.path(reviewer_out_dir, "visit_availability_by_method.csv"))

# ---------- 9. Covariables site ----------
site_alt <- dat_eff %>%
  dplyr::group_by(site_id) %>%
  dplyr::summarise(ALTITUDE = stats::median(ALTITUDE, na.rm = TRUE), .groups = "drop")

site_doy <- dat_eff %>%
  dplyr::group_by(site_id) %>%
  dplyr::summarise(DOY = stats::median(as.integer(format(DATE, "%j")), na.rm = TRUE), .groups = "drop")

# NDVI optionnel: ne bloque pas le pipeline si le fichier manque.
ndvi_path <- "data/derived-data/ndvi.csv"
if (file.exists(ndvi_path)) {
  site_ndvi <- read.csv(ndvi_path, h = TRUE, sep = ";") %>%
    dplyr::rename(site_id = site) %>%
    dplyr::mutate(site_id = as.character(site_id))
} else {
  site_ndvi <- tibble::tibble(site_id = sites, ndvi_m = NA_real_, ndvi_sd = NA_real_)
}

site_cov <- tibble::tibble(site_id = sites) %>%
  dplyr::left_join(site_alt, by = "site_id") %>%
  dplyr::left_join(site_doy, by = "site_id") %>%
  dplyr::left_join(site_ndvi, by = "site_id") %>%
  dplyr::mutate(
    ALTITUDE_z = zstd(ALTITUDE),
    DOY_z      = zstd(DOY),
    NDVI_m_z   = zstd(ndvi_m),
    NDVI_sd_z  = zstd(ndvi_sd)
  ) %>%
  dplyr::select(site_id, ALTITUDE_z, DOY_z, NDVI_m_z, NDVI_sd_z)

# Habitat si disponible dans dat0
if ("HABITAT_OC" %in% names(dat_eff)) {
  site_hab <- dat_eff %>%
    dplyr::filter(!is.na(HABITAT_OC)) %>%
    dplyr::count(site_id, HABITAT_OC, name = "n") %>%
    dplyr::group_by(site_id) %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(site_id, HABITAT_OC)
  site_cov <- site_cov %>% dplyr::left_join(site_hab, by = "site_id")
}

# ---------- 10. Contrôles console ----------
message("04_detection_effort_reviewed_v2: ", length(sites), " sites, ", length(spp), " species, reps = ", paste(Rnames, collapse = ", "))
message("Visit availability:")
print(visit_summary)
message("NA in Y: ", sum(is.na(Y)), " / ", length(Y))
message("NA in eff_z: ", sum(is.na(eff_z)), " / ", length(eff_z))
