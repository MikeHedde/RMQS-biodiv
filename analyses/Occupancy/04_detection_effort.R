# ============================================================
# DÉTECTIONS PAR RÉPLICAT & MATRICES D’EFFORT
# Sorties: det_wide, det_cols, Y, eff_z, site_cov, sites, spp, Rnames
# ============================================================

# 4.1 Détections
det_tabs <- list()

if ("Pitfall" %in% methods_use) {
  make_det_pit <- function(n) {
    nm <- paste0("Pitfall", n)
    dat0 %>%
      filter(method == "Pitfall", !is.na(trap), trap %in% 1:n) %>%
      group_by(site_id, species) %>%
      summarise(!!nm := as.integer(any(det == 1, na.rm = TRUE)), .groups = "drop")
  }
  det_tabs <- c(det_tabs, setNames(lapply(pit_reps, make_det_pit), paste0("Pitfall", pit_reps)))
}
if ("GPD" %in% methods_use) {
  det_tabs$GPD <- dat0 %>% filter(method == "GPD") %>%
    group_by(site_id, species) %>% summarise(GPD = as.integer(any(det==1, na.rm=TRUE)), .groups="drop")
}
if ("DVAC" %in% methods_use) {
  det_tabs$DVAC <- dat0 %>% filter(method == "DVAC") %>%
    group_by(site_id, species) %>% summarise(DVAC = as.integer(any(det==1, na.rm=TRUE)), .groups="drop")
}
if ("TRI" %in% methods_use) {
  det_tabs$TM <- dat0 %>% filter(method == "TRI") %>%
    group_by(site_id, species) %>% summarise(TM = as.integer(any(det==1, na.rm=TRUE)), .groups="drop")
}

keys <- dat0 %>% distinct(site_id, species)
det_wide0 <- keys
if (length(det_tabs)) det_wide0 <- reduce(det_tabs, ~left_join(.x, .y, by=c("site_id","species")))
det_cols <- setdiff(names(det_wide0), c("site_id","species"))

# 4.2 Présence des méthodes par site
effort_parts <- list(site_id = tibble(site_id = sort(unique(dat0$site_id))))

if ("Pitfall" %in% methods_use) {
  eff_pit_site <- dat0 %>%
    filter(method == "Pitfall", !is.na(trap), trap %in% 1:10) %>%
    group_by(site_id) %>%
    summarise(
      days_pit = if (all(is.na(days))) 1 else stats::median(days, na.rm = TRUE),
      n_trap_1_10 = n_distinct(trap[trap %in% 1:10]),
      n_trap_1_8  = n_distinct(trap[trap %in% 1:8]),
      n_trap_1_6  = n_distinct(trap[trap %in% 1:6]),
      n_trap_1_4  = n_distinct(trap[trap %in% 1:4]),
      n_trap_1_2  = n_distinct(trap[trap %in% 1:2]),
      .groups = "drop"
    ) %>%
    mutate(
      eff_p10 = pit_circ_m * days_pit * n_trap_1_10,
      eff_p8  = pit_circ_m * days_pit * n_trap_1_8,
      eff_p6  = pit_circ_m * days_pit * n_trap_1_6,
      eff_p4  = pit_circ_m * days_pit * n_trap_1_4,
      eff_p2  = pit_circ_m * days_pit * n_trap_1_2
    ) %>% select(site_id, starts_with("eff_p"))
  effort_parts$pit <- eff_pit_site
}

if ("GPD" %in% methods_use) {
  eff_gpd_site <- dat0 %>% filter(method == "GPD") %>%
    group_by(site_id) %>%
    summarise(
      days_gpd = if (all(is.na(days))) 1 else stats::median(days, na.rm = TRUE),
      n_active_subtraps = n_distinct(gpd_sub[!is.na(gpd_sub)]),
      .groups = "drop"
    ) %>%
    mutate(
      n_active_subtraps = ifelse(is.na(n_active_subtraps), 0L, n_active_subtraps),
      eff_gpd = gpd_unit_intensity * days_gpd * pmax(0, n_active_subtraps)
    ) %>% select(site_id, eff_gpd, n_active_subtraps)
  effort_parts$gpd <- eff_gpd_site
}

if ("DVAC" %in% methods_use) {
  effort_parts$dvac <- dat0 %>% filter(method == "DVAC") %>% distinct(site_id) %>% mutate(eff_dvac = 1)
}
if ("TRI" %in% methods_use) {
  effort_parts$tm <- dat0 %>% filter(method == "TRI") %>% distinct(site_id) %>% mutate(eff_tm = 6 * tm_unit)
}

effort_wide <- reduce(effort_parts, ~left_join(.x, .y, by = "site_id")) %>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))

sampling_map <- dat0 %>%
  distinct(site_id, method) %>%
  tidyr::pivot_wider(names_from  = method, values_from = method, values_fill = NA, names_expand = FALSE) %>%
  mutate(across(-site_id, ~ as.integer(!is.na(.))))

if ("GPD" %in% methods_use) {
  sampling_map <- sampling_map %>%
    left_join(effort_wide %>% select(site_id, eff_gpd, n_active_subtraps), by="site_id") %>%
    mutate(GPD = ifelse(eff_gpd > 0 & n_active_subtraps > 0, 1L, 0L)) %>%
    select(-eff_gpd, -n_active_subtraps)
}

det_wide <- det_wide0 %>% mutate(across(all_of(det_cols), ~ tidyr::replace_na(., 0L)))

for (cl in det_cols) {
  if (!cl %in% names(det_wide)) next
  method_name <- sub("^Pitfall\\d+$", "Pitfall", cl)
  if (!method_name %in% names(sampling_map)) next
  idx <- match(det_wide$site_id, sampling_map$site_id)
  present_vec <- sampling_map[[method_name]][idx]; present_vec[is.na(present_vec)] <- 0L
  x <- as.integer(det_wide[[cl]]); x[present_vec != 1L] <- NA_integer_
  det_wide[[cl]] <- x
}

det_wide <- det_wide %>% select(site_id, species, all_of(det_cols))

# 4.3 Y & Effort standardisé
sites <- sort(unique(det_wide$site_id))
spp   <- sort(unique(det_wide$species))
Rnames <- det_cols

Y <- array(NA_integer_, dim = c(length(sites), length(Rnames), length(spp)),
           dimnames = list(site = sites, rep = Rnames, species = spp))

det_long <- det_wide %>%
  mutate(site_id = factor(site_id, levels = sites),
         species = factor(species, levels = spp)) %>%
  arrange(species, site_id)

for (k in seq_along(spp)) {
  dk <- det_long %>% filter(species == spp[k]) %>% arrange(site_id)
  if (nrow(dk) < length(sites)) dk <- tibble(site_id = factor(sites, levels = sites)) %>% left_join(dk, by="site_id")
  for (r in seq_along(Rnames)) Y[, r, k] <- as.numeric(dk[[Rnames[r]]])
}

eff_mat <- matrix(0, nrow = length(sites), ncol = length(Rnames), dimnames = list(sites, Rnames))
map_eff_name <- c(setNames(paste0("eff_p", pit_reps), paste0("Pitfall", pit_reps)),
                  "GPD" = "eff_gpd", "DVAC" = "eff_dvac", "TM" = "eff_tm")
for (rn in Rnames) {
  nm <- map_eff_name[rn]
  if (!is.na(nm) && nm %in% names(effort_wide)) {
    eff_mat[, rn] <- effort_wide[[nm]][match(sites, effort_wide$site_id)]
  }
}
eff_vec <- as.vector(eff_mat); eff_log <- log1p(eff_vec)
eff_mean <- mean(eff_log, na.rm = TRUE); eff_sd <- sd(eff_log, na.rm = TRUE); if (!is.finite(eff_sd) || eff_sd==0) eff_sd <- 1
eff_z <- matrix((eff_log - eff_mean)/eff_sd, nrow = length(sites), ncol = length(Rnames), byrow = FALSE, dimnames = dimnames(eff_mat))
mask_allNA_sr <- apply(Y, c(1,2), function(a) all(is.na(a)))
eff_z[mask_allNA_sr] <- NA

# 4.4 Covariables site
site_alt <- dat0 %>% group_by(site_id) %>% summarise(ALTITUDE = median(ALTITUDE, na.rm=TRUE), .groups="drop")
site_doy <- dat0 %>% group_by(site_id) %>% summarise(DOY = median(as.integer(format(DATE, "%j")),  na.rm=TRUE), .groups="drop")
site_ndvi <- read.csv("data/derived-data/ndvi.csv", h = T, sep = ";") %>%
  rename(site_id = site)%>%
  mutate(site_id = as.character(site_id))

site_cov <- tibble(site_id = sites) %>%
  left_join(site_alt, by="site_id") %>%
  left_join(site_doy, by="site_id") %>%
  rename(site_id_full = site_id) %>%
  separate(site_id_full,
           into = c("programme", "year", "locality", "site_id"),
           sep = "_",
           remove = FALSE) %>%
  left_join(site_hab, by="site_id") %>%
  left_join(site_ndvi) %>% 
  mutate(ALTITUDE_z = zstd(ALTITUDE), 
         DOY_z = zstd(DOY),
         NDVI_m_z = zstd(ndvi_m),
         NDVI_sd_z = zstd(ndvi_sd)
         ) %>%
  select(site_id_full, ALTITUDE_z, DOY_z, NDVI_m_z, NDVI_sd_z) %>%
  rename(site_id = site_id_full)
