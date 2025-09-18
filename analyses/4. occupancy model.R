# ============================================================
# Occupancy espèce-par-espèce (unmarked::occu) avec effort
# Réplicats = Pitfall10 / Pitfall6 / Pitfall3 / GPD (GPD = 1 colonne)
# Effort GPD ajuste le nb de sous-pièges actifs (1..4)
# ψ ~ ALTITUDE_z + DOY_z ; p ~ méthode + effort_z
# ============================================================

## 0) Paramètres ---------------------------------------------------------------
in_path   <- "data/raw-data/lab files/databases update/update_araneae.csv"
out_dir   <- "data/derived-data/occupancy"
target_yr <- "2025"
target_project <- "RMQS_2025"
min_sites <- 2
taxo_level <- "ES"
require_contrast <- TRUE

# Constantes effort
pit_diam_m         <- 0.05
pit_circ_m         <- pi * pit_diam_m
gpd_open_circ_m    <- pi * pit_diam_m
gpd_cross_m        <- 1.0
gpd_unit_intensity <- gpd_open_circ_m + gpd_cross_m
tm_unit            <- 0.0625

## 1) Packages ----------------------------------------------------------------
pkgs <- c()
librarian::shelf(dplyr, tidyr, ggplot2, vegan, forcats, divent, ggdist,
                 readr,stringr,tibble,purrr,unmarked, lme4,
                 pbkrtest, lmerTest)

## 2) Lecture CSV et détection colonnes ---------------------------------------
raw <- readr::read_delim(in_path, delim = ";", locale = readr::locale(encoding = "Latin1"),
                         show_col_types = FALSE)

col_exists <- function(x) x[x %in% names(raw)][1]
col_station  <- col_exists(c("NOM_STATION"))
col_method   <- col_exists(c("METHODE"))
col_species  <- col_exists(c("LB_NOM"))
col_abund    <- col_exists(c("ABONDANCE_TOTALE"))
col_sample   <- col_exists(c("ID_ECHANTILLON"))
col_project  <- col_exists(c("PROJET"))
col_DATE     <- col_exists(c("DATE"))
col_DATEFIN  <- col_exists(c("DATE_FIN"))
col_alt      <- col_exists(c("ALTITUDE"))
col_rang     <- col_exists(c("RANG"))

stopifnot(!is.na(col_station), !is.na(col_method), !is.na(col_species),
          !is.na(col_abund), !is.na(col_sample), !is.na(col_rang))

## 3) Normalisation / parsing --------------------------------------------------
norm_method <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x, "(?i)gpd|directionnel|croisillon") ~ "GPD",
    stringr::str_detect(x, "(?i)barb|pitfall|piege\\s*barber") ~ "Pitfall",
    stringr::str_detect(x, "(?i)dvac|aspir") ~ "DVAC",          # TODO futur
    stringr::str_detect(x, "(?i)tri\\s*manuel|tri\\s*manu") ~ "TRI",# TODO futur
    TRUE ~ x
  )
}
parse_abund <- function(x) suppressWarnings(as.numeric(x))
parse_date_fr <- function(x) as.Date(suppressWarnings(readr::parse_date(as.character(x), "%d/%m/%Y")))

dat0 <- raw %>%
  filter(PROJET %in% c(target_project),
         RANG %in% c(taxo_level)) %>%
  transmute(
    site_id = as.character(.data[[col_station]]),
    method  = norm_method(.data[[col_method]]),
    species = as.character(.data[[col_species]]),
    abund   = parse_abund(.data[[col_abund]]),
    sample  = as.character(.data[[col_sample]]),
    year    = if (!is.na(col_project))
      stringr::str_extract(as.character(.data[[col_project]]), "(?<=RMQS_)\\d{4}")
    else NA_character_,
    DATE     = if (!is.na(col_DATE))    parse_date_fr(.data[[col_DATE]])    else as.Date(NA),
    DATE_FIN = if (!is.na(col_DATEFIN)) parse_date_fr(.data[[col_DATEFIN]]) else as.Date(NA),
    ALTITUDE = if (!is.na(col_alt)) suppressWarnings(as.numeric(.data[[col_alt]])) else NA_real_
  ) 

stopifnot(nrow(dat0) > 0)

# Pièges individuels Pitfall (PB1..PB10)
dat0 <- dat0 %>%
  mutate(
    trap = ifelse(method == "Pitfall",
                  suppressWarnings(readr::parse_integer(stringr::str_extract(sample, "(?<=PB)\\d+$"))),
                  NA_integer_),
    det  = as.integer(abund > 0),
    days = {
      d <- as.numeric(DATE_FIN - DATE) + 1
      ifelse(is.na(d) | d < 1, 1, d)
    }
  )

message("Comptes par méthode (", target_yr, ") :"); print(table(dat0$method, useNA = "ifany"))

## 4) Détections agrégées par réplicat -----------------------------------------
det_pit10 <- dat0 %>%
  filter(method == "Pitfall", !is.na(trap), trap %in% 1:10) %>%
  group_by(site_id, species) %>%
  summarise(Pitfall10 = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

det_pit8 <- dat0 %>%
  filter(method == "Pitfall", !is.na(trap), trap %in% 1:8) %>%
  group_by(site_id, species) %>%
  summarise(Pitfall8 = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

det_pit6 <- dat0 %>%
  filter(method == "Pitfall", !is.na(trap), trap %in% 1:6) %>%
  group_by(site_id, species) %>%
  summarise(Pitfall6 = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

det_pit4 <- dat0 %>%
  filter(method == "Pitfall", !is.na(trap), trap %in% 1:4) %>%
  group_by(site_id, species) %>%
  summarise(Pitfall4 = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

det_pit2 <- dat0 %>%
  filter(method == "Pitfall", !is.na(trap), trap %in% 1:2) %>%
  group_by(site_id, species) %>%
  summarise(Pitfall2 = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

# GPD : 1 réplicat (détection si ≥1 sous-piège a détecté)
det_gpd <- dat0 %>%
  filter(method == "GPD") %>%
  group_by(site_id, species) %>%
  summarise(GPD = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

# DVAC : 1 réplicat 
det_dvac <- dat0 %>%
  filter(method == "DVAC") %>%
  group_by(site_id, species) %>%
  summarise(DVAC = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

# DVAC : 1 réplicat 
det_tm <- dat0 %>%
  filter(method == "Tri manuel") %>%
  group_by(site_id, species) %>%
  summarise(TM = as.integer(any(det == 1, na.rm=TRUE)), .groups="drop")

## 5) Effort par site × méthode -----------------------------------------------

# --- Pitfall 
eff_pit_site <- dat0 %>%
  filter(method == "Pitfall", !is.na(trap), trap %in% 1:10) %>%
  group_by(site_id) %>%
  summarise(
    days_pit    = if (all(is.na(days))) 1 else stats::median(days, na.rm = TRUE),
    n_trap_1_10 = n_distinct(trap[trap %in% 1:10]),
    n_trap_1_8 = n_distinct(trap[trap %in% 1:8]),
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
  )

# --- GPD : nb de sous-pièges ACTIFS (ayant au moins une identification) ------
# On récupère l’ID de sous-piège en lisant le dernier nombre de ID_ECHANTILLON.
extract_last_int <- function(x) suppressWarnings(readr::parse_integer(stringr::str_extract(x, "(\\d+)\\s*$")))
gpd_subtrap <- dat0 %>%
  filter(method == "GPD") %>%
  mutate(sub_id = extract_last_int(sample))

eff_gpd_site <- gpd_subtrap %>%
  group_by(site_id) %>%
  summarise(
    days_gpd = if (all(is.na(days))) 1 else stats::median(days, na.rm = TRUE),
    n_active_subtraps = n_distinct(sub_id[!is.na(sub_id)]),  # 0..4
    .groups = "drop"
  ) %>%
  mutate(
    n_active_subtraps = ifelse(is.na(n_active_subtraps), 0L, n_active_subtraps),
    eff_gpd = gpd_unit_intensity * days_gpd * pmax(0, n_active_subtraps)
  )

# DVAC : 100 aspirations=1 sample. 1 seul sample par site
eff_dvac_site <- dat0 %>%
  filter(method == "DVAC" )%>%
  select(site_id) %>%
  mutate(eff_dvac = 1)

# --- TM : nb de repetitions = 6 ------
eff_tm_site <- dat0 %>%
  filter(method == "TRI") %>%
  select(site_id) %>%
  mutate(eff_tm = 6 * tm_unit) %>%
  unique()

#--- Fusion effort (une ligne par site, colonnes par réplicat)
sites_all <- sort(unique(dat0$site_id))
effort_wide <- tibble(site_id = sites_all) %>%
  left_join(eff_pit_site %>% select(site_id, eff_p10, eff_p8, eff_p6, eff_p4, eff_p2), by="site_id") %>%
  left_join(eff_gpd_site %>% select(site_id, eff_gpd, n_active_subtraps), by="site_id") %>%
  left_join(eff_dvac_site %>% select(site_id, eff_dvac), by="site_id") %>%
  left_join(eff_tm_site %>% select(site_id, eff_tm), by="site_id") %>%
  mutate(across(c(eff_p10, eff_p8, eff_p6, eff_p4, eff_p2, eff_gpd, eff_dvac, eff_tm), ~ tidyr::replace_na(., 0)),  
         n_active_subtraps = tidyr::replace_na(n_active_subtraps, 0L))

effort_long <- effort_wide %>%
  pivot_longer(cols = 2:ncol(.), names_to = "method", values_to = "values") %>%
  filter(method %in% c("eff_p10", "eff_p8", "eff_p6", "eff_p4","eff_p2", "eff_gpd", "eff_dvac", "eff_tm")) %>%  
  filter(values>0)

effort_long$method <- factor(effort_long$method, levels = c("eff_gpd", "eff_p10", "eff_p8", "eff_p6", "eff_p4","eff_p2", "eff_dvac"))

#--- Visualisation des efforts
ggplot(effort_long, aes(x=method, y =values))+
  geom_boxplot()+
  geom_jitter()+
  labs(y="Sampling effort (meter.days)")+
  theme_bw()

## 6) Table des méthodes disponibles (GPD absent si 0 sous-piège actif) --------
sampling_map <- dat0 %>%
  distinct(site_id, method) %>%
  tidyr::pivot_wider(names_from = method, values_from = method, values_fill = NA) %>%
  transmute(
    site_id,
    has_pit_raw = as.integer(!is.na(Pitfall)),
    has_gpd_raw = as.integer(!is.na(GPD)),
    has_dvac_raw = as.integer(!is.na(DVAC)),
    has_tm_raw = as.integer(!is.na(TRI))
  ) %>%
  left_join(effort_wide %>% select(site_id, eff_gpd, n_active_subtraps), by="site_id") %>%
  mutate(
    has_pit = has_pit_raw,
    has_gpd = ifelse(eff_gpd > 0 & n_active_subtraps > 0, 1L, 0L),
    has_dvac = has_dvac_raw,
    has_tm = has_tm_raw
  ) %>%
  select(site_id, has_pit, has_gpd, has_dvac, has_tm) %>%
  unique()

## 7) Détections + NA si méthode absente ---------------------------------------
keys <- dat0 %>% distinct(site_id, species)
det_wide <- keys %>%
  left_join(det_pit10, by = c("site_id","species")) %>%
  left_join(det_pit8,  by = c("site_id","species")) %>%
  left_join(det_pit6,  by = c("site_id","species")) %>%
  left_join(det_pit4,  by = c("site_id","species")) %>%
  left_join(det_pit2,  by = c("site_id","species")) %>%
  left_join(det_gpd,   by = c("site_id","species")) %>%
  left_join(det_tm,    by = c("site_id","species")) %>%
  left_join(det_dvac,  by = c("site_id","species")) %>%
  mutate(across(c(Pitfall10,Pitfall8,Pitfall6,Pitfall4,Pitfall2, GPD, DVAC, TM), ~ tidyr::replace_na(., 0L))) %>%
  left_join(sampling_map, by="site_id") %>%
  mutate(
    Pitfall10  = ifelse(has_pit==1L, Pitfall10, NA),
    Pitfall8   = ifelse(has_pit==1L, Pitfall8,  NA),
    Pitfall6   = ifelse(has_pit==1L, Pitfall6,  NA),
    Pitfall4   = ifelse(has_pit==1L, Pitfall4,  NA),
    Pitfall2   = ifelse(has_pit==1L, Pitfall2,  NA),
    GPD        = ifelse(has_gpd==1L, GPD,       NA),
    TM        = ifelse(has_tm==1L,  TM,       NA),
    DVAC       = ifelse(has_dvac==1L,DVAC,      NA)
  ) %>%
  select(site_id, species, Pitfall10, Pitfall8,Pitfall6,Pitfall4,Pitfall2, GPD, DVAC, TM)

## 8) Construire Y (sites × réplicats × espèces) + matrices d’effort ------------
sites <- sort(unique(det_wide$site_id))
spp   <- sort(unique(det_wide$species))
J <- length(sites); R <- 8; N <- length(spp)

Y <- array(NA_integer_, dim = c(J, R, N),
           dimnames = list(site = sites,
                           rep  = c("Pitfall10","Pitfall8","Pitfall6","Pitfall4","Pitfall2","GPD", "DVAC", "TM"),
                           species = spp))

det_long <- det_wide %>%
  mutate(site_id = factor(site_id, levels = sites),
         species = factor(species, levels = spp)) %>%
  arrange(species, site_id)

for (k in seq_along(spp)) {
  dk <- det_long %>% filter(species == spp[k]) %>% arrange(site_id)
  if (nrow(dk) < J) {
    dk <- tibble(site_id = factor(sites, levels = sites)) %>% left_join(dk, by="site_id")
  }
  Y[,1,k] <- as.numeric(dk$Pitfall10)
  Y[,2,k] <- as.numeric(dk$Pitfall8)  
  Y[,3,k] <- as.numeric(dk$Pitfall6)
  Y[,4,k] <- as.numeric(dk$Pitfall4)
  Y[,5,k] <- as.numeric(dk$Pitfall2)
  Y[,6,k] <- as.numeric(dk$GPD)
  Y[,7,k] <- as.numeric(dk$DVAC)
  Y[,8,k] <- as.numeric(dk$TM)
  }

# Effort matrices alignées (J × R)
eff_tab <- effort_wide %>% mutate(site_id = factor(site_id, levels = sites)) %>% arrange(site_id) %>% unique()
eff_mat <- matrix(0, nrow=J, ncol=R, 
                  dimnames = list(sites, c("Pitfall10","Pitfall8","Pitfall6","Pitfall4","Pitfall2","GPD", "DVAC", "TM")))
eff_mat[, "Pitfall10"] <- eff_tab$eff_p10
eff_mat[, "Pitfall8"]  <- eff_tab$eff_p8
eff_mat[, "Pitfall6"]  <- eff_tab$eff_p6
eff_mat[, "Pitfall4"]  <- eff_tab$eff_p4
eff_mat[, "Pitfall2"]  <- eff_tab$eff_p2
eff_mat[, "GPD"]       <- eff_tab$eff_gpd
eff_mat[, "DVAC"]      <- eff_tab$eff_dvac
eff_mat[, "TM"]        <- eff_tab$eff_tm

# Standardisation effort (log1p + z)
eff_vec  <- as.vector(eff_mat)
eff_log  <- log1p(eff_vec)
eff_mean <- mean(eff_log, na.rm = TRUE)
eff_sd   <- stats::sd(eff_log, na.rm = TRUE)
eff_zvec <- (eff_log - eff_mean) / ifelse(eff_sd > 0, eff_sd, 1)
eff_z    <- matrix(eff_zvec, nrow = J, ncol = R, byrow = FALSE,
                   dimnames = dimnames(eff_mat))

# Masque site×réplicat tout-NA (toutes espèces NA => pas de données)
mask_allNA_sr <- apply(Y, c(1,2), function(a) all(is.na(a)))
eff_z[mask_allNA_sr] <- NA

## 9) Covariables de site : ALTITUDE_z + DOY_z ---------------------------------
# ALTITUDE par site (médiane si répété)
site_alt <- dat0 %>% group_by(site_id) %>%
  summarise(ALTITUDE = stats::median(ALTITUDE, na.rm = TRUE), .groups="drop")
# DOY = jour de l’année basé sur DATE (médiane par site)
to_doy <- function(d) ifelse(is.na(d), NA, as.integer(format(d, "%j")))
site_doy <- dat0 %>% group_by(site_id) %>%
  summarise(DOY = stats::median(to_doy(DATE), na.rm = TRUE), .groups="drop")

site_cov <- tibble(site_id = sites) %>%
  left_join(site_alt, by="site_id") %>%
  left_join(site_doy, by="site_id")

# z-score (remplace NA par 0 après centrage)
zstd <- function(v) {
  m <- mean(v, na.rm = TRUE); s <- stats::sd(v, na.rm = TRUE)
  z <- (v - m) / ifelse(s > 0, s, 1)
  z[is.na(z)] <- 0
  z
}
site_cov <- site_cov %>%
  mutate(ALTITUDE_z = zstd(ALTITUDE),
         DOY_z      = zstd(DOY)) %>%
  select(site_id, ALTITUDE_z, DOY_z)

row.names_site <- sites  # pour alignement rapide

## 10) Fonction unmarked::occu (ψ ~ ALT + DOY ; p ~ méthode + effort) ----------
method_eff_medians <- function(eff_z_matrix) {
  apply(eff_z_matrix, 2, function(col) if (all(is.na(col))) 0 else stats::median(col, na.rm=TRUE))
}

fit_one_species_occu <- function(sp_name) {
  k <- which(dimnames(Y)$species == sp_name)
  if (length(k) != 1) return(NULL)
  Yk <- Y[,,k, drop=FALSE][,,1]  # J × R
  
  # Sites informatifs
  keep_site <- apply(Yk, 1, function(a) any(!is.na(a)))
  Yk <- Yk[keep_site, , drop = FALSE]
  if (nrow(Yk) < min_sites) return(NULL)
  
  # Réplicats informatifs
  keep_rep <- apply(Yk, 2, function(a) any(!is.na(a)))
  Yk <- Yk[, keep_rep, drop = FALSE]
  if (ncol(Yk) < 2) return(NULL)
  
  # Contraste (≥1 détection & ≥1 non-détection)
  yvec <- as.vector(Yk[!is.na(Yk)])
  if (require_contrast && !(any(yvec == 1) && any(yvec == 0))) return(NULL)
  
  rep_names <- colnames(Yk)
  I8_mat   <- matrix(as.numeric(rep_names=="Pitfall8"), nrow(Yk), ncol(Yk), byrow=TRUE)
  I6_mat   <- matrix(as.numeric(rep_names=="Pitfall6"), nrow(Yk), ncol(Yk), byrow=TRUE)
  I4_mat   <- matrix(as.numeric(rep_names=="Pitfall4"), nrow(Yk), ncol(Yk), byrow=TRUE)
  I2_mat   <- matrix(as.numeric(rep_names=="Pitfall2"), nrow(Yk), ncol(Yk), byrow=TRUE)
  IGPD_mat <- matrix(as.numeric(rep_names=="GPD"),      nrow(Yk), ncol(Yk), byrow=TRUE)
  IDVAC_mat<- matrix(as.numeric(rep_names=="DVAC"),     nrow(Yk), ncol(Yk), byrow=TRUE)
  ITM_mat  <- matrix(as.numeric(rep_names=="TM"),       nrow(Yk), ncol(Yk), byrow=TRUE)
  
  eff_sp <- eff_z[rownames(Yk), rep_names, drop = FALSE]
  I8_mat[is.na(Yk)] <- NA;I6_mat[is.na(Yk)] <- NA;I4_mat[is.na(Yk)] <- NA; I2_mat[is.na(Yk)] <- NA; 
  IGPD_mat[is.na(Yk)] <- NA; IDVAC_mat[is.na(Yk)] <- NA; ITM_mat[is.na(Yk)]<- NA; eff_sp[is.na(Yk)] <- NA
  
  # Covariables de site alignées + variance locale
  sc <- site_cov %>%
    dplyr::filter(site_id %in% rownames(Yk)) %>%
    dplyr::arrange(match(site_id, rownames(Yk))) %>%
    dplyr::select(ALTITUDE_z, DOY_z)
  stopifnot(nrow(sc) == nrow(Yk))
  
  has_var_alt <- stats::sd(sc$ALTITUDE_z, na.rm = TRUE) > 0
  has_var_doy <- stats::sd(sc$DOY_z,      na.rm = TRUE) > 0
  
  # Formule de l’état (ψ) adaptée
  state_terms <- c()
  if (has_var_alt) state_terms <- c(state_terms, "ALTITUDE_z")
  if (has_var_doy) state_terms <- c(state_terms, "DOY_z")
  state_rhs <- if (length(state_terms)) paste(state_terms, collapse = " + ") else "1"
  form_occu <- as.formula(paste("~ I8 + I6 + I4 + I2 + IGPD + IDVAC + ITM + eff_z ~", state_rhs))
  
  umf <- unmarked::unmarkedFrameOccu(
    y = Yk,
    siteCovs = sc,
    obsCovs  = list(I8 = I8_mat, I6 = I6_mat,
                    I4 = I4_mat, I2 = I2_mat, 
                    IGPD = IGPD_mat, IDVAC = IDVAC_mat,
                    ITM = ITM_mat, eff_z = eff_sp)
  )
  
  fm <- try(unmarked::occu(form_occu, data = umf, 
                           control = list(maxit = 200)), 
            silent = TRUE)
  if (inherits(fm, "try-error")) return(NULL)
  
  # Prédictions p(method) à l'effort médian global
  eff_med_all <- method_eff_medians(eff_z)
  eff_for <- function(m) 
    if (m %in% names(eff_med_all)) eff_med_all[[m]] 
    else 0
  mk_pred <- function(I8v, I6v, I4v, I2v, IGPDv, ITMv, IDVACv, effv) {
    # Pour la prédiction de p, l’occupation est conditionnée par (ALT=0, DOY=0) si présents
    nd <- data.frame(I8 = I8v, I6 = I6v, I4 = I4v, I2 = I2v, IGPD = IGPDv, IDVAC = IDVACv, ITM = ITMv, eff_z = effv)
    if (has_var_alt) nd$ALTITUDE_z <- 0
    if (has_var_doy) nd$DOY_z      <- 0
    as.data.frame(unmarked::predict(fm, type = "det", newdata = nd))
  }
  p_p10 <- mk_pred(0,0,0,0,0,0,0, eff_for("Pitfall10"))
  p_p8  <- mk_pred(1,0,0,0,0,0,0, eff_for("Pitfall8"))
  p_p6  <- mk_pred(0,1,0,0,0,0,0, eff_for("Pitfall6"))
  p_p4  <- mk_pred(0,0,1,0,0,0,0, eff_for("Pitfall4"))
  p_p2  <- mk_pred(0,0,0,1,0,0,0, eff_for("Pitfall2"))
  p_gpd <- mk_pred(0,0,0,0,1,0,0, eff_for("GPD"))
  p_dvac<- mk_pred(0,0,0,0,0,1,0, eff_for("DVAC"))
  p_tm  <- mk_pred(0,0,0,0,0,0,1, eff_for("TM"))
  
  tab_p <- tibble::tibble(
    species  = sp_name,
    method   = c("Pitfall10","Pitfall8","Pitfall6","Pitfall4","Pitfall2","GPD", "DVAC", "TM"),
    p_hat    = c(p_p10$Predicted, p_p8$Predicted, p_p6$Predicted, p_p4$Predicted, p_p2$Predicted, p_gpd$Predicted, p_dvac$Predicted, p_tm$Predicted),
    lcl      = c(p_p10$lower,     p_p8$lower,     p_p6$lower,     p_p4$lower,     p_p2$lower,     p_gpd$lower,     p_dvac$lower,     p_tm$lower),
    ucl      = c(p_p10$upper,     p_p8$upper,     p_p6$upper,     p_p4$upper,     p_p2$upper,     p_gpd$upper,     p_dvac$upper,     p_tm$upper),
    n_sites  = nrow(Yk),
    n_rep    = ncol(Yk)
  )
  
  # Coefs ψ (ALT, DOY)
  co <- summary(fm)  # matrix with Estimate/SE
  co_state <- co$state
  beta_ALT <- if ("ALTITUDE_z" %in% rownames(co_state)) co_state[rownames(co_state)=="ALTITUDE_z", colnames(co_state) == "Estimate"] else NA_real_
  se_ALT   <- if ("ALTITUDE_z" %in% rownames(co_state)) co_state[rownames(co_state)=="ALTITUDE_z", colnames(co_state) == "SE"] else NA_real_
  beta_DOY <- if ("DOY_z" %in% rownames(co_state)) co_state[rownames(co_state)=="DOY_z", colnames(co_state) == "Estimate"] else NA_real_
  se_DOY   <- if ("DOY_z" %in% rownames(co_state)) co_state[rownames(co_state)=="DOY_z", colnames(co_state) == "SE"]       else NA_real_
  
  # ψ̂ par site (pour richesse attendue Σψ̂)
  psi_pred <- as.data.frame(unmarked::predict(fm, type="state")) # ordre = sites de Yk
  psi_hat  <- psi_pred$Predicted
  names(psi_hat) <- rownames(Yk)
  
  list(p_table = tab_p,
       beta = tibble::tibble(species=sp_name,
                             beta_alt=beta_ALT, se_alt=se_ALT,
                             beta_doy=beta_DOY, se_doy=se_DOY),
       psi_hat = tibble::tibble(site_id=rownames(Yk),
                                species=sp_name, psi_hat=psi_hat))
}

## 11) Lancer sur toutes les espèces -------------------------------------------
spp_all <- dimnames(Y)$species
res_list <- purrr::map(spp_all, ~try(fit_one_species_occu(.x), silent=TRUE))

order <- raw %>%
  select(LB_NOM, ORDRE)%>%
  unique() %>%
  rename(species = LB_NOM, order = ORDRE)

family <- raw %>%
  select(LB_NOM, FAMILLE)%>%
  unique() %>%
  rename(species = LB_NOM, family = FAMILLE)

tab_p <- res_list %>%
  purrr::keep(~!inherits(.x,"try-error") && !is.null(.x$p_table)) %>%
  purrr::map("p_table") %>% dplyr::bind_rows() %>%
  left_join(family) %>%
  left_join(order)

tab_beta <- res_list %>%
  purrr::keep(~!inherits(.x,"try-error") && !is.null(.x$beta)) %>%
  purrr::map("beta") %>% dplyr::bind_rows()%>%
  left_join(family) %>%
  left_join(order)

tab_psi <- res_list %>%
  purrr::keep(~!inherits(.x,"try-error") && !is.null(.x$psi_hat)) %>%
  purrr::map("psi_hat") %>% dplyr::bind_rows()

readr::write_csv(tab_p,   file.path(out_dir, "p_hat_by_method_unmarked_full.csv"))
readr::write_csv(tab_beta,file.path(out_dir, "psi_coef_ALT_DOY_by_species.csv"))
readr::write_csv(tab_psi, file.path(out_dir, "psi_hat_by_site_species.csv"))

message("Espèces retenues : ", dplyr::n_distinct(tab_p$species), "/", length(spp_all))

# Richesse attendue par site (corrigée détectabilité) : Σ_s ψ̂_is
psi_rich_site <- tab_psi %>%
  group_by(site_id) %>%
  summarise(S_exp = sum(psi_hat, na.rm=TRUE), .groups="drop")

readr::write_csv(psi_rich_site, file.path(out_dir, "expected_richness_by_site.csv"))


## 12) Synthèse communautaire & figures ----------------------------------------
# Détectabilité par méthode
comm_sum <- tab_p %>%
  group_by(method, order, family) %>%
  summarise(p_med = median(p_hat, na.rm=TRUE),
            p_mean = mean(p_hat, na.rm=TRUE),
            lcl_med = median(lcl, na.rm=TRUE),
            ucl_med = median(ucl, na.rm=TRUE),
            n_sp = n_distinct(species), .groups="drop")
print(comm_sum)
readr::write_csv(comm_sum, file.path(out_dir, "community_detection_unmarked_effort_psi_cov.csv"))

# Distribution des effets ALT & DOY sur ψ (par espèce)
tab_beta <- tab_beta %>%
  distinct(species, order, family, beta_alt, se_alt, beta_doy, se_doy)

readr::write_csv(tab_beta, file.path(out_dir, "psi_coefficients_by_species.csv"))

# Petits graphiques de synthèse
g1 <- ggplot(tab_beta, aes(x = beta_alt)) + geom_density(na.rm=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Effect on occupancy (β_ALTITUDE_z)", y="Density", title="Community distribution of altitude effects on ψ") +
  facet_wrap(order+family~.)+
  theme_minimal()

g2 <- ggplot(tab_beta, aes(x = beta_doy)) + geom_density(na.rm=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Effect on occupancy (β_DOY_z)", y="Density", title="Community distribution of date (DOY) effects on ψ") +
  facet_wrap(order+family~.)+
  theme_minimal()

ggsave(file.path(out_dir, "psi_effects_altitude_density.png"), g1, width=6, height=4, dpi=200)
ggsave(file.path(out_dir, "psi_effects_doy_density.png"), g2, width=6, height=4, dpi=200)

# Violin plot p_hat par méthode
tab_p$method <- factor(tab_p$method, levels = c("GPD", "Pitfall10", "Pitfall8", "Pitfall6", "Pitfall4", "Pitfall2", "DVAC", "TM"))
p_fig <- ggplot(tab_p, aes(x = method, y = p_hat)) +
  geom_violin(fill = "grey92") +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  geom_point(data = comm_sum, aes(y = p_med), color = "red", size = 3) +
  facet_wrap(order+family~.)+
  labs(
    y = "p̂ (détectabilité, occupancy unmarked;\n effort médian par méthode)",
    x = NULL,
    title = "Détectabilité par méthode —\najustée pour l’effort (unmarked::occu)"
  ) +
  theme_minimal()

ggplot(tab_p, aes(x = method, y = p_hat)) +
  geom_violin(fill = "grey92") +
  #geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  geom_point(data = comm_sum, aes(y = p_med, color = family), size = 3) +
  geom_point(data=subset(tab_p, species == "Trechus quadristriatus"), 
             aes(y=p_hat), color = "blue", size=2)+
  facet_wrap(order~.)+
  labs(
    y = "p̂ (détectabilité, occupancy unmarked;\n effort médian par méthode)",
    x = NULL,
    title = "Détectabilité par méthode —\najustée pour l’effort (unmarked::occu)"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "p_by_method_unmarked_effort_psi_cov.png"), p_fig, width = 8, height = 5, dpi = 200)

message("=== FIN ===  Résultats : ", normalizePath(out_dir))


################ DO TRAITS EXPLAIN DETECTABILITY ##############################
# Traits mesurés sur les animaux du dataset
trait <- raw %>%
  filter(PROJET == "RMQS_2025") %>%
  select(LB_NOM, `MASSE (mg)`) %>%
  rename(species = LB_NOM) %>%
  group_by(species) %>%
  summarise(mass = median(`MASSE (mg)`, na.rm = T))

dat <- tab_p %>%
  left_join(trait)

# traits issus de BETSI
devtools::install_git("https://src.koda.cnrs.fr/cefe/bioflux/betsir.git")
library(betsir)

h <- get_coded_value("Carabidae", "diet", include_subtaxa = TRUE)
############### /!\ à terminer #################


dat <- tab_p %>%
  left_join(trait)

ggplot(subset(dat, method %in% c("GPD", "Pitfall10", "Pitfall6", "Pitfall4")), 
       aes(x=log(mass), y=p_hat, colour=method))+
  geom_point()+
  stat_smooth(method='lm')+
  facet_grid(method~.)

# ===== Passer p_hat sur l’échelle logit + poids méta ====================
# Calcul SE(p) approx. depuis l'IC95% si possible : se_p ≈ (ucl - lcl) / (2 * 1.96)
# Puis delta-method : se_logit = se_p / (p * (1 - p))
# On tronque p pour éviter 0/1 exacts.

logit <- function(p) log(p/(1-p))

dat2 <- dat %>%
  mutate(
    p_clip  = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    logit_p = logit(p_clip),
    se_p    = ifelse(!is.na(lcl) & !is.na(ucl),
                     (ucl - lcl) / (2 * 1.96),
                     NA_real_),
    se_logit = ifelse(!is.na(se_p),
                      se_p / (p_clip * (1 - p_clip)),
                      NA_real_),
    weight   = 1 / (se_logit^2)
  )

# Si certains n'ont pas d'IC, on leur met un poids neutre (1)
dat2 <- dat2 %>%
  mutate(weight = ifelse(is.finite(weight), weight, 1))

# ===== 4) Ajuster un modèle : logit(p_hat) ~ method * log(mass) =============
# Comme mass est au niveau espèce et qu’on a répété par méthode,
# on met un intercept aléatoire par espèce pour la corrélation intra-espèce.

dat2 <- dat2 %>% filter(!is.na(mass) & mass > 0)

dat2 <- dat2 %>%
  mutate(
    log_mass = log10(mass),
    species  = factor(species)
  )

# lmer accepte des poids observationnels
m_mass_int <- lmer(logit_p ~ method * log_mass + (1 | species),
                   data = dat2, weights = weight, REML = TRUE)

summary(m_mass_int)

# ===== 5) Pentes par méthode et comparaisons (emmeans) =====================

# Pente de logit(p) par rapport à log_mass pour chaque méthode
trends <- emtrends(m_mass_int, ~ method, var = "log_mass")
trends
pairs(trends)  # comparer les pentes entre méthodes

# ===== 6) Visualisation simple =============================================

# Prédictions (lissées) par méthode sur une grille de masses
newgrid <- dat2 %>%
  group_by(method) %>%
  summarise(m_min = quantile(log_mass, 0.05),
            m_max = quantile(log_mass, 0.95), .groups="drop") %>%
  rowwise() %>%
  do({
    tibble(
      method  = .$method,
      log_mass = seq(.$m_min, .$m_max, length.out = 100)
    )
  }) %>%
  ungroup()

pred <- cbind(newgrid,
              predict(m_mass_int, newdata = newgrid, re.form = NA, se.fit = TRUE)) %>%
  as_tibble() %>%
  mutate(
    p_fit  = plogis(fit),
    p_lcl  = plogis(fit - 1.96 * se.fit),
    p_ucl  = plogis(fit + 1.96 * se.fit)
  )

dat2_r <- subset(dat2, method %in% c("GPD", "Pitfall10", "Pitfall6"))
pred_r <- subset(pred, method %in% c("GPD", "Pitfall10", "Pitfall6"))

ggplot(dat2_r,
       aes(x = log_mass, y = p_clip, colour = method)) +
  geom_point(alpha = 0.5) +
  geom_line(data = pred_r, aes(x = log_mass, y = p_fit, colour = method), linewidth = 1) +
  geom_ribbon(data = pred_r,
              aes(x = log_mass, ymin = p_lcl, ymax = p_ucl, fill = method),
              alpha = 0.15, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(method~.)+
  labs(x = "log10(body mass, mg)",
       y = "Detectability p̂\n (occupancy, logit back-transformed)",
       title = "Effect of body mass on detectability by method") +
  theme_minimal()+
  theme(legend.position = "none")


###################### DIversité ~##################################
# =========================
# Diversité spécifique par méthode
# + Profils de Hill q∈[0,2]
# =========================



# -------------------------
# 0) Données d'entrée
# -------------------------
# dat0 : site_id, method, species, abund  (déjà construit dans ton script précédent)
stopifnot(all(c("site_id","method","species","abund") %in% names(dat0)))

# (Option) Détectabilités par espèce×méthode (si dispo) pour correction HT
# Décommente si tu veux activer la correction
# p_tab <- readr::read_csv("data/derived-data/occupancy/species_detection_by_method_unmarked_effort.csv")
# head(p_tab)
# -> colonnes attendues : species, method, p_hat

# Normalise/filtre les méthodes que tu veux comparer
method_levels <- c("GPD","Pitfall")
dat_use <- dat0 %>%
  mutate(method = factor(method, levels = method_levels)) %>%
  filter(!is.na(method)) %>%
  # garde seulement les méthodes présentes
  filter(method %in% levels(method)[!is.na(levels(method))])

# -------------------------
# 1) Matrices site × espèce par méthode (abondances & présence)
# -------------------------
# Abondance (somme par site×espèce×méthode)
abund_long <- dat_use %>%
  group_by(site_id, method, species) %>%
  summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop")

# Présence/absence dérivée
pa_long <- abund_long %>%
  mutate(det = as.integer(abund > 0))

# Pour calculs par méthode, on construit des matrices site×espèce séparées par méthode
make_comm_matrix <- function(df_long, value = c("abund","det")) {
  val <- match.arg(value)
  df_wide <- df_long %>%
    select(site_id, method, species, !!sym(val)) %>%
    pivot_wider(names_from = species, values_from = !!sym(val), values_fill = 0)
  df_wide
}

comm_abund <- make_comm_matrix(abund_long, "abund")
comm_pa    <- make_comm_matrix(pa_long,    "det")

# -------------------------
# 2) Indices de diversité par site × méthode
# -------------------------
#Calculer les profils de diversité à l'aide de 0.1. divent_multitaxa
#extraire les valeurs de 3 indices

m_S  <- lmer(richness ~ method + (1|site_id), data = div_hill)
m_q1 <- lmer(Hill_q1 ~ method + (1|site_id), data = div_hill)
m_q2 <- lmer(Hill_q2 ~ method + (1|site_id), data = div_hill)

summary(m_S)
summary(m_q1)
summary(m_q2)

# Visualisation simple
ggplot(div_hill, aes(x = method, y = richness)) +
  geom_violin(fill = "grey92") +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  ggtitle("Richesse spécifique par méthode (abondances)") +
  theme_minimal()

ggplot(div_hill, aes(x = method, y = Hill_q1)) +
  geom_violin(fill = "grey92") +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  #ylab(expression(^1*D==e^{H})) +
  ggtitle("Diversité de Hill q=1 (exp(Shannon)) par méthode") +
  theme_minimal()

ggplot(div_hill, aes(x = method, y = Hill_q2)) +
  geom_violin(fill = "grey92") +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  #ylab(expression(^2*D)) +
  ggtitle("Diversité de Hill q=2 (Inverse Simpson) par méthode") +
  theme_minimal()

# -------------------------
# 5) (Option) Correction HT par détectabilité p̂ (espèce×méthode)
# -------------------------
# Si 'tab_p' est disponible : species, method, p_hat
# -> On corrige les abondances par 1/p_hat (borné), puis on refait les indices et Hill
#    (Rappel : c’est une correction grossière, p pas spécifique au site)

# Clean p_hat
   p_tab2 <- tab_p %>%
     select(species, method, p_hat) %>%
     mutate(method = factor(method, levels = method_levels)) %>%
     filter(!is.na(p_hat) & p_hat > 0) %>%
     mutate(p_hat = pmin(p_hat, 0.999))  # borne

   abund_corr <- abund_long %>%
     left_join(p_tab2, by = c("species","method")) %>%
     mutate(abund_corr = ifelse(!is.na(p_hat) & p_hat > 0,
                                abund / p_hat, NA_real_)) %>%
     # fallback: si p absent, on garde l’abondance observée (ou NA si tu préfères)
     mutate(abund_corr = ifelse(is.na(abund_corr), abund, abund_corr))

   comm_abund_corr <- abund_corr %>%
     select(site_id, method, species, abund_corr) %>%
     pivot_wider(names_from = species, values_from = abund_corr, values_fill = 0)

   div_abund_corr <- calc_div(comm_abund_corr, is_abundance = TRUE) %>%
     mutate(Hill_q1 = exp(shannon), Hill_q2 = simpson)

   # Compare modèles “observé” vs “corrigé”
   m_S_corr  <- lmer(richness ~ method + (1|site_id), data = div_abund_corr)
   m_q1_corr <- lmer(Hill_q1 ~ method + (1|site_id), data = div_abund_corr)
   m_q2_corr <- lmer(Hill_q2 ~ method + (1|site_id), data = div_abund_corr)

   summary(m_S_corr)
   summary(m_q1_corr)
   summary(m_q2_corr)

   # Profils de Hill “corrigés”
   hill_profiles_corr <- comm_abund_corr %>%
     group_by(site_id, method) %>%
     group_modify(~{
       abund_vec <- as.numeric(.x %>% select(-1:2))
       den <- sum(abund_vec)
       if (den == 0) tibble(q = q_grid, hill = 0) else {
         pvec <- abund_vec / den
         tibble(q = q_grid, hill = sapply(q_grid, function(q) hill_D(pvec, q)))
       }
     }) %>%
     ungroup() %>%
     group_by(method, q) %>%
     summarise(hill_med = median(hill, na.rm=TRUE),
               hill_lo  = quantile(hill, 0.25, na.rm=TRUE),
               hill_hi  = quantile(hill, 0.75, na.rm=TRUE),
               .groups = "drop")

   ggplot(hill_profiles_corr, aes(x=q, y=hill_med, color=method)) +
     geom_ribbon(aes(ymin = hill_lo, ymax = hill_hi, fill = method), alpha = 0.15, color = NA) +
     geom_line(linewidth = 1) +
     scale_x_continuous("q (ordre de Hill)", breaks = seq(0,2,0.5), limits = c(0,2)) +
     #ylab(expression(^q*D)) +
     ggtitle("Profils de diversité (Hill) par méthode — abondances corrigées par p̂") +
     theme_minimal()

