# ============================================================
# OCCUPANCY espèce-par-espèce (unmarked::occu) avec effort
# ψ ~ ALTITUDE_z + DOY_z ; p ~ méthode + effort_z
# Méthodes sélectionnables : Pitfall, GPD, DVAC, TRI
# + Composition (Venn/PCoA/PERMANOVA/betapart) seulement pour methods_use
# ============================================================

## 0) Paramètres ---------------------------------------------------------------
in_path        <- "data/raw-data/lab files/databases update/update_carabidae.csv"
out_dir        <- "data/derived-data/occupancy"
target_project <- "RMQS_2025"
taxo_level     <- "ES"
min_sites      <- 2
require_contrast <- TRUE

# ---- Choix des méthodes à étudier (seules celles listées seront traitées)
methods_use <- c("Pitfall","GPD")     # <--- édite ici (parmi "Pitfall","GPD","DVAC","TRI")

# Sous-ensembles Pitfall à comparer (réplicats)
pit_reps <- c(10, 8, 6, 4, 2)         # ignoré si "Pitfall" n’est pas dans methods_use

# Constantes effort
pit_diam_m         <- 0.05
pit_circ_m         <- pi * pit_diam_m
gpd_open_circ_m    <- pi * pit_diam_m
gpd_cross_m        <- 1.0
gpd_unit_intensity <- gpd_open_circ_m + gpd_cross_m
tm_unit            <- 0.0625

## 1) Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  librarian::shelf(
    dplyr, tidyr, ggplot2, readr, stringr, tibble, purrr,
    unmarked, vegan, betapart, ggvenn
  )
})
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 2) Lecture & colonnes -------------------------------------------------------
raw <- readr::read_delim(in_path, delim = ";",
                         locale = readr::locale(encoding = "Latin1"),
                         show_col_types = FALSE)

col_exists <- function(x) x[x %in% names(raw)][1]
col_station <- col_exists(c("NOM_STATION"))
col_method  <- col_exists(c("METHODE"))
col_species <- col_exists(c("LB_NOM"))
col_abund   <- col_exists(c("ABONDANCE_TOTALE"))
col_sample  <- col_exists(c("ID_ECHANTILLON"))
col_project <- col_exists(c("PROJET"))
col_DATE    <- col_exists(c("DATE"))
col_DATEFIN <- col_exists(c("DATE_FIN"))
col_alt     <- col_exists(c("ALTITUDE"))
col_rang    <- col_exists(c("RANG"))
stopifnot(!is.na(col_station), !is.na(col_method), !is.na(col_species),
          !is.na(col_abund), !is.na(col_sample), !is.na(col_rang))

## 3) Normalisation / parsing --------------------------------------------------
norm_method <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x, "(?i)gpd|directionnel|croisillon") ~ "GPD",
    stringr::str_detect(x, "(?i)barb|pitfall|piege\\s*barber") ~ "Pitfall",
    stringr::str_detect(x, "(?i)dvac|aspir") ~ "DVAC",
    stringr::str_detect(x, "(?i)tri\\s*manuel|tri\\s*manu") ~ "TRI",
    TRUE ~ x
  )
}
parse_abund   <- function(x) suppressWarnings(as.numeric(x))
parse_date_fr <- function(x) as.Date(suppressWarnings(readr::parse_date(as.character(x), "%d/%m/%Y")))
to_doy        <- function(d) ifelse(is.na(d), NA, as.integer(format(d, "%j")))
extract_last_int <- function(x) suppressWarnings(readr::parse_integer(stringr::str_extract(x, "(\\d+)\\s*$")))

dat0 <- raw %>%
  filter(.data[[col_project]] %in% target_project,
         .data[[col_rang]]    %in% taxo_level) %>%
  transmute(
    site_id = as.character(.data[[col_station]]),
    method  = norm_method(.data[[col_method]]),
    species = as.character(.data[[col_species]]),
    abund   = parse_abund(.data[[col_abund]]),
    sample  = as.character(.data[[col_sample]]),
    DATE     = parse_date_fr(.data[[col_DATE]]),
    DATE_FIN = parse_date_fr(.data[[col_DATEFIN]]),
    ALTITUDE = suppressWarnings(as.numeric(.data[[col_alt]]))
  ) %>%
  filter(method %in% methods_use)  # <- filtrage fort ici

stopifnot(nrow(dat0) > 0)

# identifiants pièges
dat0 <- dat0 %>%
  mutate(
    trap    = ifelse(method == "Pitfall",
                     suppressWarnings(readr::parse_integer(stringr::str_extract(sample, "(?<=PB)\\d+$"))),
                     NA_integer_),
    gpd_sub = ifelse(method == "GPD", extract_last_int(sample), NA_integer_),
    det     = as.integer(abund > 0),
    days    = {
      d <- as.numeric(DATE_FIN - DATE) + 1
      ifelse(is.na(d) | d < 1, 1, d)
    }
  )

## 4) Détections (dynamiques selon methods_use) --------------------------------
det_tabs <- list()

# Pitfall (réplicats simulés)
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

# GPD (unique réplicat)
if ("GPD" %in% methods_use) {
  det_tabs$GPD <- dat0 %>% filter(method == "GPD") %>%
    group_by(site_id, species) %>% summarise(GPD = as.integer(any(det==1, na.rm=TRUE)), .groups="drop")
}

# DVAC
if ("DVAC" %in% methods_use) {
  det_tabs$DVAC <- dat0 %>% filter(method == "DVAC") %>%
    group_by(site_id, species) %>% summarise(DVAC = as.integer(any(det==1, na.rm=TRUE)), .groups="drop")
}

# TRI
if ("TRI" %in% methods_use) {
  det_tabs$TM <- dat0 %>% filter(method == "TRI") %>%
    group_by(site_id, species) %>% summarise(TM = as.integer(any(det==1, na.rm=TRUE)), .groups="drop")
}

## 5) Effort (dynamiques) ------------------------------------------------------
effort_parts <- list(site_id = tibble(site_id = sort(unique(dat0$site_id))))

# Pitfall
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

# GPD (nb sous-pièges actifs)
if ("GPD" %in% methods_use) {
  eff_gpd_site <- dat0 %>% filter(method == "GPD") %>%
    group_by(site_id) %>%
    summarise(
      days_gpd = if (all(is.na(days))) 1 else stats::median(days, na.rm = TRUE),
      n_active_subtraps = n_distinct(gpd_sub[!is.na(gpd_sub)]),  # 0..4
      .groups = "drop"
    ) %>%
    mutate(
      n_active_subtraps = ifelse(is.na(n_active_subtraps), 0L, n_active_subtraps),
      eff_gpd = gpd_unit_intensity * days_gpd * pmax(0, n_active_subtraps)
    ) %>% select(site_id, eff_gpd, n_active_subtraps)
  effort_parts$gpd <- eff_gpd_site
}

# DVAC
if ("DVAC" %in% methods_use) {
  effort_parts$dvac <- dat0 %>% filter(method == "DVAC") %>% distinct(site_id) %>% mutate(eff_dvac = 1)
}

# TRI
if ("TRI" %in% methods_use) {
  effort_parts$tm <- dat0 %>% filter(method == "TRI") %>% distinct(site_id) %>% mutate(eff_tm = 6 * tm_unit)
}

# Fusion dynamique des efforts
effort_wide <- reduce(effort_parts, ~left_join(.x, .y, by = "site_id"))
effort_wide <- effort_wide %>% mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))

## 6) Méthodes présentes par site (sans colonnes fantômes) --------------------
sampling_map <- dat0 %>%
  distinct(site_id, method) %>%
  tidyr::pivot_wider(
    names_from  = method,
    values_from = method,
    values_fill = NA,
    names_expand = FALSE
  ) %>% mutate(across(-site_id, ~ as.integer(!is.na(.))))

# Ajout conditionnel de has_gpd selon effort (si GPD utilisé)
if ("GPD" %in% methods_use) {
  sampling_map <- sampling_map %>%
    left_join(effort_wide %>% select(site_id, eff_gpd, n_active_subtraps), by="site_id") %>%
    mutate(GPD = ifelse(eff_gpd > 0 & n_active_subtraps > 0, 1L, 0L)) %>%
    select(-eff_gpd, -n_active_subtraps)
}

## 7) Détections + NA si méthode absente ---------------------------------------
keys <- dat0 %>% distinct(site_id, species)
det_wide0 <- keys
if (length(det_tabs)) det_wide0 <- reduce(det_tabs, 
                                         ~left_join(.x, .y, by=c("site_id","species")))

# Toutes colonnes présentes actuellement (= seules méthodes utilisées)
det_cols <- setdiff(names(det_wide0), c("site_id","species"))
# Remplacer NA par 0 (détection), puis remettre NA si méthode absente au site
det_wide <- det_wide0 %>%
  mutate(across(all_of(det_cols), ~ tidyr::replace_na(., 0L))) #%>%
  #left_join(sampling_map, by="site_id")

# pour chaque colonne méthode, masque si absent au site
# Masquer chaque colonne de détection si la méthode est absente au site
# - det_cols : noms de colonnes dans det_wide à masquer (ex. "Pitfall10","Pitfall6","GPD", ...)
stopifnot(exists("det_cols"))
for (cl in det_cols) {
  # ignorer si la colonne n'existe pas (méthode non sélectionnée / non créée)
  if (!cl %in% names(det_wide)) next
  
  # "Pitfall10" -> "Pitfall" pour retrouver la clé de présence dans sampling_map
  method_name <- sub("^Pitfall\\d+$", "Pitfall", cl)
  
  # ignorer si sampling_map n'a pas cette méthode (méthode non sélectionnée)
  if (!method_name %in% names(sampling_map)) next
  
  # aligne la présence méthode/site sur l'ordre de det_wide
  idx <- match(det_wide$site_id, sampling_map$site_id)
  present_vec <- sampling_map[[method_name]][idx]
  present_vec[is.na(present_vec)] <- 0L
  
  # vecteur de détection ; force en entier pour éviter les soucis de type
  x <- det_wide[[cl]]
  x <- as.integer(x)
  
  # masque : si méthode absente (0), mettre NA
  x[present_vec != 1L] <- NA_integer_
  det_wide[[cl]] <- x
}

det_wide <- det_wide %>% select(site_id, species, all_of(det_cols))

## 8) Y et Effort (uniquement colonnes existantes) -----------------------------
sites <- sort(unique(det_wide$site_id))
spp   <- sort(unique(det_wide$species))
Rnames <- det_cols  # réplicats réellement présents
J <- length(sites); N <- length(spp); R <- length(Rnames)

Y <- array(NA_integer_, dim = c(J, R, N),
           dimnames = list(site = sites, rep = Rnames, species = spp))

det_long <- det_wide %>%
  mutate(site_id = factor(site_id, levels = sites),
         species = factor(species, levels = spp)) %>%
  arrange(species, site_id)

for (k in seq_along(spp)) {
  dk <- det_long %>% filter(species == spp[k]) %>% arrange(site_id)
  if (nrow(dk) < J) dk <- tibble(site_id = factor(sites, levels = sites)) %>% left_join(dk, by="site_id")
  for (r in seq_along(Rnames)) Y[, r, k] <- as.numeric(dk[[Rnames[r]]])
}

# Effort
eff_mat <- matrix(0, nrow = J, ncol = R, dimnames = list(sites, Rnames))
map_eff_name <- c(setNames(paste0("eff_p", pit_reps), paste0("Pitfall", pit_reps)),
                  "GPD" = "eff_gpd", "DVAC" = "eff_dvac", "TM" = "eff_tm")
for (rn in Rnames) {
  nm <- map_eff_name[rn]
  if (!is.na(nm) && nm %in% names(effort_wide)) {
    eff_mat[, rn] <- effort_wide[[nm]][match(sites, effort_wide$site_id)]
  }
}

# Standardisation effort (log1p + z) et masque all-NA
eff_vec <- as.vector(eff_mat); eff_log <- log1p(eff_vec)
eff_mean <- mean(eff_log, na.rm = TRUE); eff_sd <- sd(eff_log, na.rm = TRUE); if (!is.finite(eff_sd) || eff_sd==0) eff_sd <- 1
eff_z <- matrix((eff_log - eff_mean)/eff_sd, nrow = J, ncol = R, byrow = FALSE, dimnames = dimnames(eff_mat))
mask_allNA_sr <- apply(Y, c(1,2), function(a) all(is.na(a)))
eff_z[mask_allNA_sr] <- NA

## 9) Covariables site ---------------------------------------------------------
site_alt <- dat0 %>% group_by(site_id) %>% summarise(ALTITUDE = median(ALTITUDE, na.rm=TRUE), .groups="drop")
site_doy <- dat0 %>% group_by(site_id) %>% summarise(DOY = median(to_doy(DATE),  na.rm=TRUE), .groups="drop")
zstd <- function(v){ m <- mean(v,na.rm=TRUE); s <- sd(v,na.rm=TRUE); z <- (v-m)/ifelse(s>0,s,1); z[is.na(z)] <- 0; z }
site_cov <- tibble(site_id = sites) %>% left_join(site_alt, by="site_id") %>% left_join(site_doy, by="site_id") %>%
  mutate(ALTITUDE_z = zstd(ALTITUDE), DOY_z = zstd(DOY)) %>% select(site_id, ALTITUDE_z, DOY_z)

## 10) Fit unmarked::occu ------------------------------------------------------
method_eff_medians <- function(eff_z_matrix)
  apply(eff_z_matrix, 2, function(col) if (all(is.na(col))) 0 else median(col, na.rm=TRUE))

fit_one_species_occu <- function(sp_name) {
  k <- which(dimnames(Y)$species == sp_name); if (length(k) != 1) return(NULL)
  Yk <- Y[,,k, drop=FALSE][,,1]
  keep_site <- apply(Yk, 1, function(a) any(!is.na(a))); Yk <- Yk[keep_site,, drop=FALSE]
  if (nrow(Yk) < min_sites) return(NULL)
  keep_rep  <- apply(Yk, 2, function(a) any(!is.na(a))); Yk <- Yk[, keep_rep, drop=FALSE]
  if (ncol(Yk) < 2) return(NULL)
  yvec <- as.vector(Yk[!is.na(Yk)]); if (require_contrast && !(any(yvec==1) && any(yvec==0))) return(NULL)
  
  rep_names_k <- colnames(Yk)
  mkI <- function(tag) matrix(as.numeric(rep_names_k == tag), nrow(Yk), ncol(Yk), byrow=TRUE)
  I_list <- setNames(lapply(rep_names_k, mkI), rep_names_k)
  eff_sp <- eff_z[rownames(Yk), rep_names_k, drop=FALSE]
  for (nm in rep_names_k) I_list[[nm]][is.na(Yk)] <- NA
  eff_sp[is.na(Yk)] <- NA
  
  sc <- site_cov %>% filter(site_id %in% rownames(Yk)) %>% arrange(match(site_id, rownames(Yk))) %>% select(ALTITUDE_z, DOY_z)
  has_var_alt <- sd(sc$ALTITUDE_z, na.rm=TRUE) > 0; has_var_doy <- sd(sc$DOY_z, na.rm=TRUE) > 0
  state_terms <- c(); if (has_var_alt) state_terms <- c(state_terms,"ALTITUDE_z"); if (has_var_doy) state_terms <- c(state_terms,"DOY_z")
  state_rhs <- if (length(state_terms)) paste(state_terms, collapse=" + ") else "1"
  
  obsCovs <- c(I_list, list(eff_z = eff_sp))
  umf <- unmarked::unmarkedFrameOccu(y = Yk, siteCovs = sc, obsCovs = obsCovs)
  
  det_terms <- paste(c(rep_names_k, "eff_z"), collapse = " + ")
  form_occu <- as.formula(paste("~", det_terms, "~", state_rhs))
  
  fm <- try(unmarked::occu(form_occu, data = umf, control = list(maxit = 200)), silent = TRUE)
  if (inherits(fm, "try-error")) return(NULL)
  
  eff_med_all <- method_eff_medians(eff_z)
  eff_for <- function(m) if (m %in% names(eff_med_all)) eff_med_all[[m]] else 0
  mk_pred <- function(tag) {
    nd <- data.frame(eff_z = eff_for(tag))
    for (nm in rep_names_k) nd[[nm]] <- as.numeric(nm == tag)
    if (has_var_alt) nd$ALTITUDE_z <- 0
    if (has_var_doy) nd$DOY_z      <- 0
    as.data.frame(unmarked::predict(fm, type = "det", newdata = nd))
  }
  pred_list <- lapply(rep_names_k, mk_pred)
  
  co_state <- summary(fm)$state
  getcoef <- function(nm) if (nm %in% rownames(co_state)) co_state[nm,"Estimate"] else NA_real_
  getse   <- function(nm) if (nm %in% rownames(co_state)) co_state[nm,"SE"]       else NA_real_
  
  list(
    p_table = tibble::tibble(
      species = sp_name, method = rep_names_k,
      p_hat = vapply(pred_list, function(z) z$Predicted, numeric(1)),
      lcl   = vapply(pred_list, function(z) z$lower,     numeric(1)),
      ucl   = vapply(pred_list, function(z) z$upper,     numeric(1)),
      n_sites = nrow(Yk), n_rep = ncol(Yk)
    ),
    beta = tibble::tibble(
      species  = sp_name,
      beta_alt = getcoef("ALTITUDE_z"), se_alt = getse("ALTITUDE_z"),
      beta_doy = getcoef("DOY_z"),      se_doy = getse("DOY_z")
    ),
    psi_hat = tibble::tibble(
      site_id = rownames(Yk), species = sp_name,
      psi_hat = as.data.frame(unmarked::predict(fm, type="state"))$Predicted
    )
  )
}

spp_all  <- dimnames(Y)$species
res_list <- purrr::map(spp_all, ~try(fit_one_species_occu(.x), silent = TRUE))

tab_p   <- res_list %>% keep(~!inherits(.x,"try-error") && !is.null(.x$p_table)) %>% map("p_table") %>% bind_rows()
tab_beta<- res_list %>% keep(~!inherits(.x,"try-error") && !is.null(.x$beta))     %>% map("beta")     %>% bind_rows()
tab_psi <- res_list %>% keep(~!inherits(.x,"try-error") && !is.null(.x$psi_hat))  %>% map("psi_hat")  %>% bind_rows()

readr::write_csv(tab_p,   file.path(out_dir, "p_hat_by_method_unmarked_full.csv"))
readr::write_csv(tab_beta,file.path(out_dir, "psi_coef_ALT_DOY_by_species.csv"))
readr::write_csv(tab_psi, file.path(out_dir, "psi_hat_by_site_species.csv"))

## 11) Figure p(method) --------------------------------------------------------
if (nrow(tab_p) > 0) {
  comm_sum <- tab_p %>% group_by(method) %>%
    summarise(p_med = median(p_hat, na.rm=TRUE), n_sp = n_distinct(species), .groups="drop")
  tab_p$method <- factor(tab_p$method, levels = unique(tab_p$method))
  p_fig <- ggplot(tab_p, aes(x = method, y = p_hat)) +
    geom_violin(fill = "grey92") +
    geom_point(aes(colour = species), position = position_jitter(width = 0.08), alpha = 0.5) +
    geom_point(data = comm_sum, aes(y = p_med), color = "red", size = 3) +
    labs(y = "p̂ (détectabilité, occupancy; effort médian par méthode)",
         x = NULL, title = "Détectabilité par méthode — ajustée pour l’effort") +
    theme_minimal()
  ggsave(file.path(out_dir, "p_by_method_unmarked_effort.png"), p_fig, width = 9, height = 5, dpi = 200)
}

## 12) Synthèse communautaire & figures ----------------------------------------
# Détectabilité par méthode
comm_sum <- tab_p %>%
  group_by(method) %>%
  summarise(p_med = median(p_hat, na.rm=TRUE),
            p_mean = mean(p_hat, na.rm=TRUE),
            lcl_med = median(lcl, na.rm=TRUE),
            ucl_med = median(ucl, na.rm=TRUE),
            n_sp = n_distinct(species), .groups="drop")
print(comm_sum)
readr::write_csv(comm_sum, file.path(out_dir, "community_detection_unmarked_effort_psi_cov.csv"))

# Distribution des effets ALT & DOY sur ψ (par espèce)
tab_beta <- tab_beta %>%
  distinct(species, beta_alt, se_alt, beta_doy, se_doy)

readr::write_csv(tab_beta, file.path(out_dir, "psi_coefficients_by_species.csv"))

# Petits graphiques de synthèse
g1 <- ggplot(tab_beta, aes(x = beta_alt)) + geom_density(na.rm=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Effect on occupancy (β_ALTITUDE_z)", y="Density", title="Community distribution of altitude effects on ψ") +
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
  #facet_wrap(order+family~.)+
  labs(
    y = "p̂ (détectabilité, occupancy unmarked;\n effort médian par méthode)",
    x = NULL,
    title = "Détectabilité par méthode —\najustée pour l’effort (unmarked::occu)"
  ) +
  theme_minimal()

ggplot(tab_p, aes(x = method, y = p_hat)) +
  geom_violin(fill = "grey92") +
  #geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  geom_point(data = comm_sum, aes(y = p_med), size = 3) +
  geom_point(data=subset(tab_p, species == "Trechus quadristriatus"), 
             aes(y=p_hat), color = "blue", size=2)+
  labs(
    y = "p̂ (détectabilité, occupancy unmarked;\n effort médian par méthode)",
    x = NULL,
    title = "Détectabilité par méthode —\najustée pour l’effort (unmarked::occu)"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "p_by_method_unmarked_effort_psi_cov.png"), p_fig, width = 8, height = 5, dpi = 200)

message("=== FIN ===  Résultats : ", normalizePath(out_dir))


# =====================================================================
# PARTIE B — DO TRAITS EXPLAIN DETECTABILITY
# =====================================================================

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

ggplot(dat, 
       aes(x=log(mass), y=p_hat, colour=method))+
  geom_point()+
  stat_smooth(method='lm')+
  labs(
    y = "p̂\n (détectabilité, occupancy unmarked;\n effort médian par méthod et par espèce)",
    x = "Masse moyenne des espèces (mg) \n-échelle logarithmique-",
    title = "Effet de la masse moyenne des espèces \nsur leur détectabilité par méthode —\najustée pour l’effort (unmarked::occu)"
  )+
  facet_wrap(method~.)+
  theme_bw()

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



# =====================================================================
# PARTIE C — COMPOSITION (UNIQUEMENT methods_use)
# =====================================================================

# Abondance site×méthode×espèce (filtrée methods_use dès le début)
abund_long <- dat0 %>%
  group_by(site_id, method, species) %>%
  summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop")

# Présence/absence
pa_long <- abund_long %>% mutate(det = as.integer(abund > 0))

# Matrice site×espèce par méthode (chaque ligne = un site×méthode)
comm_pa <- pa_long %>%
  select(site_id, method, species, det) %>%
  tidyr::pivot_wider(names_from = species, values_from = det, values_fill = 0)

# ----- Venn (si ≤4 méthodes et site par site) -----
if (length(methods_use) <= 4) {
  methods_for_venn <- methods_use
  site_list <- sort(unique(dat0$site_id))[1:min(12, length(unique(dat0$site_id)))]
  for (sid in site_list) {
    sets <- lapply(methods_for_venn, function(m) {
      pa_long %>% filter(site_id == sid, method == m, det == 1L) %>% pull(species) %>% unique()
    })
    names(sets) <- methods_for_venn
    if (all(lengths(sets) == 0)) next
    p_venn <- ggvenn::ggvenn(sets, fill_alpha = .3, stroke_size = .5) + ggtitle(paste("Venn —", sid))
    ggsave(file.path(out_dir, paste0("venn_", gsub("[^A-Za-z0-9]+","_",sid), ".png")),
           p_venn, width = 5, height = 4, dpi = 200)
  }
}

# ----- PCoA + PERMANOVA (strata = site) -----
comm_pa_all <- comm_pa %>% relocate(site_id, method)
X <- comm_pa_all %>% select(-site_id, -method) %>% as.matrix()
if (nrow(X) >= 3 && ncol(X) >= 1) {
  d_jac <- vegan::vegdist(X, method = "jaccard", binary = TRUE)
  pco <- stats::cmdscale(d_jac, k = 2, eig = TRUE)
  scores <- as.data.frame(pco$points); colnames(scores) <- c("PCoA1","PCoA2")
  scores$site_id <- comm_pa_all$site_id; scores$method <- comm_pa_all$method
  p_pcoa <- ggplot(scores, aes(PCoA1, PCoA2, color = method)) +
    geom_point(size = 2, alpha = .8) +
    stat_ellipse(type = "norm", level = .67, linewidth = .6, alpha = .4) +
    theme_minimal() + ggtitle("PCoA (Jaccard, présence/absence) — coloré par méthode")
  ggsave(file.path(out_dir, "PCoA_methods.png"), p_pcoa, width = 7, height = 5, dpi = 200)
  
  adon <- vegan::adonis2(d_jac ~ method, data = comm_pa_all,
                         permutations = 9999, strata = comm_pa_all$site_id, by = "margin")
  capture.output(adon, file = file.path(out_dir, "PERMANOVA_methods.txt"))
}

# ----- Turnover vs Nestedness (betapart) -----
pair_turnover_nestedness_by_site <- function(sid, methods_vec) {
  df <- pa_long %>% filter(site_id == sid, method %in% methods_vec) %>%
    select(method, species, det) %>% pivot_wider(names_from = species, values_from = det, values_fill = 0)
  if (nrow(df) < 2) return(NULL)
  rownames(df) <- df$method
  B  <- betapart::betapart.core(as.data.frame(df %>% select(-method)))
  bp <- betapart::beta.pair(B, index.family = "jaccard")
  mat_to_long <- function(M, comp) {
    as.data.frame(as.table(as.matrix(M))) %>%
      setNames(c("m1","m2","value")) %>%
      filter(as.character(m1) < as.character(m2)) %>%
      mutate(component = comp, site_id = sid)
  }
  dplyr::bind_rows(
    mat_to_long(bp$beta.jtu, "turnover"),
    mat_to_long(bp$beta.jne, "nestedness"),
    mat_to_long(bp$beta.jac, "total")
  )
}
beta_pairs <- purrr::map_dfr(sort(unique(dat0$site_id)),
                             ~ pair_turnover_nestedness_by_site(.x, methods_use))
if (nrow(beta_pairs) > 0) {
  beta_summary <- beta_pairs %>% group_by(m1, m2, component) %>%
    summarise(median = median(value, na.rm=TRUE),
              q25    = quantile(value, .25, na.rm=TRUE),
              q75    = quantile(value, .75, na.rm=TRUE),
              n_site = dplyr::n(), .groups="drop")
  readr::write_csv(beta_pairs,  file.path(out_dir, "beta_pairs_by_site_turnover_nestedness.csv"))
  readr::write_csv(beta_summary,file.path(out_dir, "beta_pairs_summary_turnover_nestedness.csv"))
  p_beta <- beta_summary %>%
    mutate(pair = paste(m1, "vs", m2)) %>%
    ggplot(aes(x = pair, y = median, fill = component)) +
    geom_col(position = position_dodge(width=.8)) +
    geom_errorbar(aes(ymin = q25, ymax = q75),
                  position = position_dodge(width=.8), width = .2) +
    coord_flip() +
    labs(x = "Paire de méthodes", y = "β (Jaccard) — médiane [IQR]",
         title = "Partition de la β-diversité : turnover vs nestedness (entre méthodes, au sein des sites)") +
    theme_minimal()
  ggsave(file.path(out_dir, "beta_partition_methods.png"), p_beta, width = 8, height = 5, dpi = 200)
}

message("=== FIN ===  Résultats : ", normalizePath(out_dir))
