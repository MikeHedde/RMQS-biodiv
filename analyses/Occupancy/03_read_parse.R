# ============================================================
# LECTURE & PARSING
# Sorties: raw, dat0 (site_id, method, species, abund, sample, DATE, DATE_FIN, ALTITUDE)
# ============================================================

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

dat00 <- raw %>%
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
  filter(method %in% methods_use)

stopifnot(nrow(dat00) > 0)

dat0 <- dat00 %>%
  mutate(
    trap    = ifelse(method == "Pitfall",
                     suppressWarnings(readr::parse_integer(stringr::str_extract(sample, "(?<=PB)\\d+$"))),
                     NA_integer_),
    gpd_sub = ifelse(method == "GPD", extract_last_int(sample), NA_integer_),
    dvac = ifelse(method == "DVAC", extract_last_int(sample), NA_integer_),
    det     = as.integer(abund > 0),
    days    = { d <- as.numeric(DATE_FIN - DATE) + 1; ifelse(is.na(d) | d < 1, 1, d) }
  )

##################################################
# Manage habitat names
##################################################
# === HABITAT -> Open / Closed  ===============

# ==== HABITAT -> Open / Closed (adapté à ta liste) ==========================
# Dépend : stringr, stringi
if (!requireNamespace("stringi", quietly = TRUE)) install.packages("stringi")

classify_habitat_open_closed <- function(x, unknown = NA_character_) {
  # Normalisation : minuscules, espaces, sans accents
  y <- tolower(trimws(as.character(x)))
  y <- stringr::str_replace_all(y, "\\s+", " ")
  y <- stringi::stri_trans_general(y, "Latin-ASCII")  # hetraie, chenaie, foret, etc.
  
  res <- rep(unknown, length(y))
  
  # ---- Overrides (demandés) -> Closed ----
  idx_force_closed <- stringr::str_detect(
    y,
    paste0(
      "\\bchenaie\\s*claire\\b",                                   # Chênaie claire
      "|\\bancienne\\s*chataigneraie\\s*sur\\s*terrasse\\b"        # Ancienne châtaigneraie sur terrasse
    )
  )
  res[idx_force_closed] <- "Closed"
  
  # ---- Fermé (forêts/boisements) ----
  # couvre : Forêt/mixte/feuillus ; Hêtraie(-sapinière, -Épicéa) ; Sapinière ; Épicéa ;
  # Chênaie (y compris verte, claire déjà forcée ci-dessus) ; Châtaigneraie ; Érablaie ; combinaisons
  re_closed <- paste(
    "\\bforet(s)?\\b",                         # Foret / Forets
    "foret\\s*mixte", "foret\\s*feuillus",
    "\\bhetraie\\b", "hetraie-?sapiniere", "hetraie-?epicea",
    "\\bsapiniere\\b", "\\bepicea(s)?\\b",
    "\\bchenaie(?!\\s*claire)\\b",             # 'claire' déjà forcée
    "chenaie\\s*verte",
    "\\bchataigneraie\\b",
    "\\berablaie\\b",
    "chenaie\\s*/\\s*erablaie",                # "Chênaie/ Érablaie"
    sep = "|"
  )
  
  # ---- Ouvert (prairies, pelouses, landes, friches, cultures, vignes, etc.) ----
  # couvre : Prairie (de fauche, permanente, rase, pâturée/pâture), Pelouse (alpine, écorchée...),
  # Lande (à bruyère/fougeraie/rhododendrons), Friche (humide), Jachère,
  # Cultures/Monoculture/Grande culture/Labour + cultures nommées,
  # Verger, Vigne(s), Oliveraie, Parc urbain, 'ouvert'
  re_open <- paste(
    "prairie(\\s*(de\\s*fauche|permanente|rase|de\\s*pature|paturee|pature))?",
    "\\bpelouse\\b(\\s*(alpine|ecorchee\\s*a\\s*salicorne))?",
    "lande(\\s*a\\s*bruyere\\s*/\\s*fougeraie|\\s*a\\s*rhododendrons)?",
    "friche(\\s*humide)?", "\\bjachere\\b",
    "grande\\s*culture", "\\blabour\\b",
    "monoculture(\\s*\\((ble|colza|feverolle)\\))?",
    "\\bculture\\b(\\s*de\\s*(ble|coriande))?",
    "\\bvigne(s)?\\b", "\\bverger(s)?\\b", "\\boliveraie(s)?\\b",
    "parc\\s*urbain",
    "\\bouvert(e)?\\b",
    sep = "|"
  )
  
  # Application hiérarchique : (1) overrides Closed, (2) closed, (3) open
  idx_closed <- stringr::str_detect(y, re_closed)
  res[idx_closed & is.na(res)] <- "Closed"
  
  idx_open <- stringr::str_detect(y, re_open)
  res[idx_open & is.na(res)] <- "Open"
  
  # Retour : facteur ordonné
  factor(res, levels = c("Closed", "Open"))
}


# utilisation
raw$HABITAT_OC <- classify_habitat_open_closed(raw$HABITAT)

# ---- HABITAT_OC (site-level) puis jointure sur dat0 ----
# 1) valeur "canonique" par site (modalité la plus fréquente)
# après avoir calculé raw$HABITAT_OC
hab_site <- raw %>%
  filter(.data[[col_project]] %in% target_project,
         .data[[col_rang]]    %in% taxo_level) %>%
  transmute(
    site_id    = as.character(.data[[col_station]]),
    HABITAT_OC = as.character(.data[["HABITAT_OC"]])
  ) %>%
  filter(!is.na(HABITAT_OC)) %>%
  count(site_id, HABITAT_OC, name = "n") %>%
  group_by(site_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site_id, HABITAT_OC)

# avant de joindre, s’assurer qu’aucune ancienne colonne ne traîne :
dat0 <- dat0 %>% select(-matches("^HABITAT_OC(\\.|$)"), -any_of("HABITAT_OC"))

dat0 <- dat0 %>%
  left_join(hab_site, by = "site_id") %>%
  mutate(HABITAT_OC = factor(HABITAT_OC, levels = c("Closed","Open")))
