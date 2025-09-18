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
  filter(method %in% methods_use)

stopifnot(nrow(dat0) > 0)

dat0 <- dat0 %>%
  mutate(
    trap    = ifelse(method == "Pitfall",
                     suppressWarnings(readr::parse_integer(stringr::str_extract(sample, "(?<=PB)\\d+$"))),
                     NA_integer_),
    gpd_sub = ifelse(method == "GPD", extract_last_int(sample), NA_integer_),
    det     = as.integer(abund > 0),
    days    = { d <- as.numeric(DATE_FIN - DATE) + 1; ifelse(is.na(d) | d < 1, 1, d) }
  )
