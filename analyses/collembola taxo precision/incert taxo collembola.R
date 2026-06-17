
# taxonomic_uncertainty_collembola_v4_divent.R
# ------------------------------------------------------------
# Objectif : tester comment l'incertitude taxonomique modifie les inférences
# écologiques dans un jeu de données de collemboles.
#
# 
# Fichiers attendus dans le working directory :
#   collembola.csv
#   all_env_variables.csv
#   taxref_collembola.csv
#
# Packages :
# install.packages(c("tidyverse", "janitor", "vegan", "broom", "patchwork", "divent"))
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(vegan)      # distances, PCoA/Procrustes, PERMANOVA
  library(divent)     # Hill numbers + sample coverage
  library(broom)
  library(patchwork)
})

# -----------------------------
# 0. PARAMETRES
# -----------------------------

PATH_COLLEMBOLA <- "data/raw-data/1.faune/collembola.csv"
PATH_ENV        <- "data/derived-data/all_env_variables.csv"
PATH_TAXREF     <- "data/raw-data/taxref_collembola.csv"

OUT_DIR <- "outputs_taxonomic_uncertainty_v9"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


set.seed(123)

ANALYSIS_UNIT <- "station"   # recommandé avec drivers environnementaux stationnels
N_SIM <- 10

# PERMANOVA rapide, univariée par driver
RUN_DRIVER_BLOCK <- TRUE
PERMANOVA_N_PERM <- 99
PERMANOVA_MAX_STOCHASTIC_ITER <- 10
PERMANOVA_MAX_NUMERIC_DRIVERS <- 3

# Scénarios d'erreur entre espèces congénériques
GENUS_CONFUSION_RATES <- c(0.03, 0.05, 0.10)
RARE_WEIGHTED_MEAN_ERROR <- c(0.03, 0.05, 0.10)
RARE_WEIGHTED_MAX_ERROR <- 0.25

# Définition simple des espèces rares
RARE_MAX_TOTAL_ABUNDANCE <- 3

# TAXREF mainland : garder les espèces présentes/endémiques/cryptogènes en France métropolitaine.
# D = douteux exclu par défaut.
MAINLAND_FR_STATUS_KEEP <- c("P", "E", "C")

# Drivers environnementaux a priori pour le papier.
# On évite ici la sélection automatique de dizaines de variables :
# le test H2 porte explicitement sur trois gradients écologiquement interprétables.
FOCAL_ENV_DRIVERS <- c("t360_mean", "mos", "p_h")
DRIVER_CANDIDATES_NUMERIC <- FOCAL_ENV_DRIVERS
DRIVER_CANDIDATES_FACTOR <- character(0)
MAX_MISSING_PROP_DRIVER <- 0.25
MAX_ABS_COR_DRIVER <- 0.95  # utilisé seulement pour l'audit, pas pour exclure les trois drivers a priori

# divent
DIVENT_HILL_ESTIMATOR <- "naive"
DIVENT_COVERAGE_ESTIMATORS <- c("ZhangHuang", "Chao", "Turing")
COVERAGE_STANDARDIZATION_LEVEL <- NULL

# -----------------------------
# 1. OUTILS
# -----------------------------

message_header <- function(x) {
  cat("\n", strrep("=", 72), "\n", x, "\n", strrep("=", 72), "\n", sep = "")
}

safe_read_semicolon <- function(path, encoding = "Latin1") {
  readr::read_delim(
    path,
    delim = ";",
    locale = readr::locale(encoding = encoding),
    show_col_types = FALSE,
    name_repair = "unique"
  )
}

extract_binomial <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  out <- stringr::str_extract(
    x,
    "^[A-Z][A-Za-zÀ-ÖØ-öø-ÿ\\-]+\\s+[a-z][A-Za-zÀ-ÖØ-öø-ÿ\\-]+"
  )
  dplyr::na_if(out, "")
}

clean_chr <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[1]
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

safe_spearman <- function(x, y) {
  safe_cor(x, y, method = "spearman")
}

safe_top_driver <- function(term, r2) {
  ok <- !is.na(r2)
  if (!any(ok)) return(NA_character_)
  as.character(term[ok][which.max(r2[ok])])
}

safe_divent_hill <- function(x, q, level = NULL, estimator = DIVENT_HILL_ESTIMATOR) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0 || sum(x) <= 0) return(NA_real_)
  out <- tryCatch(
    divent::div_hill(x, q = q, estimator = estimator, level = level, as_numeric = TRUE),
    error = function(e) NA_real_
  )
  as.numeric(out[1])
}

safe_divent_coverage <- function(x, estimator = "Chao") {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0 || sum(x) <= 0) return(NA_real_)
  out <- tryCatch(
    divent::coverage(x, estimator = estimator, as_numeric = TRUE),
    error = function(e) NA_real_
  )
  as.numeric(out[1])
}

make_comm_matrix <- function(dat, unit_col = ANALYSIS_UNIT) {
  dat %>%
    filter(!is.na(.data[[unit_col]]), !is.na(taxon_unit), !is.na(abundance), abundance > 0) %>%
    group_by(unit = .data[[unit_col]], taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = taxon_unit, values_from = abundance, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("unit") %>%
    as.matrix()
}

hill_metrics_from_matrix <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) == 0) return(tibble())
  
  purrr::map_dfr(seq_len(nrow(mat)), function(i) {
    x <- as.numeric(mat[i, ])
    x_pos <- x[is.finite(x) & x > 0]
    
    coverage_vals <- purrr::map_dbl(
      DIVENT_COVERAGE_ESTIMATORS,
      ~ safe_divent_coverage(x_pos, estimator = .x)
    )
    names(coverage_vals) <- paste0("coverage_", tolower(DIVENT_COVERAGE_ESTIMATORS))
    
    out <- tibble(
      unit = rownames(mat)[i],
      total_abundance = sum(x_pos),
      q0 = safe_divent_hill(x_pos, q = 0),
      q1 = safe_divent_hill(x_pos, q = 1),
      q2 = safe_divent_hill(x_pos, q = 2),
      coverage_zhanghuang = as.numeric(coverage_vals["coverage_zhanghuang"]),
      coverage_chao = as.numeric(coverage_vals["coverage_chao"]),
      coverage_turing = as.numeric(coverage_vals["coverage_turing"])
    )
    
    if (!is.null(COVERAGE_STANDARDIZATION_LEVEL)) {
      out <- out %>% mutate(
        q0_covstd = safe_divent_hill(x_pos, q = 0, level = COVERAGE_STANDARDIZATION_LEVEL),
        q1_covstd = safe_divent_hill(x_pos, q = 1, level = COVERAGE_STANDARDIZATION_LEVEL),
        q2_covstd = safe_divent_hill(x_pos, q = 2, level = COVERAGE_STANDARDIZATION_LEVEL)
      )
    }
    out
  })
}

community_distance_vector <- function(mat, method = "bray", binary = FALSE) {
  mat <- as.matrix(mat)
  if (nrow(mat) < 3 || ncol(mat) < 2) return(NA_real_)
  keep_rows <- rowSums(mat) > 0
  keep_cols <- colSums(mat) > 0
  mat <- mat[keep_rows, keep_cols, drop = FALSE]
  if (nrow(mat) < 3 || ncol(mat) < 2) return(NA_real_)
  as.vector(vegan::vegdist(mat, method = method, binary = binary))
}

bray_distance_vector <- function(mat) {
  community_distance_vector(mat, method = "bray", binary = FALSE)
}

# Présence/absence : important pour tester la sensibilité aux espèces rares.
jaccard_distance_vector <- function(mat) {
  community_distance_vector(mat, method = "jaccard", binary = TRUE)
}

# vegan::vegdist(method = "bray", binary = TRUE) correspond à l'équivalent Sørensen/Bray en présence-absence.
sorensen_distance_vector <- function(mat) {
  community_distance_vector(mat, method = "bray", binary = TRUE)
}

safe_mean_ratio <- function(x_scenario, x_base) {
  if (all(is.na(x_scenario)) || all(is.na(x_base))) return(NA_real_)
  mb <- mean(x_base, na.rm = TRUE)
  ms <- mean(x_scenario, na.rm = TRUE)
  if (!is.finite(mb) || mb <= 0) return(NA_real_)
  ms / mb
}

pcoa_scores_2d <- function(mat) {
  mat <- as.matrix(mat)
  keep_rows <- rowSums(mat) > 0
  keep_cols <- colSums(mat) > 0
  mat <- mat[keep_rows, keep_cols, drop = FALSE]
  if (nrow(mat) < 4 || ncol(mat) < 2) return(NULL)
  d <- vegan::vegdist(mat, method = "bray")
  fit <- tryCatch(cmdscale(d, k = 2, eig = TRUE), error = function(e) NULL)
  if (is.null(fit) || is.null(fit$points) || ncol(fit$points) < 2) return(NULL)
  out <- as.data.frame(fit$points)
  names(out) <- c("PCoA1", "PCoA2")
  out$unit <- rownames(out)
  out
}

procrustes_r2 <- function(mat_base, mat_scen) {
  common <- intersect(rownames(mat_base), rownames(mat_scen))
  if (length(common) < 4) return(NA_real_)
  mb <- mat_base[common, , drop = FALSE]
  ms <- mat_scen[common, , drop = FALSE]
  sb <- pcoa_scores_2d(mb)
  ss <- pcoa_scores_2d(ms)
  if (is.null(sb) || is.null(ss)) return(NA_real_)
  common2 <- intersect(sb$unit, ss$unit)
  if (length(common2) < 4) return(NA_real_)
  sbm <- sb %>% filter(unit %in% common2) %>% arrange(match(unit, common2)) %>% select(PCoA1, PCoA2) %>% as.matrix()
  ssm <- ss %>% filter(unit %in% common2) %>% arrange(match(unit, common2)) %>% select(PCoA1, PCoA2) %>% as.matrix()
  fit <- tryCatch(vegan::procrustes(sbm, ssm, symmetric = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  as.numeric(1 - fit$ss)
}

stability_one <- function(mat_base, mat_scen) {
  common <- intersect(rownames(mat_base), rownames(mat_scen))
  if (length(common) < 3) {
    return(tibble(
      abundance_cor = NA_real_,
      q0_cor = NA_real_, q1_cor = NA_real_, q2_cor = NA_real_,
      coverage_chao_cor = NA_real_,
      bray_curtis_cor = NA_real_,
      jaccard_pa_cor = NA_real_,
      sorensen_pa_cor = NA_real_,
      bray_mean_ratio = NA_real_,
      jaccard_mean_ratio = NA_real_,
      sorensen_mean_ratio = NA_real_,
      procrustes_r2 = NA_real_,
      n_common_units = length(common)
    ))
  }
  
  mb <- mat_base[common, , drop = FALSE]
  ms <- mat_scen[common, , drop = FALSE]
  
  alpha_b <- hill_metrics_from_matrix(mb)
  alpha_s <- hill_metrics_from_matrix(ms)
  alpha_join <- alpha_b %>%
    select(unit, total_abundance, q0, q1, q2, coverage_chao) %>%
    rename_with(~ paste0(.x, "_base"), -unit) %>%
    inner_join(
      alpha_s %>%
        select(unit, total_abundance, q0, q1, q2, coverage_chao) %>%
        rename_with(~ paste0(.x, "_scenario"), -unit),
      by = "unit"
    )
  
  db <- bray_distance_vector(mb)
  ds <- bray_distance_vector(ms)
  jb <- jaccard_distance_vector(mb)
  js <- jaccard_distance_vector(ms)
  sob <- sorensen_distance_vector(mb)
  sos <- sorensen_distance_vector(ms)
  
  bray_cor <- if (length(db) == length(ds) && length(db) > 3) safe_cor(db, ds) else NA_real_
  jaccard_cor <- if (length(jb) == length(js) && length(jb) > 3) safe_cor(jb, js) else NA_real_
  sorensen_cor <- if (length(sob) == length(sos) && length(sob) > 3) safe_cor(sob, sos) else NA_real_
  
  tibble(
    abundance_cor = safe_cor(alpha_join$total_abundance_base, alpha_join$total_abundance_scenario),
    q0_cor = safe_cor(alpha_join$q0_base, alpha_join$q0_scenario),
    q1_cor = safe_cor(alpha_join$q1_base, alpha_join$q1_scenario),
    q2_cor = safe_cor(alpha_join$q2_base, alpha_join$q2_scenario),
    coverage_chao_cor = safe_cor(alpha_join$coverage_chao_base, alpha_join$coverage_chao_scenario),
    bray_curtis_cor = bray_cor,
    jaccard_pa_cor = jaccard_cor,
    sorensen_pa_cor = sorensen_cor,
    bray_mean_ratio = safe_mean_ratio(ds, db),
    jaccard_mean_ratio = safe_mean_ratio(js, jb),
    sorensen_mean_ratio = safe_mean_ratio(sos, sob),
    procrustes_r2 = procrustes_r2(mb, ms),
    n_common_units = length(common)
  )
}

select_uncorrelated <- function(dat, vars, max_abs_cor = 0.70) {
  vars <- vars[vars %in% names(dat)]
  selected <- character(0)
  for (v in vars) {
    if (length(selected) == 0) {
      selected <- c(selected, v)
    } else {
      cors <- sapply(selected, function(s) {
        suppressWarnings(cor(dat[[v]], dat[[s]], use = "pairwise.complete.obs"))
      })
      if (all(is.na(cors)) || max(abs(cors), na.rm = TRUE) < max_abs_cor) {
        selected <- c(selected, v)
      }
    }
  }
  selected
}

# -----------------------------
# 2. IMPORT ET NETTOYAGE
# -----------------------------

message_header("Import des donnees")

col_raw <- safe_read_semicolon(PATH_COLLEMBOLA, encoding = "Latin1") %>%
  janitor::clean_names(ascii = TRUE)

env_raw <- readr::read_csv(PATH_ENV, show_col_types = FALSE, name_repair = "unique") %>%
  janitor::clean_names(ascii = TRUE)

taxref_raw <- safe_read_semicolon(PATH_TAXREF, encoding = "UTF-8") %>%
  janitor::clean_names(ascii = TRUE)

# Alias robustes pour quelques variables dont le nettoyage des noms peut varier.

required_col <- c("station", "repetition", "cd_nom", "lb_nom", "nom_valide", "rank", "famille", "ordre", "abondance_totale")
missing_col <- setdiff(required_col, names(col_raw))
if (length(missing_col) > 0) {
  stop("Colonnes manquantes dans collembola.csv : ", paste(missing_col, collapse = ", "))
}

# TAXREF est utilisé comme référentiel d'identifiants :
# la base d'observation contient CD_NOM ; on récupère le CD_REF accepté via TAXREF.
# Les noms restent seulement des labels lisibles.
taxref_key_by_cdnom <- taxref_raw %>%
  mutate(
    cd_nom = clean_chr(cd_nom),
    cd_ref = clean_chr(cd_ref),
    rang = toupper(clean_chr(rang)),
    fr = clean_chr(fr),
    lb_nom = clean_chr(lb_nom),
    nom_valide = clean_chr(nom_valide),
    famille = clean_chr(famille),
    lb_binomial = extract_binomial(lb_nom),
    valid_binomial = extract_binomial(nom_valide)
  ) %>%
  transmute(
    cd_nom,
    taxref_cd_ref = cd_ref,
    taxref_rang = rang,
    taxref_fr = fr,
    taxref_lb_nom = lb_nom,
    taxref_nom_valide = nom_valide,
    taxref_lb_binomial = lb_binomial,
    taxref_valid_binomial = valid_binomial,
    taxref_famille = famille
  ) %>%
  distinct(cd_nom, .keep_all = TRUE)

col_dat <- col_raw %>%
  mutate(
    station = as.character(station),
    repetition = as.character(repetition),
    sample_id = paste(station, repetition, sep = "_"),
    abundance = as.numeric(abondance_totale),
    cd_nom = clean_chr(cd_nom),
    lb_nom = clean_chr(lb_nom),
    nom_valide = clean_chr(nom_valide),
    rank = toupper(clean_chr(rank)),
    famille = clean_chr(famille),
    ordre = clean_chr(ordre),
    binomial_lb = extract_binomial(lb_nom),
    binomial_valid = extract_binomial(nom_valide)
  ) %>%
  left_join(taxref_key_by_cdnom, by = "cd_nom") %>%
  mutate(
    # Binôme lisible : priorité au nom valide TAXREF quand le CD_NOM matche.
    binomial = if_else(
      rank %in% c("S", "ES"),
      coalesce(taxref_valid_binomial, binomial_lb, binomial_valid),
      NA_character_
    ),
    
    # Clé espèce robuste : CD_REF TAXREF si disponible, sinon fallback explicite.
    species_cd_ref = if_else(rank %in% c("S", "ES") & !is.na(taxref_cd_ref), taxref_cd_ref, NA_character_),
    species_key = case_when(
      rank %in% c("S", "ES") & !is.na(species_cd_ref) ~ paste0("cdref:", species_cd_ref),
      rank %in% c("S", "ES") & !is.na(binomial) ~ paste0("name_fallback:", binomial),
      TRUE ~ NA_character_
    ),
    species_label = if_else(rank %in% c("S", "ES"), coalesce(taxref_valid_binomial, binomial), NA_character_),
    
    genus = case_when(
      !is.na(species_label) ~ stringr::word(species_label, 1),
      rank %in% c("G", "GN") & !is.na(lb_nom) ~ stringr::word(lb_nom, 1),
      TRUE ~ NA_character_
    ),
    taxo_level = case_when(
      rank %in% c("S", "ES") & !is.na(species_key) ~ "species",
      rank %in% c("G", "GN") & !is.na(genus) ~ "genus",
      rank %in% c("F", "FM") & !is.na(famille) ~ "family",
      rank %in% c("O", "OR") & !is.na(ordre) ~ "order",
      TRUE ~ "unusable"
    ),
    
    # RTU = reported taxonomic unit. Les espèces sont identifiées par CD_REF quand possible.
    rtu_best = case_when(
      taxo_level == "species" ~ paste0("species:", species_key),
      taxo_level == "genus" ~ paste0("genus:", genus),
      taxo_level == "family" ~ paste0("family:", famille),
      TRUE ~ NA_character_
    ),
    genus_unit = if_else(taxo_level %in% c("species", "genus") & !is.na(genus), paste0("genus:", genus), NA_character_),
    family_unit = if_else(taxo_level %in% c("species", "genus", "family") & !is.na(famille), paste0("family:", famille), NA_character_),
    species_unit = if_else(taxo_level == "species" & !is.na(species_key), paste0("species:", species_key), NA_character_)
  ) %>%
  filter(!is.na(abundance), abundance > 0)

taxref_cdref_match_audit <- col_dat %>%
  mutate(
    taxref_match_status = case_when(
      taxo_level == "species" & !is.na(species_cd_ref) ~ "species_matched_to_cd_ref",
      taxo_level == "species" & is.na(species_cd_ref) ~ "species_unmatched_name_fallback",
      taxo_level %in% c("genus", "family", "order") & !is.na(taxref_cd_ref) ~ paste0(taxo_level, "_matched_to_cd_ref"),
      taxo_level %in% c("genus", "family", "order") & is.na(taxref_cd_ref) ~ paste0(taxo_level, "_unmatched"),
      TRUE ~ "unusable"
    )
  ) %>%
  count(taxo_level, rank, taxref_match_status, name = "n_rows") %>%
  arrange(taxo_level, rank, taxref_match_status)

readr::write_csv(taxref_cdref_match_audit, file.path(OUT_DIR, "taxref_cdref_match_audit.csv"))

taxref_unmatched_records <- col_dat %>%
  filter(is.na(taxref_cd_ref)) %>%
  distinct(cd_nom, lb_nom, nom_valide, rank, taxo_level, binomial, species_key, species_label) %>%
  arrange(rank, lb_nom)

readr::write_csv(taxref_unmatched_records, file.path(OUT_DIR, "taxref_unmatched_records.csv"))

# Audit : les individus seulement à l'ordre ne sont pas utilisés dans les matrices communautés.
taxo_audit <- col_dat %>%
  group_by(taxo_level) %>%
  summarise(
    n_rows = n(),
    total_abundance = sum(abundance),
    n_lb_nom = n_distinct(lb_nom, na.rm = TRUE),
    n_species_binomial = n_distinct(binomial, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(taxo_level, c("species", "genus", "family", "order", "unusable")))

readr::write_csv(taxo_audit, file.path(OUT_DIR, "taxonomic_resolution_audit.csv"))

excluded_too_coarse <- col_dat %>%
  filter(taxo_level %in% c("order", "unusable"))
readr::write_csv(excluded_too_coarse, file.path(OUT_DIR, "excluded_order_or_unusable_records.csv"))

col_used <- col_dat %>%
  filter(taxo_level %in% c("species", "genus", "family"))

message("Lignes totales collemboles : ", nrow(col_dat))
message("Lignes utilisées dans les analyses communauté : ", nrow(col_used))
message("Lignes exclues car ordre ou inutilisables : ", nrow(excluded_too_coarse))
message("Abondance exclue car ordre ou inutilisable : ", sum(excluded_too_coarse$abundance, na.rm = TRUE))
message("Stations : ", n_distinct(col_dat$station))
message("Echantillons station x repetition : ", n_distinct(col_dat$sample_id))

# -----------------------------
# 3. TAXREF MAINLAND + POOLS DE CONFUSION
# -----------------------------

message_header("Preparation TAXREF mainland et pools de confusion")

taxref_mainland <- taxref_raw %>%
  mutate(
    cd_nom = clean_chr(cd_nom),
    cd_ref = clean_chr(cd_ref),
    rang = toupper(clean_chr(rang)),
    fr = clean_chr(fr),
    lb_nom = clean_chr(lb_nom),
    nom_valide = clean_chr(nom_valide),
    famille = clean_chr(famille),
    lb_binomial = extract_binomial(lb_nom),
    valid_binomial = extract_binomial(nom_valide),
    # Pour les candidats TAXREF, on travaille sur le taxon valide accepté.
    species_cd_ref = cd_ref,
    species_key = if_else(!is.na(species_cd_ref), paste0("cdref:", species_cd_ref), NA_character_),
    species_label = coalesce(valid_binomial, lb_binomial),
    genus = stringr::word(species_label, 1)
  ) %>%
  filter(
    rang %in% c("ES", "S"),
    fr %in% MAINLAND_FR_STATUS_KEEP,
    !is.na(species_key),
    !is.na(species_label),
    !is.na(genus)
  ) %>%
  distinct(species_key, species_cd_ref, species_label, genus, famille, fr)

observed_species <- col_used %>%
  filter(taxo_level == "species") %>%
  group_by(species_key, species_cd_ref, species_label, binomial, genus, famille) %>%
  summarise(total_abundance = sum(abundance), n_rows = n(), .groups = "drop") %>%
  arrange(genus, species_label)

rare_species_keys <- observed_species %>%
  filter(total_abundance <= RARE_MAX_TOTAL_ABUNDANCE) %>%
  pull(species_key)

observed_same_genus_pool <- observed_species %>%
  transmute(
    candidate_key = species_key,
    candidate_label = species_label,
    genus,
    source_pool = "observed"
  ) %>%
  distinct()

taxref_same_genus_pool <- taxref_mainland %>%
  transmute(
    candidate_key = species_key,
    candidate_label = species_label,
    genus,
    source_pool = "taxref_mainland"
  ) %>%
  distinct()

confusion_pool_summary <- bind_rows(observed_same_genus_pool, taxref_same_genus_pool) %>%
  group_by(source_pool, genus) %>%
  summarise(n_candidate_species = n_distinct(candidate_key), .groups = "drop") %>%
  arrange(source_pool, genus)

readr::write_csv(observed_species, file.path(OUT_DIR, "observed_species_summary.csv"))
readr::write_csv(confusion_pool_summary, file.path(OUT_DIR, "confusion_pool_summary.csv"))

# Audit TAXREF des espèces observées : appariement prioritaire par CD_NOM -> CD_REF.
taxref_audit <- observed_species %>%
  mutate(
    in_taxref_by_cdref = !is.na(species_cd_ref),
    in_mainland_kept_status = species_key %in% taxref_mainland$species_key
  ) %>%
  left_join(
    taxref_mainland %>%
      transmute(
        species_key,
        taxref_mainland_label = species_label,
        taxref_direct_fr_status = fr,
        taxref_nom_valide = species_label
      ) %>%
      distinct(species_key, .keep_all = TRUE),
    by = "species_key"
  )

readr::write_csv(taxref_audit, file.path(OUT_DIR, "taxref_cdref_audit_observed_species.csv"))

message("Espèces observées au rang espèce : ", nrow(observed_species))
message("Espèces rares <= ", RARE_MAX_TOTAL_ABUNDANCE, " individus : ", length(rare_species_keys))
message("Genres observés avec au moins deux espèces observées : ",
        confusion_pool_summary %>% filter(source_pool == "observed", n_candidate_species >= 2) %>% nrow())

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

# -----------------------------
# 6. METRIQUES ALPHA + COUVERTURE
# -----------------------------

message_header("Calcul diversité alpha et couverture avec divent")

alpha_diversity_by_unit <- scenario_data %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
  group_split() %>%
  purrr::map_dfr(function(df) {
    mat <- make_comm_matrix(df, unit_col = ANALYSIS_UNIT)
    hill_metrics_from_matrix(mat) %>%
      mutate(
        scenario_family = df$scenario_family[1],
        scenario = df$scenario[1],
        baseline_scenario = df$baseline_scenario[1],
        unit_type = df$unit_type[1],
        iter = df$iter[1]
      )
  }) %>%
  relocate(scenario_family, scenario, baseline_scenario, unit_type, iter, unit)

readr::write_csv(alpha_diversity_by_unit, file.path(OUT_DIR, "alpha_diversity_and_coverage_by_unit.csv"))

gamma_by_scenario_iter <- scenario_data %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
  summarise(gamma_taxon_units = n_distinct(taxon_unit), .groups = "drop")

scenario_summary_by_iter <- alpha_diversity_by_unit %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
  summarise(
    n_units = n(),
    total_abundance = sum(total_abundance, na.rm = TRUE),
    mean_local_q0 = mean(q0, na.rm = TRUE),
    mean_local_q1 = mean(q1, na.rm = TRUE),
    mean_local_q2 = mean(q2, na.rm = TRUE),
    mean_coverage_chao = mean(coverage_chao, na.rm = TRUE),
    p10_coverage_chao = quantile(coverage_chao, probs = 0.10, na.rm = TRUE, names = FALSE),
    median_coverage_chao = median(coverage_chao, na.rm = TRUE),
    p90_coverage_chao = quantile(coverage_chao, probs = 0.90, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  left_join(gamma_by_scenario_iter, by = c("scenario_family", "scenario", "baseline_scenario", "unit_type", "iter"))

readr::write_csv(scenario_summary_by_iter, file.path(OUT_DIR, "scenario_summary_by_iter.csv"))

coverage_summary_by_iter <- alpha_diversity_by_unit %>%
  select(scenario_family, scenario, baseline_scenario, unit_type, iter, unit, starts_with("coverage_")) %>%
  pivot_longer(starts_with("coverage_"), names_to = "coverage_estimator", values_to = "coverage") %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter, coverage_estimator) %>%
  summarise(
    mean = mean(coverage, na.rm = TRUE),
    p10 = quantile(coverage, 0.10, na.rm = TRUE, names = FALSE),
    median = median(coverage, na.rm = TRUE),
    p90 = quantile(coverage, 0.90, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )
readr::write_csv(coverage_summary_by_iter, file.path(OUT_DIR, "coverage_summary_by_iter.csv"))

# -----------------------------
# 7. STABILITE VS BASELINE PROPRE A CHAQUE FAMILLE DE SCENARIOS
# -----------------------------

message_header("Calcul stabilité des inférences")

# Matrices par scénario/itération
matrix_index <- scenario_data %>%
  distinct(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
  arrange(scenario_family, scenario, iter)

make_matrix_for_key <- function(fam, scen, iter_i) {
  scenario_data %>%
    filter(scenario_family == fam, scenario == scen, iter == iter_i) %>%
    make_comm_matrix(unit_col = ANALYSIS_UNIT)
}

stability_by_iter <- matrix_index %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
  group_split() %>%
  purrr::map_dfr(function(key) {
    fam <- key$scenario_family[1]
    scen <- key$scenario[1]
    base <- key$baseline_scenario[1]
    iter_i <- key$iter[1]
    
    # Baseline déterministe : iter = 1.
    mat_base <- make_matrix_for_key(fam, base, 1)
    mat_scen <- make_matrix_for_key(fam, scen, iter_i)
    
    stability_one(mat_base, mat_scen) %>%
      mutate(
        scenario_family = fam,
        scenario = scen,
        baseline_scenario = base,
        unit_type = key$unit_type[1],
        iter = iter_i
      )
  }) %>%
  relocate(scenario_family, scenario, baseline_scenario, unit_type, iter)

readr::write_csv(stability_by_iter, file.path(OUT_DIR, "stability_by_iter.csv"))

stability_summary <- stability_by_iter %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type) %>%
  summarise(
    across(c(abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
             bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor,
             bray_mean_ratio, jaccard_mean_ratio, sorensen_mean_ratio,
             procrustes_r2, n_common_units),
           safe_mean, .names = "{.col}_mean"),
    across(c(abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
             bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor,
             bray_mean_ratio, jaccard_mean_ratio, sorensen_mean_ratio,
             procrustes_r2, n_common_units),
           safe_sd, .names = "{.col}_sd"),
    n_iter = n(),
    .groups = "drop"
  )
readr::write_csv(stability_summary, file.path(OUT_DIR, "stability_summary.csv"))

# Format long pour figures
stability_long <- stability_by_iter %>%
  select(scenario_family, scenario, baseline_scenario, unit_type, iter,
         abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
         bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2) %>%
  pivot_longer(
    cols = c(abundance_cor, q0_cor, q1_cor, q2_cor, coverage_chao_cor,
             bray_curtis_cor, jaccard_pa_cor, sorensen_pa_cor, procrustes_r2),
    names_to = "inference",
    values_to = "stability"
  ) %>%
  mutate(
    inference_label = recode(
      inference,
      abundance_cor = "total abundance",
      q0_cor = "q0 richness / number of units",
      q1_cor = "Hill q1 / effective units",
      q2_cor = "Hill q2 / dominant units",
      coverage_chao_cor = "sample coverage Chao",
      bray_curtis_cor = "Bray-Curtis abundance distance",
      jaccard_pa_cor = "Jaccard presence-absence distance",
      sorensen_pa_cor = "Sørensen presence-absence distance",
      procrustes_r2 = "ordination Procrustes R2"
    )
  )
readr::write_csv(stability_long, file.path(OUT_DIR, "stability_long.csv"))

# -----------------------------
# 8. DRIVERS ENVIRONNEMENTAUX : PERMANOVA UNIVARIEE
# -----------------------------

message_header("Bloc drivers environnementaux")

permanova_by_iter <- tibble()
permanova_summary <- tibble()
driver_rank_stability <- tibble()

if (RUN_DRIVER_BLOCK && length(selected_drivers) > 0) {
  
  meta_unit <- if (ANALYSIS_UNIT == "station") {
    meta_station %>% transmute(unit = station, across(all_of(selected_drivers)))
  } else {
    col_dat %>%
      distinct(sample_id, station) %>%
      left_join(meta_station, by = "station") %>%
      transmute(unit = sample_id, across(all_of(selected_drivers)))
  }
  
  num_to_scale <- intersect(selected_numeric_drivers_driverblock, names(meta_unit))
  meta_unit <- meta_unit %>%
    mutate(across(all_of(num_to_scale), ~ as.numeric(scale(.x))))
  
  run_univariate_permanova_one <- function(df_one, scenario_name, scenario_family_name, baseline_name, iter_id) {
    mat <- make_comm_matrix(df_one, unit_col = ANALYSIS_UNIT)
    common <- intersect(rownames(mat), meta_unit$unit)
    mat <- mat[common, , drop = FALSE]
    meta <- meta_unit %>% filter(unit %in% common) %>% arrange(match(unit, rownames(mat)))
    
    purrr::map_dfr(selected_drivers, function(term_i) {
      complete <- complete.cases(meta %>% select(all_of(term_i)))
      mat_i <- mat[complete, , drop = FALSE]
      meta_i <- meta[complete, , drop = FALSE]
      
      keep_cols <- if (ncol(mat_i) > 0) colSums(mat_i) > 0 else logical(0)
      mat_i <- mat_i[, keep_cols, drop = FALSE]
      
      if (
        nrow(mat_i) < 10 ||
        ncol(mat_i) < 2 ||
        dplyr::n_distinct(meta_i[[term_i]], na.rm = TRUE) < 2
      ) {
        return(tibble(
          scenario_family = scenario_family_name,
          scenario = scenario_name,
          baseline_scenario = baseline_name,
          iter = iter_id,
          term = term_i,
          df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
          n_units = nrow(mat_i), n_taxa = ncol(mat_i)
        ))
      }
      
      meta_tmp <- tibble(x = meta_i[[term_i]])
      fit <- tryCatch(
        vegan::adonis2(mat_i ~ x, data = meta_tmp, method = "bray", permutations = PERMANOVA_N_PERM),
        error = function(e) NULL
      )
      
      if (is.null(fit)) {
        return(tibble(
          scenario_family = scenario_family_name,
          scenario = scenario_name,
          baseline_scenario = baseline_name,
          iter = iter_id,
          term = term_i,
          df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
          n_units = nrow(mat_i), n_taxa = ncol(mat_i)
        ))
      }
      
      as.data.frame(fit) %>%
        rownames_to_column("adonis_term") %>%
        filter(adonis_term == "x") %>%
        janitor::clean_names() %>%
        transmute(
          scenario_family = scenario_family_name,
          scenario = scenario_name,
          baseline_scenario = baseline_name,
          iter = iter_id,
          term = term_i,
          df = df,
          sumsq = sum_of_sqs,
          r2 = r2,
          f = f,
          p = pr_f,
          n_units = nrow(mat_i),
          n_taxa = ncol(mat_i)
        )
    })
  }
  
  deterministic_scenarios <- scenario_data %>%
    group_by(scenario_family, scenario) %>%
    summarise(n_iter = n_distinct(iter), .groups = "drop") %>%
    filter(n_iter == 1) %>%
    pull(scenario)
  
  driver_scenario_data <- scenario_data %>%
    filter(
      scenario %in% deterministic_scenarios |
        iter <= PERMANOVA_MAX_STOCHASTIC_ITER
    )
  
  permanova_by_iter <- driver_scenario_data %>%
    distinct(scenario_family, scenario, baseline_scenario, iter) %>%
    group_by(scenario_family, scenario, baseline_scenario, iter) %>%
    group_split() %>%
    purrr::map_dfr(function(key) {
      df <- driver_scenario_data %>%
        filter(
          scenario_family == key$scenario_family[1],
          scenario == key$scenario[1],
          baseline_scenario == key$baseline_scenario[1],
          iter == key$iter[1]
        )
      run_univariate_permanova_one(
        df_one = df,
        scenario_name = key$scenario[1],
        scenario_family_name = key$scenario_family[1],
        baseline_name = key$baseline_scenario[1],
        iter_id = key$iter[1]
      )
    })
  
  permanova_summary <- permanova_by_iter %>%
    group_by(scenario_family, scenario, baseline_scenario, term) %>%
    summarise(
      r2_mean = safe_mean(r2),
      r2_sd = safe_sd(r2),
      p_median = safe_median(p),
      p_prop_lt_005 = if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE),
      n_iter = n(),
      n_valid = sum(!is.na(r2)),
      .groups = "drop"
    )
  
  baseline_rank <- permanova_summary %>%
    filter(scenario == baseline_scenario, !is.na(r2_mean)) %>%
    group_by(scenario_family, baseline_scenario) %>%
    mutate(rank_baseline = min_rank(desc(r2_mean))) %>%
    ungroup() %>%
    select(scenario_family, baseline_scenario, term, rank_baseline, r2_baseline = r2_mean)
  
  driver_rank_stability <- permanova_summary %>%
    group_by(scenario_family, scenario, baseline_scenario) %>%
    mutate(rank_scenario = min_rank(desc(r2_mean))) %>%
    ungroup() %>%
    left_join(baseline_rank, by = c("scenario_family", "baseline_scenario", "term")) %>%
    group_by(scenario_family, scenario, baseline_scenario) %>%
    summarise(
      n_complete_rank_pairs = sum(complete.cases(rank_baseline, rank_scenario)),
      driver_rank_spearman = safe_spearman(rank_baseline, rank_scenario),
      top_driver = safe_top_driver(term, r2_mean),
      baseline_top_driver = safe_top_driver(term, r2_baseline),
      top_driver_changed = case_when(
        is.na(top_driver) | is.na(baseline_top_driver) ~ NA,
        TRUE ~ top_driver != baseline_top_driver
      ),
      .groups = "drop"
    )
  
  # PERMANOVA multivariée marginale : T360_mean + MOS + pH dans le même modèle.
  run_multivariate_permanova_one <- function(df_one, scenario_name, scenario_family_name, baseline_name, iter_id) {
    mat <- make_comm_matrix(df_one, unit_col = ANALYSIS_UNIT)
    common <- intersect(rownames(mat), meta_unit$unit)
    mat <- mat[common, , drop = FALSE]
    meta <- meta_unit %>% filter(unit %in% common) %>% arrange(match(unit, rownames(mat)))
    
    complete <- complete.cases(meta %>% select(all_of(selected_drivers)))
    mat_i <- mat[complete, , drop = FALSE]
    meta_i <- meta[complete, , drop = FALSE]
    
    keep_cols <- if (ncol(mat_i) > 0) colSums(mat_i) > 0 else logical(0)
    mat_i <- mat_i[, keep_cols, drop = FALSE]
    
    if (nrow(mat_i) < 10 || ncol(mat_i) < 2 || length(selected_drivers) < 1) {
      return(tibble(
        scenario_family = scenario_family_name,
        scenario = scenario_name,
        baseline_scenario = baseline_name,
        iter = iter_id,
        term = selected_drivers,
        df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
        n_units = nrow(mat_i), n_taxa = ncol(mat_i), model = "multivariate_margin"
      ))
    }
    
    form <- as.formula(paste("mat_i ~", paste(selected_drivers, collapse = " + ")))
    fit <- tryCatch(
      vegan::adonis2(form, data = meta_i, method = "bray", permutations = PERMANOVA_N_PERM, by = "margin"),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      return(tibble(
        scenario_family = scenario_family_name,
        scenario = scenario_name,
        baseline_scenario = baseline_name,
        iter = iter_id,
        term = selected_drivers,
        df = NA_real_, sumsq = NA_real_, r2 = NA_real_, f = NA_real_, p = NA_real_,
        n_units = nrow(mat_i), n_taxa = ncol(mat_i), model = "multivariate_margin"
      ))
    }
    
    as.data.frame(fit) %>%
      rownames_to_column("adonis_term") %>%
      filter(adonis_term %in% selected_drivers) %>%
      janitor::clean_names() %>%
      transmute(
        scenario_family = scenario_family_name,
        scenario = scenario_name,
        baseline_scenario = baseline_name,
        iter = iter_id,
        term = adonis_term,
        df = df,
        sumsq = sum_of_sqs,
        r2 = r2,
        f = f,
        p = pr_f,
        n_units = nrow(mat_i),
        n_taxa = ncol(mat_i),
        model = "multivariate_margin"
      )
  }
  
  permanova_multivar_by_iter <- driver_scenario_data %>%
    distinct(scenario_family, scenario, baseline_scenario, iter) %>%
    group_by(scenario_family, scenario, baseline_scenario, iter) %>%
    group_split() %>%
    purrr::map_dfr(function(key) {
      df <- driver_scenario_data %>%
        filter(
          scenario_family == key$scenario_family[1],
          scenario == key$scenario[1],
          baseline_scenario == key$baseline_scenario[1],
          iter == key$iter[1]
        )
      run_multivariate_permanova_one(
        df_one = df,
        scenario_name = key$scenario[1],
        scenario_family_name = key$scenario_family[1],
        baseline_name = key$baseline_scenario[1],
        iter_id = key$iter[1]
      )
    })
  
  permanova_multivar_summary <- permanova_multivar_by_iter %>%
    group_by(scenario_family, scenario, baseline_scenario, term) %>%
    summarise(
      r2_mean = safe_mean(r2),
      r2_sd = safe_sd(r2),
      p_median = safe_median(p),
      p_prop_lt_005 = if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE),
      n_iter = n(),
      n_valid = sum(!is.na(r2)),
      .groups = "drop"
    )
  
  # Modèles alpha : stabilité des coefficients directionnels pour q0, q1, q2 et couverture.
  alpha_driver_input <- alpha_diversity_by_unit %>%
    left_join(meta_unit, by = "unit") %>%
    pivot_longer(
      cols = c(q0, q1, q2, coverage_chao),
      names_to = "alpha_metric",
      values_to = "alpha_value"
    )
  
  fit_alpha_driver_one <- function(df) {
    fam <- df$scenario_family[1]
    scen <- df$scenario[1]
    base <- df$baseline_scenario[1]
    iter_i <- df$iter[1]
    metric_i <- df$alpha_metric[1]
    
    dd <- df %>%
      select(alpha_value, all_of(selected_drivers)) %>%
      mutate(y = log1p(alpha_value)) %>%
      filter(complete.cases(.))
    
    if (nrow(dd) < 10 || n_distinct(dd$y) < 3 || length(selected_drivers) < 1) {
      return(tibble(
        scenario_family = fam, scenario = scen, baseline_scenario = base,
        iter = iter_i, alpha_metric = metric_i,
        term = selected_drivers, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p = NA_real_,
        adj_r2 = NA_real_, n_units = nrow(dd)
      ))
    }
    
    form <- as.formula(paste("y ~", paste(selected_drivers, collapse = " + ")))
    fit <- tryCatch(lm(form, data = dd), error = function(e) NULL)
    if (is.null(fit)) {
      return(tibble(
        scenario_family = fam, scenario = scen, baseline_scenario = base,
        iter = iter_i, alpha_metric = metric_i,
        term = selected_drivers, estimate = NA_real_, std_error = NA_real_, statistic = NA_real_, p = NA_real_,
        adj_r2 = NA_real_, n_units = nrow(dd)
      ))
    }
    
    adj_r2 <- broom::glance(fit)$adj.r.squared[1]
    broom::tidy(fit) %>%
      filter(term %in% selected_drivers) %>%
      transmute(
        scenario_family = fam, scenario = scen, baseline_scenario = base,
        iter = iter_i, alpha_metric = metric_i,
        term = term,
        estimate = estimate,
        std_error = std.error,
        statistic = statistic,
        p = p.value,
        adj_r2 = adj_r2,
        n_units = nrow(dd)
      )
  }
  
  alpha_driver_by_iter <- alpha_driver_input %>%
    filter(scenario %in% deterministic_scenarios | iter <= PERMANOVA_MAX_STOCHASTIC_ITER) %>%
    group_by(scenario_family, scenario, baseline_scenario, iter, alpha_metric) %>%
    group_split() %>%
    purrr::map_dfr(fit_alpha_driver_one)
  
  alpha_driver_summary <- alpha_driver_by_iter %>%
    group_by(scenario_family, scenario, baseline_scenario, alpha_metric, term) %>%
    summarise(
      estimate_mean = safe_mean(estimate),
      estimate_sd = safe_sd(estimate),
      p_median = safe_median(p),
      p_prop_lt_005 = if (all(is.na(p))) NA_real_ else mean(p < 0.05, na.rm = TRUE),
      adj_r2_mean = safe_mean(adj_r2),
      n_iter = n(),
      n_valid = sum(!is.na(estimate)),
      .groups = "drop"
    )
  
  alpha_driver_baseline <- alpha_driver_summary %>%
    filter(scenario == baseline_scenario) %>%
    transmute(
      scenario_family, baseline_scenario, alpha_metric, term,
      baseline_estimate = estimate_mean,
      baseline_abs_estimate_rank = min_rank(desc(abs(estimate_mean)))
    )
  
  alpha_driver_stability <- alpha_driver_summary %>%
    left_join(alpha_driver_baseline, by = c("scenario_family", "baseline_scenario", "alpha_metric", "term")) %>%
    mutate(
      estimate_sign_changed = case_when(
        is.na(estimate_mean) | is.na(baseline_estimate) ~ NA,
        estimate_mean == 0 | baseline_estimate == 0 ~ FALSE,
        TRUE ~ sign(estimate_mean) != sign(baseline_estimate)
      ),
      estimate_ratio = if_else(!is.na(baseline_estimate) & baseline_estimate != 0, estimate_mean / baseline_estimate, NA_real_)
    )
  
  readr::write_csv(permanova_multivar_by_iter, file.path(OUT_DIR, "permanova_multivar_by_iter.csv"))
  readr::write_csv(permanova_multivar_summary, file.path(OUT_DIR, "permanova_multivar_summary.csv"))
  readr::write_csv(alpha_driver_by_iter, file.path(OUT_DIR, "alpha_driver_by_iter.csv"))
  readr::write_csv(alpha_driver_summary, file.path(OUT_DIR, "alpha_driver_summary.csv"))
  readr::write_csv(alpha_driver_stability, file.path(OUT_DIR, "alpha_driver_stability.csv"))
  
  readr::write_csv(permanova_by_iter, file.path(OUT_DIR, "permanova_by_iter.csv"))
  readr::write_csv(permanova_summary, file.path(OUT_DIR, "permanova_summary.csv"))
  readr::write_csv(driver_rank_stability, file.path(OUT_DIR, "driver_rank_stability.csv"))
  
  message("PERMANOVA univariée : ", nrow(permanova_by_iter), " lignes écrites.")
  message("PERMANOVA multivariée : ", nrow(permanova_multivar_by_iter), " lignes écrites.")
  message("Modèles alpha : ", nrow(alpha_driver_by_iter), " lignes écrites.")
  
} else {
  warning("Bloc drivers ignoré : aucun driver sélectionné ou RUN_DRIVER_BLOCK = FALSE.")
  readr::write_csv(permanova_by_iter, file.path(OUT_DIR, "permanova_by_iter.csv"))
  readr::write_csv(permanova_summary, file.path(OUT_DIR, "permanova_summary.csv"))
  readr::write_csv(driver_rank_stability, file.path(OUT_DIR, "driver_rank_stability.csv"))
  readr::write_csv(tibble(), file.path(OUT_DIR, "permanova_multivar_by_iter.csv"))
  readr::write_csv(tibble(), file.path(OUT_DIR, "permanova_multivar_summary.csv"))
  readr::write_csv(tibble(), file.path(OUT_DIR, "alpha_driver_by_iter.csv"))
  readr::write_csv(tibble(), file.path(OUT_DIR, "alpha_driver_summary.csv"))
  readr::write_csv(tibble(), file.path(OUT_DIR, "alpha_driver_stability.csv"))
}

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
# 10. FIGURES
# -----------------------------

message_header("Figures")

theme_set(theme_bw(base_size = 12))

# Ordre d'affichage par logique analytique.
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

# Figure 1 : audit résolution initiale
p_rank <- col_dat %>%
  group_by(taxo_level) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  mutate(taxo_level = factor(taxo_level, levels = c("species", "genus", "family", "order", "unusable"))) %>%
  ggplot(aes(x = taxo_level, y = total_abundance)) +
  geom_col() +
  labs(
    title = "Initial taxonomic resolution in the dataset",
    subtitle = "Order-level and unusable records are audited but excluded from community analyses",
    x = "Reported taxonomic level",
    y = "Total abundance"
  )
ggsave(file.path(OUT_DIR, "fig_01_taxonomic_resolution_audit.png"), p_rank, width = 9, height = 5, dpi = 300)

# Figure 2 : stabilité par famille de scénario
p_stability <- stability_long %>%
  group_by(scenario_family, scenario, baseline_scenario, inference_label) %>%
  summarise(stability = mean(stability, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    scenario = factor(scenario, levels = scenario_order),
    inference_label = factor(
      inference_label,
      levels = c(
        "total abundance",
        "q0 richness / number of units",
        "Hill q1 / effective units",
        "Hill q2 / dominant units",
        "sample coverage Chao",
        "Bray-Curtis abundance distance",
        "Jaccard presence-absence distance",
        "Sørensen presence-absence distance",
        "ordination Procrustes R2"
      )
    )
  ) %>%
  ggplot(aes(x = scenario, y = stability, group = inference_label, linetype = inference_label, shape = inference_label)) +
  geom_hline(yintercept = 1, linewidth = 0.2) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  facet_wrap(~ scenario_family, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Robustness of ecological inference to taxonomic uncertainty",
    subtitle = "Each scenario is compared to the appropriate baseline for its analytical family",
    x = "Scenario",
    y = "Stability vs appropriate baseline",
    linetype = "Inference",
    shape = "Inference"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "fig_02_stability_clean_taxounits.png"), p_stability, width = 14, height = 8, dpi = 300)

# Figure 3 : couverture, avec quantiles plutôt que minimum instable
p_cov <- coverage_summary_by_iter %>%
  filter(coverage_estimator == "coverage_chao") %>%
  group_by(scenario_family, scenario, unit_type) %>%
  summarise(
    mean = mean(mean, na.rm = TRUE),
    p10 = mean(p10, na.rm = TRUE),
    median = mean(median, na.rm = TRUE),
    p90 = mean(p90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(scenario = factor(scenario, levels = scenario_order)) %>%
  ggplot(aes(x = scenario, y = mean)) +
  geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.15) +
  geom_point() +
  facet_wrap(~ scenario_family, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Sample coverage under taxonomic uncertainty",
    subtitle = "Mean Chao coverage; error bars = mean 10th and 90th percentiles across units",
    x = "Scenario",
    y = "Sample coverage"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "fig_03_sample_coverage_clean.png"), p_cov, width = 13, height = 7, dpi = 300)

# Figure 4 : homogénéisation / synthèse de métriques d'état
p_state <- scenario_summary_by_iter %>%
  group_by(scenario_family, scenario, unit_type) %>%
  summarise(
    gamma_taxon_units = mean(gamma_taxon_units, na.rm = TRUE),
    mean_local_q0 = mean(mean_local_q0, na.rm = TRUE),
    mean_local_q1 = mean(mean_local_q1, na.rm = TRUE),
    mean_local_q2 = mean(mean_local_q2, na.rm = TRUE),
    mean_coverage_chao = mean(mean_coverage_chao, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(gamma_taxon_units, mean_local_q0, mean_local_q1, mean_local_q2, mean_coverage_chao),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_order),
    metric = recode(
      metric,
      gamma_taxon_units = "Gamma number of units",
      mean_local_q0 = "Mean local q0",
      mean_local_q1 = "Mean local Hill q1",
      mean_local_q2 = "Mean local Hill q2",
      mean_coverage_chao = "Mean sample coverage Chao"
    )
  ) %>%
  ggplot(aes(x = scenario, y = value)) +
  geom_point() +
  facet_grid(metric ~ scenario_family, scales = "free", space = "free_x") +
  labs(
    title = "Consequences of taxonomic filtering and coarsening",
    x = "Scenario",
    y = "Value"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "fig_04_state_metrics_clean.png"), p_state, width = 14, height = 10, dpi = 300)

# Figure 4b : homogénéisation explicite des distances communautaires
p_dist_ratio <- stability_by_iter %>%
  select(scenario_family, scenario, baseline_scenario, iter,
         bray_mean_ratio, jaccard_mean_ratio, sorensen_mean_ratio) %>%
  pivot_longer(
    cols = c(bray_mean_ratio, jaccard_mean_ratio, sorensen_mean_ratio),
    names_to = "distance_metric",
    values_to = "mean_distance_ratio"
  ) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_order),
    distance_metric = recode(
      distance_metric,
      bray_mean_ratio = "Bray-Curtis abundance",
      jaccard_mean_ratio = "Jaccard presence-absence",
      sorensen_mean_ratio = "Sørensen presence-absence"
    )
  ) %>%
  group_by(scenario_family, scenario, distance_metric) %>%
  summarise(
    mean_ratio = mean(mean_distance_ratio, na.rm = TRUE),
    p10 = quantile(mean_distance_ratio, 0.10, na.rm = TRUE, names = FALSE),
    p90 = quantile(mean_distance_ratio, 0.90, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = scenario, y = mean_ratio, shape = distance_metric, linetype = distance_metric, group = distance_metric)) +
  geom_hline(yintercept = 1, linewidth = 0.2) +
  geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.12, na.rm = TRUE) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  facet_wrap(~ scenario_family, scales = "free_x") +
  labs(
    title = "Taxonomic uncertainty and artificial homogenisation",
    subtitle = "Mean pairwise distance ratio vs the appropriate baseline; values < 1 indicate homogenisation",
    x = "Scenario",
    y = "Mean distance ratio",
    shape = "Distance",
    linetype = "Distance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "fig_04b_distance_homogenization_clean.png"), p_dist_ratio, width = 14, height = 8, dpi = 300)

# Figure 5 : drivers PERMANOVA
if (nrow(permanova_summary) > 0) {
  baseline_terms_to_show <- permanova_summary %>%
    filter(scenario == baseline_scenario, !is.na(r2_mean)) %>%
    group_by(scenario_family) %>%
    arrange(desc(r2_mean), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    distinct(scenario_family, term)
  
  if (nrow(baseline_terms_to_show) > 0) {
    p_drivers <- permanova_summary %>%
      semi_join(baseline_terms_to_show, by = c("scenario_family", "term")) %>%
      mutate(scenario = factor(scenario, levels = scenario_order)) %>%
      ggplot(aes(x = scenario, y = r2_mean, group = term, linetype = term, shape = term)) +
      geom_line(na.rm = TRUE) +
      geom_point(na.rm = TRUE) +
      facet_wrap(~ scenario_family, scales = "free_x") +
      labs(
        title = "Instability of environmental drivers under taxonomic uncertainty",
        subtitle = "Univariate marginal PERMANOVA R²; top baseline drivers shown per analytical family",
        x = "Scenario",
        y = "Marginal PERMANOVA R²",
        linetype = "Driver",
        shape = "Driver"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, "fig_05_driver_stability_permanova_clean.png"), p_drivers, width = 14, height = 8, dpi = 300)
  }
}

# Figure 6 : ranking des drivers
if (nrow(driver_rank_stability) > 0) {
  p_rank_driver <- driver_rank_stability %>%
    mutate(scenario = factor(scenario, levels = scenario_order)) %>%
    ggplot(aes(x = scenario, y = driver_rank_spearman)) +
    geom_hline(yintercept = 1, linewidth = 0.2) +
    geom_col(na.rm = TRUE) +
    facet_wrap(~ scenario_family, scales = "free_x") +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = "Stability of environmental-driver ranking",
      subtitle = "Spearman correlation of driver ranking vs the appropriate baseline",
      x = "Scenario",
      y = "Driver-rank stability"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_DIR, "fig_06_driver_rank_stability_clean.png"), p_rank_driver, width = 13, height = 7, dpi = 300)
}

# Figure 6b : stabilité des coefficients des drivers sur les métriques alpha
if (exists("alpha_driver_stability") && nrow(alpha_driver_stability) > 0) {
  p_alpha_drivers <- alpha_driver_stability %>%
    filter(alpha_metric %in% c("q0", "q1", "q2", "coverage_chao")) %>%
    mutate(
      scenario = factor(scenario, levels = scenario_order),
      alpha_metric = recode(
        alpha_metric,
        q0 = "q0 / number of units",
        q1 = "Hill q1",
        q2 = "Hill q2",
        coverage_chao = "Chao coverage"
      )
    ) %>%
    ggplot(aes(x = scenario, y = estimate_mean, shape = term, group = term, linetype = term)) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    facet_grid(alpha_metric ~ scenario_family, scales = "free_x") +
    labs(
      title = "Stability of alpha-diversity driver coefficients",
      subtitle = "Linear models on log1p(alpha metric); predictors are scaled T360_mean, MOS and pH",
      x = "Scenario",
      y = "Standardised coefficient",
      shape = "Driver",
      linetype = "Driver"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_DIR, "fig_06b_alpha_driver_coefficients_clean.png"), p_alpha_drivers, width = 15, height = 10, dpi = 300)
}

# Figure 7 : heatmap stabilité
p_heat <- stability_long %>%
  group_by(scenario_family, scenario, inference_label) %>%
  summarise(stability = mean(stability, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    scenario = factor(scenario, levels = scenario_order),
    inference_label = factor(
      inference_label,
      levels = rev(c(
        "total abundance",
        "q0 richness / number of units",
        "Hill q1 / effective units",
        "Hill q2 / dominant units",
        "sample coverage Chao",
        "Bray-Curtis abundance distance",
        "Jaccard presence-absence distance",
        "Sørensen presence-absence distance",
        "ordination Procrustes R2"
      ))
    )
  ) %>%
  ggplot(aes(x = scenario, y = inference_label, fill = stability)) +
  geom_tile(color = "white") +
  facet_wrap(~ scenario_family, scales = "free_x") +
  scale_fill_viridis_c(option = "D", limits = c(0, 1), na.value = "grey90") +
  labs(
    title = "Robustness map of ecological inference",
    subtitle = "Grey = not calculable or not comparable",
    x = "Scenario",
    y = "Inference type",
    fill = "Stability"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "fig_07_robustness_heatmap_clean.png"), p_heat, width = 15, height = 8, dpi = 300)

# Figure 8 : espèces rares
p_rare <- rare_contribution %>%
  pivot_longer(cols = c(abundance, n_rows, n_lb_nom, n_species_binomial), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = rare_status, y = value)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Contribution of rare species and non-species RTUs to the dataset",
    x = NULL,
    y = "Count"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(OUT_DIR, "fig_08_rare_species_contribution_clean.png"), p_rare, width = 12, height = 7, dpi = 300)

# -----------------------------
# 11. FIN
# -----------------------------

message_header("Terminé")
message("Dossier de sortie : ", OUT_DIR)
message("Objets disponibles notamment :")
message("  - taxo_audit")
message("  - excluded_too_coarse")
message("  - driver_selection_used")
message("  - scenario_definitions")
message("  - alpha_diversity_by_unit")
message("  - stability_summary")
message("  - permanova_summary")
message("  - permanova_multivar_summary")
message("  - alpha_driver_summary")
message("  - alpha_driver_stability")
message("  - driver_rank_stability")

