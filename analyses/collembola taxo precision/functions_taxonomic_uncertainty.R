# ============================================================
# functions_taxonomic_uncertainty.R
# Fonctions partagées par les étapes du workflow
# ============================================================

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

