# ============================================================
# 08_diversity_indices.R — Diversité via divent::accum_hill
#   - Conserve la segmentation Pitfall : Pitfall10/8/6/4/2
#   - Indices ^qD pour q = 0, 1, 2 au point "Sample" (observé)
#   - Option : correction grossière par p̂ (species × method_seg)
# Entrées : dat0, methods_use, pit_reps, out_dir
# Dépend de : 01_config.R, 02_packages.R, 03_read_parse.R
# ============================================================

# -- Sécurité packages (si 02_packages.R ne charge pas 'divent') -------------
if (!requireNamespace("divent", quietly = TRUE)) {
  stop("Le package 'divent' est requis. Ajoute-le à R/02_packages.R (librarian::shelf(divent)).")
}

# ---------- Segmentation des méthodes (site × method_seg × species → abund) --
make_pit_subset <- function(n){
  dat0 %>%
    dplyr::filter(method == "Pitfall", !is.na(trap), trap %in% 1:n) %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(method_seg = paste0("Pitfall", n))
}

abund_long_parts <- list()

if ("Pitfall" %in% methods_use) {
  pit_parts <- lapply(pit_reps, make_pit_subset)
  names(pit_parts) <- paste0("Pitfall", pit_reps)
  abund_long_parts <- c(abund_long_parts, pit_parts)
}
if ("GPD" %in% methods_use) {
  abund_long_parts$GPD <- dat0 %>%
    dplyr::filter(method == "GPD") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(method_seg = "GPD")
}
if ("DVAC" %in% methods_use) {
  abund_long_parts$DVAC <- dat0 %>%
    dplyr::filter(method == "DVAC") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(method_seg = "DVAC")
}
if ("TRI" %in% methods_use) {
  abund_long_parts$TM <- dat0 %>%
    dplyr::filter(method == "TRI") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(method_seg = "TM")
}

stopifnot(length(abund_long_parts) > 0)
abund_long_seg <- dplyr::bind_rows(abund_long_parts)

# Ordre des méthodes (affichage)
method_order <- c(
  if ("GPD" %in% methods_use) "GPD",
  if ("Pitfall" %in% methods_use) paste0("Pitfall", pit_reps),
  if ("DVAC" %in% methods_use) "DVAC",
  if ("TRI" %in% methods_use) "TM"
)
method_levels <- intersect(method_order, unique(abund_long_seg$method_seg))
stopifnot(length(method_levels) >= 1)

# ---------- Matrice wide (une ligne = site × method_seg) ---------------------
comm_abund <- abund_long_seg %>%
  dplyr::transmute(site_id,
                   method = factor(method_seg, levels = method_levels),
                   species, abund) %>%
  tidyr::pivot_wider(names_from = species, values_from = abund, values_fill = 0) %>%
  dplyr::arrange(method, site_id)

  # --- helpers numériques & diversité "Sample" ---
  to_numeric_counts <- function(df) {
    df %>%
      dplyr::mutate(dplyr::across(dplyr::everything(),
                                  ~ suppressWarnings(as.numeric(.)))) %>%
      tidyr::replace_na(list())
  }
row_coverage <- function(v) {
  n <- sum(v, na.rm = TRUE)
  if (n <= 0) return(0)
  f1 <- sum(v == 1, na.rm = TRUE)
  1 - f1 / n
}
row_sample_div <- function(v) {
  n <- sum(v, na.rm = TRUE)
  if (n <= 0) return(c(q0 = 0, q1 = 0, q2 = 0))
  p <- v[v > 0] / n
  q0 <- length(p)
  H  <- -sum(p * log(p))
  q1 <- exp(H)
  q2 <- 1 / sum(p^2)
  c(q0 = q0, q1 = q1, q2 = q2)
}
get_q_sample <- function(abd_obj, q) {
  divent::accum_hill(abd_obj, q = q, n_simulations = 0) %>%
    dplyr::filter(.data$estimator == "Sample") %>%
    dplyr::select(site, diversity) %>%
    dplyr::rename(!!paste0("q", q) := diversity)
}

# --- matrice de comptes, id des lignes utiles ---
X_counts <- comm_abund %>% dplyr::select(-site_id, -method) %>% to_numeric_counts() %>% as.data.frame()
keep_npos <- rowSums(X_counts, na.rm = TRUE) > 0
X_counts2 <- X_counts[keep_npos, , drop = FALSE]
meta2 <- comm_abund[keep_npos, c("site_id","method")]
meta2$row_id <- seq_len(nrow(meta2))  # index interne pour re-coller

# --- couverture par ligne & split ---
cov_vec <- apply(X_counts2, 1, row_coverage)
idx_cov_pos   <- which(cov_vec > 0)
idx_cov_zero  <- which(cov_vec == 0)

# --- 1) lignes à couverture > 0 : accum_hill (Sample) ---
res_pos <- NULL
if (length(idx_cov_pos) > 0) {
  abd_pos <- divent::as_abundances(X_counts2[idx_cov_pos, , drop = FALSE])
  meta_pos <- tibble::tibble(site = abd_pos$site,
                             row_id = meta2$row_id[idx_cov_pos])
  q0_pos <- get_q_sample(abd_pos, 0)
  q1_pos <- get_q_sample(abd_pos, 1)
  q2_pos <- get_q_sample(abd_pos, 2)
  res_pos <- meta_pos %>%
    dplyr::left_join(q0_pos, by = "site") %>%
    dplyr::left_join(q1_pos, by = "site") %>%
    dplyr::left_join(q2_pos, by = "site") %>%
    dplyr::select(row_id, q0, q1, q2)
}

# --- 2) lignes à couverture == 0 : formules directes ---
res_zero <- NULL
if (length(idx_cov_zero) > 0) {
  div_mat <- t(apply(X_counts2[idx_cov_zero, , drop = FALSE], 1, row_sample_div))
  res_zero <- tibble::tibble(
    row_id = meta2$row_id[idx_cov_zero],
    q0 = div_mat[, "q0"],
    q1 = div_mat[, "q1"],
    q2 = div_mat[, "q2"]
  )
}

# --- combine & restitue ---
res_all <- dplyr::bind_rows(res_pos, res_zero) %>% dplyr::arrange(row_id)
div_hill <- meta2 %>%
  dplyr::left_join(res_all, by = "row_id") %>%
  dplyr::transmute(
    site_id, method,
    richness = q0,          # ^0D
    Hill_q1  = q1,          # ^1D
    Hill_q2  = q2           # ^2D
  )

# ---------- Modèles mixtes (si ≥ 2 méthodes présentes) -----------------------
fit_if_possible <- function(formula, data, outfile){
  if (nlevels(droplevels(data$method)) < 2) {
    cat("⚠️ Moins de 2 niveaux de méthode — modèle non ajusté :", deparse(formula), "\n")
    return(invisible(NULL))
  }
  m <- lme4::lmer(formula, data = data)
  capture.output(summary(m), file = file.path(out_dir, outfile))
  invisible(m)
}

m_S  <- fit_if_possible(richness ~ method + (1|site_id), data = div_hill, outfile = "div_model_richness.txt")
m_q1 <- fit_if_possible(Hill_q1  ~ method + (1|site_id), data = div_hill, outfile = "div_model_Hill_q1.txt")
m_q2 <- fit_if_possible(Hill_q2  ~ method + (1|site_id), data = div_hill, outfile = "div_model_Hill_q2.txt")

# ---------- Figures (violins) ------------------------------------------------
div_hill <- div_hill %>% dplyr::mutate(method = factor(method, levels = method_levels))

pS <- ggplot2::ggplot(div_hill, ggplot2::aes(x = method, y = richness)) +
  ggplot2::geom_violin(fill = "grey92") +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08), alpha = 0.5) +
  ggplot2::labs(title = "Richesse spécifique par méthode (accum_hill, Sample)",
                x = NULL, y = expression({}^0*D)) +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(out_dir, "div_richness_by_method.png"), pS, width = 8, height = 5, dpi = 200)

pQ1 <- ggplot2::ggplot(div_hill, ggplot2::aes(x = method, y = Hill_q1)) +
  ggplot2::geom_violin(fill = "grey92") +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08), alpha = 0.5) +
  ggplot2::labs(title = "Diversité de Hill q=1 par méthode (accum_hill, Sample)",
                x = NULL, y = expression({}^1*D)) +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(out_dir, "div_hill_q1_by_method.png"), pQ1, width = 8, height = 5, dpi = 200)

pQ2 <- ggplot2::ggplot(div_hill, ggplot2::aes(x = method, y = Hill_q2)) +
  ggplot2::geom_violin(fill = "grey92") +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08), alpha = 0.5) +
  ggplot2::labs(title = "Diversité de Hill q=2 par méthode (accum_hill, Sample)",
                x = NULL, y = expression({}^2*D)) +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(out_dir, "div_hill_q2_by_method.png"), pQ2, width = 8, height = 5, dpi = 200)

# ---------- (Option) Correction grossière par p̂ (species × method_seg) ------
p_file <- file.path(out_dir, "p_hat_by_method_unmarked_full.csv")
if (file.exists(p_file)) {
  p_tab2 <- readr::read_csv(p_file, show_col_types = FALSE) %>%
    dplyr::transmute(species, method = factor(method, levels = method_levels), p_hat) %>%
    dplyr::filter(!is.na(method), !is.na(p_hat), p_hat > 0) %>%
    dplyr::mutate(p_hat = pmin(p_hat, 0.999))
  
  abund_corr <- abund_long_seg %>%
    dplyr::transmute(site_id, method = factor(method_seg, levels = method_levels), species, abund) %>%
    dplyr::left_join(p_tab2, by = c("species","method")) %>%
    dplyr::mutate(abund_corr = ifelse(!is.na(p_hat) & p_hat > 0, abund / p_hat, abund))
  
  comm_abund_corr <- abund_corr %>%
    dplyr::select(site_id, method, species, abund_corr) %>%
    tidyr::pivot_wider(names_from = species, values_from = abund_corr, values_fill = 0) %>%
    dplyr::arrange(method, site_id)
  
  X_counts_corr <- comm_abund_corr %>% dplyr::select(-site_id, -method)
  abd_corr <- divent::as_abundances(X_counts_corr)
  
  meta_corr <- tibble::tibble(
    site    = abd_corr$site,
    site_id = comm_abund_corr$site_id,
    method  = comm_abund_corr$method
  )
  
  q0c <- get_q_sample(abd_corr, 0)
  q1c <- get_q_sample(abd_corr, 1)
  q2c <- get_q_sample(abd_corr, 2)
  
  div_abund_corr <- meta_corr %>%
    dplyr::left_join(q0c, by = "site") %>%
    dplyr::left_join(q1c, by = "site") %>%
    dplyr::left_join(q2c, by = "site") %>%
    dplyr::transmute(
      site_id, method,
      richness = q0, Hill_q1 = q1, Hill_q2 = q2
    )
  
  fit_if_possible(richness ~ method + (1|site_id), data = div_abund_corr, outfile = "div_model_richness_corrected.txt")
  fit_if_possible(Hill_q1  ~ method + (1|site_id), data = div_abund_corr, outfile = "div_model_Hill_q1_corrected.txt")
  fit_if_possible(Hill_q2  ~ method + (1|site_id), data = div_abund_corr, outfile = "div_model_Hill_q2_corrected.txt")
  
  # (Option) Profils ^qD (médiane + IQR) via accum_hill sur une grille de q
  q_grid <- seq(0, 2, by = 0.05)
  prof_list <- lapply(q_grid, function(qq){
    acc <- divent::accum_hill(abd_corr, q = qq, n_simulations = 0) %>%
      dplyr::filter(estimator == "Sample") %>%
      dplyr::select(site, diversity) %>%
      dplyr::rename(hill = diversity)
    tibble::tibble(q = qq) %>%
      dplyr::bind_cols(acc)
  })
  prof_all <- dplyr::bind_rows(prof_list) %>%
    dplyr::left_join(meta_corr, by = "site")
  
  hill_profiles_corr <- prof_all %>%
    dplyr::group_by(method, q) %>%
    dplyr::summarise(
      hill_med = stats::median(hill, na.rm = TRUE),
      hill_lo  = stats::quantile(hill, 0.25, na.rm = TRUE),
      hill_hi  = stats::quantile(hill, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(method = factor(method, levels = method_levels))
  
  pHillC <- ggplot2::ggplot(hill_profiles_corr, ggplot2::aes(x = q, y = hill_med, color = method)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = hill_lo, ymax = hill_hi, fill = method), alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_x_continuous("q (ordre de Hill)", breaks = seq(0, 2, 0.5), limits = c(0, 2)) +
    ggplot2::labs(title = "Profils de diversité (Hill) — abondances corrigées par p̂ (accum_hill, Sample)",
                  y = expression({}^q*D)) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "div_hill_profiles_corrected.png"), pHillC, width = 9, height = 5, dpi = 200)
}

message("=== FIN 08_diversity_indices (accum_hill + segmentation Pitfall) ===  Résultats : ", normalizePath(out_dir))

