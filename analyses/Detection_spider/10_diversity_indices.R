# ============================================================
# 10_diversity_indices.R — Diversité via divent::accum_hill
#   - Conserve la segmentation Pitfall : Pitfall10/8/6/4/2
#   - Indices ^qD pour q = 0, 1, 2 au point "Sample" (observé)
# Entrées : dat0, methods_use, pit_reps, out_dir
# Dépend de : 01_config.R, 02_packages.R, 03_read_parse.R
# ============================================================


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


### Coverage
div_cov <- meta2 %>%
  mutate(N = rowSums(X_counts2, na.rm = TRUE),
         coverage = cov_vec) %>%
  left_join(res_all, by="row_id") %>%
  transmute(site_id, method,
            N, coverage,
            richness = q0, Hill_q1 = q1, Hill_q2 = q2)

hab <- read.csv("data/derived-data/liste_habitat.csv", sep = ";", h = T) %>%
  rename(id_station = "STATION") %>%
  mutate(id_station = as.numeric(id_station))

div_cov2 <- div_cov %>%
  mutate(
    cov_clip = pmin(pmax(coverage, 1e-4), 1 - 1e-4)
  ) %>%
  separate(site_id, c("pjt", "yr", "city", "id_station"), sep="_")  %>%
  mutate(id_station = as.numeric(id_station)) %>%
  left_join(hab)

write.csv(div_cov2, file.path(out_dir, "diversity/Hill_div_results.csv"))

message("=== FIN 10_diversity_indices (accum_hill + segmentation Pitfall) ===  Résultats : ", normalizePath(out_dir))

