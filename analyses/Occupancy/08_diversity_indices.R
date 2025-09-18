# ============================================================
# 08_diversity_indices.R — Diversité (Hill q=0..2)
#   -> conserve la segmentation des Pitfalls : Pitfall10/8/6/4/2
#   -> option: correction grossière par p̂ (species × method_seg)
# Entrées attendues : dat0, methods_use, pit_reps, out_dir
# Dépend de : 01_config.R, 02_packages.R, 03_read_parse.R
# ============================================================

# ---------- Helpers Hill ----------
hill_D <- function(p, q){
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  if (abs(q) < 1e-8) return(length(p))
  if (abs(q - 1) < 1e-8) return(exp(-sum(p * log(p))))
  (sum(p^q))^(1/(1 - q))
}
calc_div <- function(comm_wide, is_abundance = TRUE){
  X <- comm_wide %>% select(-site_id, -method) %>% as.matrix()
  out <- comm_wide %>% select(site_id, method)
  richness <- apply(X > 0, 1, sum, na.rm = TRUE)
  if (is_abundance) {
    den <- rowSums(X); den[den == 0] <- 1
    P <- sweep(X, 1, den, "/")
    shannon <- apply(P, 1, function(p) -sum(ifelse(p > 0, p * log(p), 0)))
    simpson <- 1 / apply(P^2, 1, sum)
    out$richness <- richness
    out$shannon  <- shannon
    out$simpson  <- simpson
  } else {
    out$richness <- richness
  }
  out
}

# ---------- Segmentation des méthodes ----------
# Construit un long-tableau d'abondances "abund_long_seg" avec method_seg
make_pit_subset <- function(n){
  dat0 %>%
    filter(method == "Pitfall", !is.na(trap), trap %in% 1:n) %>%
    group_by(site_id, species) %>%
    summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    mutate(method_seg = paste0("Pitfall", n))
}

abund_long_parts <- list()

if ("Pitfall" %in% methods_use) {
  pit_parts <- lapply(pit_reps, make_pit_subset)
  names(pit_parts) <- paste0("Pitfall", pit_reps)
  abund_long_parts <- c(abund_long_parts, pit_parts)
}

if ("GPD" %in% methods_use) {
  gpd_abund <- dat0 %>%
    filter(method == "GPD") %>%
    group_by(site_id, species) %>%
    summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    mutate(method_seg = "GPD")
  abund_long_parts$GPD <- gpd_abund
}

if ("DVAC" %in% methods_use) {
  dvac_abund <- dat0 %>%
    filter(method == "DVAC") %>%
    group_by(site_id, species) %>%
    summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    mutate(method_seg = "DVAC")
  abund_long_parts$DVAC <- dvac_abund
}

if ("TRI" %in% methods_use) {
  tm_abund <- dat0 %>%
    filter(method == "TRI") %>%
    group_by(site_id, species) %>%
    summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    mutate(method_seg = "TM")
  abund_long_parts$TM <- tm_abund
}

stopifnot(length(abund_long_parts) > 0)
abund_long_seg <- dplyr::bind_rows(abund_long_parts)

# Présence/absence dérivée
pa_long_seg <- abund_long_seg %>%
  mutate(det = as.integer(abund > 0))

# Méthodes (avec ordre souhaité) présentes dans les données
method_order <- c(
  if ("GPD" %in% methods_use) "GPD",
  if ("Pitfall" %in% methods_use) paste0("Pitfall", pit_reps),
  if ("DVAC" %in% methods_use) "DVAC",
  if ("TRI" %in% methods_use) "TM"
)
method_levels <- intersect(method_order, unique(pa_long_seg$method_seg))

# ---------- Matrices site × espèce par méthode_seg ----------
make_comm_matrix <- function(df_long, value = c("abund","det")) {
  val <- match.arg(value)
  df_long %>%
    transmute(site_id, method = factor(method_seg, levels = method_levels),
              species, value = .data[[val]]) %>%
    tidyr::pivot_wider(names_from = species, values_from = value, values_fill = 0)
}

comm_abund <- make_comm_matrix(abund_long_seg, "abund")
comm_pa    <- make_comm_matrix(pa_long_seg,    "det")

# ---------- Diversité observée (abondances) ----------
div_hill <- calc_div(comm_abund, TRUE) %>%
  mutate(Hill_q1 = exp(shannon), Hill_q2 = simpson)

# ---------- Modèles mixtes par indice ----------
fit_if_possible <- function(formula, data, outfile){
  # Ne lance le modèle que s'il y a ≥ 2 niveaux de méthode
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

# ---------- Figures (observé) ----------
div_hill <- div_hill %>% mutate(method = factor(method, levels = method_levels))

pS <- ggplot(div_hill, aes(x = method, y = richness)) +
  geom_violin(fill = "grey92") +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  labs(title = "Richesse spécifique par méthode",
       x = NULL, y = "Richesse (S)") +
  theme_minimal()
ggsave(file.path(out_dir, "div_richness_by_method.png"), pS, width = 8, height = 5, dpi = 200)

pQ1 <- ggplot(div_hill, aes(x = method, y = Hill_q1)) +
  geom_violin(fill = "grey92") +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  labs(title = "Diversité de Hill q=1 (exp(Shannon)) par méthode",
       x = NULL, y = expression({}^1*D)) +
  theme_minimal()
ggsave(file.path(out_dir, "div_hill_q1_by_method.png"), pQ1, width = 8, height = 5, dpi = 200)

pQ2 <- ggplot(div_hill, aes(x = method, y = Hill_q2)) +
  geom_violin(fill = "grey92") +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.5) +
  labs(title = "Diversité de Hill q=2 (Inverse Simpson) par méthode",
       x = NULL, y = expression({}^2*D)) +
  theme_minimal()
ggsave(file.path(out_dir, "div_hill_q2_by_method.png"), pQ2, width = 8, height = 5, dpi = 200)

# ---------- (Option) Correction grossière par p̂ (species × method_seg) ----------
p_file <- file.path(out_dir, "p_hat_by_method_unmarked_full.csv")
if (file.exists(p_file)) {
  p_tab2 <- readr::read_csv(p_file, show_col_types = FALSE) %>%
    transmute(species, method = factor(method, levels = method_levels), p_hat) %>%
    filter(!is.na(method), !is.na(p_hat), p_hat > 0) %>%
    mutate(p_hat = pmin(p_hat, 0.999))
  
  abund_corr <- abund_long_seg %>%
    transmute(site_id, method = factor(method_seg, levels = method_levels), species, abund) %>%
    left_join(p_tab2, by = c("species","method")) %>%
    mutate(abund_corr = ifelse(!is.na(p_hat) & p_hat > 0, abund / p_hat, abund))
  
  comm_abund_corr <- abund_corr %>%
    select(site_id, method, species, abund_corr) %>%
    tidyr::pivot_wider(names_from = species, values_from = abund_corr, values_fill = 0)
  
  div_abund_corr <- calc_div(comm_abund_corr, TRUE) %>%
    mutate(Hill_q1 = exp(shannon), Hill_q2 = simpson)
  
  fit_if_possible(richness ~ method + (1|site_id), data = div_abund_corr, outfile = "div_model_richness_corrected.txt")
  fit_if_possible(Hill_q1  ~ method + (1|site_id), data = div_abund_corr, outfile = "div_model_Hill_q1_corrected.txt")
  fit_if_possible(Hill_q2  ~ method + (1|site_id), data = div_abund_corr, outfile = "div_model_Hill_q2_corrected.txt")
  
  # Profils de Hill corrigés (médiane + IQR) par méthode
  q_grid <- seq(0, 2, by = 0.05)
  hill_profiles_corr <- comm_abund_corr %>%
    group_by(site_id, method) %>%
    group_modify(~{
      abund_vec <- as.numeric(.x %>% select(-1:-2))
      den <- sum(abund_vec)
      if (den == 0) tibble(q = q_grid, hill = 0) else {
        pvec <- abund_vec / den
        tibble(q = q_grid, hill = sapply(q_grid, function(q) hill_D(pvec, q)))
      }
    }) %>%
    ungroup() %>%
    group_by(method, q) %>%
    summarise(hill_med = median(hill, na.rm = TRUE),
              hill_lo  = quantile(hill, 0.25, na.rm = TRUE),
              hill_hi  = quantile(hill, 0.75, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(method = factor(method, levels = method_levels))
  
  pHillC <- ggplot(hill_profiles_corr, aes(x = q, y = hill_med, color = method)) +
    geom_ribbon(aes(ymin = hill_lo, ymax = hill_hi, fill = method), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_x_continuous("q (ordre de Hill)", breaks = seq(0, 2, 0.5), limits = c(0, 2)) +
    labs(title = "Profils de diversité (Hill) — abondances corrigées par p̂", y = expression({}^q*D)) +
    theme_minimal()
  ggsave(file.path(out_dir, "div_hill_profiles_corrected.png"), pHillC, width = 9, height = 5, dpi = 200)
}

message("=== FIN 08_diversity_indices (segmentation Pitfall OK) ===  Résultats : ", normalizePath(out_dir))
