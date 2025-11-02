# ============================================================
# 09_composition_beta.R — Composition (segmentation Pitfall OK)
#   - Venn (≤4 méthodes)
#   - PCoA (Jaccard binaire) + PERMANOVA (strata = site)
#   - Turnover vs Nestedness (betapart)
# Entrées attendues : dat0, methods_use, pit_reps, out_dir
# Dépend de : 01_config.R, 02_packages.R, 03_read_parse.R
# ============================================================

# ---------- Construction des abondances segmentées ----------
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

# Présence/absence
pa_long_seg <- abund_long_seg %>%
  dplyr::mutate(det = as.integer(abund > 0))

# Ordre lisible des méthodes segmentées
method_order <- c(
  if ("GPD" %in% methods_use) "GPD",
  if ("Pitfall" %in% methods_use) paste0("Pitfall", pit_reps),
  if ("DVAC" %in% methods_use) "DVAC",
  if ("TRI" %in% methods_use) "TM"
)
method_levels <- intersect(method_order, unique(pa_long_seg$method_seg))

# ---------- Matrice site×espèces (chaque ligne = site × method_seg) ----------
comm_pa <- pa_long_seg %>%
  dplyr::transmute(
    site_id,
    method = factor(method_seg, levels = method_levels),
    species, det
  ) %>%
  tidyr::pivot_wider(names_from = species, values_from = det, values_fill = 0)

# ---------- Venn (≤4 méthodes segmentées) ----------
# Priorité d’affichage : GPD, Pitfall10→2, DVAC, TM
prioritized <- method_order
methods_for_venn <- intersect(prioritized, method_levels)
if (length(methods_for_venn) > 4) methods_for_venn <- methods_for_venn[1:4]

if (length(methods_for_venn) >= 2) {
  site_list <- sort(unique(dat0$site_id))[1:min(12, length(unique(dat0$site_id)))]
  for (sid in site_list) {
    sets <- lapply(methods_for_venn, function(m) {
      pa_long_seg %>%
        dplyr::filter(site_id == sid, method_seg == m, det == 1L) %>%
        dplyr::pull(species) %>% unique()
    })
    names(sets) <- methods_for_venn
    if (all(lengths(sets) == 0)) next
    p_venn <- ggvenn::ggvenn(sets, fill_alpha = .3, stroke_size = .5) +
      ggplot2::ggtitle(paste("Venn —", sid))
    ggplot2::ggsave(
      file.path(out_dir, paste0("venn_", gsub("[^A-Za-z0-9]+","_", sid), ".png")),
      p_venn, width = 5, height = 4, dpi = 200
    )
  }
}

# ---------- PCoA (Jaccard binaire) + PERMANOVA (strata = site) ----------
comm_pa_all <- comm_pa %>% dplyr::relocate(site_id, method)
X <- comm_pa_all %>% dplyr::select(-site_id, -method) %>% as.matrix()

if (nrow(X) >= 3 && ncol(X) >= 1) {
  d_jac <- vegan::vegdist(sqrt(X), method = "jaccard", binary = TRUE)
  pco <- stats::cmdscale(d_jac, k = 2, eig = TRUE)
  scores <- as.data.frame(pco$points); colnames(scores) <- c("PCoA1","PCoA2")
  scores$site_id <- comm_pa_all$site_id
  scores$method  <- comm_pa_all$method
  
  p_pcoa <- ggplot2::ggplot(scores, ggplot2::aes(PCoA1, PCoA2, color = method)) +
    ggplot2::geom_point(size = 2, alpha = .8) +
    ggplot2::stat_ellipse(type = "norm", level = .67, linewidth = .6, alpha = .4) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCoA (Jaccard, présence/absence) — coloré par méthode (segmentation Pitfall)")
  ggplot2::ggsave(file.path(out_dir, "PCoA_methods.png"), p_pcoa, width = 7, height = 5, dpi = 200)
  
  set.seed(123)
  adon <- vegan::adonis2(d_jac ~ method, data = comm_pa_all,
                         permutations = 9999, strata = comm_pa_all$site_id, by = "margin")
  capture.output(adon, file = file.path(out_dir, "PERMANOVA_methods.txt"))
}

# ---------- Turnover vs Nestedness (betapart) par site ----------
pair_turnover_nestedness_by_site <- function(sid, methods_vec) {
  df <- pa_long_seg %>%
    dplyr::filter(site_id == sid, method_seg %in% methods_vec) %>%
    dplyr::select(method_seg, species, det) %>%
    tidyr::pivot_wider(names_from = species, values_from = det, values_fill = 0)
  
  if (nrow(df) < 2) return(NULL)
  rownames(df) <- df$method_seg
  mat <- as.data.frame(df %>% dplyr::select(-method_seg))
  
  B  <- betapart::betapart.core(mat)
  bp <- betapart::beta.pair(B, index.family = "jaccard")
  
  mat_to_long <- function(M, comp) {
    as.data.frame(as.table(as.matrix(M))) %>%
      stats::setNames(c("m1","m2","value")) %>%
      dplyr::filter(as.character(m1) < as.character(m2)) %>%
      dplyr::mutate(component = comp, site_id = sid)
  }
  
  dplyr::bind_rows(
    mat_to_long(bp$beta.jtu, "turnover"),
    mat_to_long(bp$beta.jne, "nestedness"),
    mat_to_long(bp$beta.jac, "total")
  )
}

beta_pairs <- purrr::map_dfr(
  sort(unique(dat0$site_id)),
  ~ pair_turnover_nestedness_by_site(.x, method_levels)
)

if (nrow(beta_pairs) > 0) {
  beta_summary <- beta_pairs %>%
    dplyr::group_by(m1, m2, component) %>%
    dplyr::summarise(
      median = median(value, na.rm = TRUE),
      q25    = stats::quantile(value, .25, na.rm = TRUE),
      q75    = stats::quantile(value, .75, na.rm = TRUE),
      n_site = dplyr::n(), .groups = "drop"
    )
  
  readr::write_csv(beta_pairs,  file.path(out_dir, "beta_pairs_by_site_turnover_nestedness.csv"))
  readr::write_csv(beta_summary,file.path(out_dir, "beta_pairs_summary_turnover_nestedness.csv"))
  
  p_beta <- beta_summary %>%
    dplyr::filter(m2 == 6) %>%
    mutate(pair = paste(as.numeric(m1)*2, " traps vs GPD"))

    
  
  p_beta <- ggplot(p_beta, aes(x = pair, y = median, fill = component)) +
    geom_col() +
    geom_errorbar(ggplot2::aes(ymin = q25, ymax = q75)) +
    coord_flip() +
    labs(
      x = "Paire de méthodes (segmentées)",
      y = "β (Jaccard) — médiane [IQR]",
      title = "Partition de la β-diversité : turnover vs nestedness\n(entre méthodes segmentées, au sein des sites)"
    ) +
    facet_wrap(component~.)+
    theme_minimal()
  
  ggplot2::ggsave(file.path(out_dir, "beta_partition_methods.png"), p_beta, width = 9, height = 5, dpi = 200)
}

message("=== FIN 09_composition_beta (segmentation Pitfall OK) ===  Résultats : ", normalizePath(out_dir))
