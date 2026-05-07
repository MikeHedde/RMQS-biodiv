# ============================================================
# 11_accumulation_curves.R — Species Accumulation (Hill q=0)
#   - S'appuie sur dat0/methods_use/pit_reps/out_dir (scripts 01→04)
#   - Conserve la segmentation Pitfall : Pitfall1-10/8/6/4/2
#   - Génére : courbes par station (facets) + résumé "Sample" par méthode
#   - Sauvegarde : CSV + PNG dans out_dir
# Dépendances : dplyr, tidyr, ggplot2, forcats, divent, ggdist
# ============================================================

# --- Sécurité : vérifier les objets requis du workflow ---
stopifnot(
  exists("dat0"), exists("methods_use"), exists("pit_reps"), exists("out_dir")
)

# Paquets (ajoutés localement si pas déjà chargés par 02_packages.R)
suppressPackageStartupMessages({
  if (!requireNamespace("librarian", quietly = TRUE)) install.packages("librarian")
  librarian::shelf(forcats, divent, ggdist)
})

# --- Contrôles & helpers ---
`%||%` <- function(x, y) if (is.null(x)) y else x

# Méthodes segmentées (ordre d'affichage harmonisé avec 08/09)
method_order_seg <- c(
  if ("GPD" %in% methods_use) "GPD",
  if ("Pitfall" %in% methods_use) paste0("Pitfall", pit_reps),
  if ("DVAC" %in% methods_use) "DVAC",
  if ("TRI" %in% methods_use)  "TM"
)

# Étiquettes jolies pour les figures (ex: "Pitfall10" -> "Pitfall 1-10")
pretty_label <- function(x){
  x <- as.character(x)
  x <- dplyr::case_when(
    x == "TM"              ~ "Tri manuel",
    grepl("^Pitfall\\d+$", x) ~ gsub("^Pitfall(\\d+)$", "Pitfall 1-\\1", x),
    TRUE ~ x
  )
  x
}

# Palette par défaut (seulem. pour niveaux présents)
method_colors_all <- c(
  "GPD"           = "black",
  "Pitfall10"     = "purple4",
  "Pitfall8"      = "mediumpurple4",
  "Pitfall6"      = "mediumpurple3",
  "Pitfall4"      = "mediumpurple2",
  "Pitfall2"      = "mediumpurple1",
  "DVAC"          = "chartreuse4",
  "TM"            = "aquamarine4"
)

# --- Construire l'abondance segmentée par site × méthode_seg × espèce ---
abund_long_parts <- list()

# Pitfall segmenté via la colonne 'trap' construite dans 03_read_parse.R
if ("Pitfall" %in% methods_use) {
  if (!"trap" %in% names(dat0)) {
    warning("Colonne 'trap' absente de dat0 : segmentation Pitfall impossible. On utilisera Pitfall global.")
    pit_global <- dat0 %>%
      dplyr::filter(method == "Pitfall") %>%
      dplyr::group_by(site_id, species) %>%
      dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(method_seg = "Pitfall10")  # fallback = 10
    abund_long_parts$Pitfall <- pit_global
  } else {
    # segmente selon pit_reps (ex: 10, 8, 6, 4, 2) en utilisant traps 1:n
    pit_parts <- lapply(pit_reps, function(n){
      dat0 %>%
        dplyr::filter(method == "Pitfall", !is.na(trap), trap %in% 1:n) %>%
        dplyr::group_by(site_id, species) %>%
        dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(method_seg = paste0("Pitfall", n))
    })
    names(pit_parts) <- paste0("Pitfall", pit_reps)
    abund_long_parts <- c(abund_long_parts, pit_parts)
  }
}

# GPD
if ("GPD" %in% methods_use) {
  abund_long_parts$GPD <- dat0 %>%
    dplyr::filter(method == "GPD") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(method_seg = "GPD")
}

# DVAC
if ("DVAC" %in% methods_use) {
  abund_long_parts$DVAC <- dat0 %>%
    dplyr::filter(method == "DVAC") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(method_seg = "DVAC")
}

# TRI
if ("TRI" %in% methods_use) {
  abund_long_parts$TM <- dat0 %>%
    dplyr::filter(method == "TRI") %>%
    dplyr::group_by(site_id, species) %>%
    dplyr::summarise(abund = sum(as.numeric(abund %||% 0), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(method_seg = "TM")
}

stopifnot(length(abund_long_parts) > 0)
abund_long_seg <- dplyr::bind_rows(abund_long_parts)

# Filtrer les méthodes réellement présentes
method_levels <- intersect(method_order_seg, unique(abund_long_seg$method_seg))
if (length(method_levels) < 1) stop("Aucune méthode segmentée détectée dans les données.")

# --- Matrice wide pour divent::as_abundances (une ligne = site×méthode_seg) ---
mymat0 <- abund_long_seg %>%
  dplyr::transmute(site_id,
                   method = factor(method_seg, levels = method_levels),
                   species, abund) %>%
  tidyr::pivot_wider(names_from = species, values_from = abund, values_fill = 0) %>%
  dplyr::arrange(method, site_id)

# Conserver les colonnes d'annotation, et conversion en matrice d'abondances
if (!all(c("site_id", "method") %in% names(mymat0))) stop("Construction de mymat0 échouée.")
X_counts <- mymat0[, !(names(mymat0) %in% c("site_id","method")), drop = FALSE]

# as_abundances : chaque ligne = 1 échantillon (site×méthode_seg), colonnes = espèces
abd <- divent::as_abundances(X_counts)

# (robuste) filtrage via la matrice de comptes
counts <- X_counts %>% as.data.frame()
keep <- rowSums(counts, na.rm = TRUE) > 0   # ou > 1 si vous souhaitez
counts_use <- counts[keep, , drop = FALSE]

abd_use <- divent::as_abundances(counts_use)

meta_use <- tibble::tibble(
  site    = abd_use$site,
  site_id = comm_abund$site_id[keep],
  method  = comm_abund$method[keep]
)

set.seed(123)
acc0 <- divent::accum_hill(abd_use, q = 0, n_simulations = 50)

acc <- dplyr::left_join(acc0, meta_use, by = "site")


# Harmoniser l'ordre et les libellés de méthode pour les figures
meth_levels_pretty <- pretty_label(method_levels)
acc$methode <- factor(acc$methode, levels = meth_levels_pretty)

# --------- 1) Courbes d'accumulation (hors "Sample") par station ----------
acc_curves <- acc %>%
  dplyr::filter(.data$estimator != "Sample")

# Palette adaptée aux niveaux présents
present_methods_raw <- as.character(method_levels)
present_colors <- method_colors_all[names(method_colors_all) %in% present_methods_raw]
names(present_colors) <- pretty_label(names(present_colors))  # clés = labels jolis

p_acc <- ggplot(acc_curves, aes(x = level, y = diversity, group = methode)) +
  geom_path(aes(colour = methode)) +
  labs(
    x = "Nb d'individus (log10)",
    y = "Nombre d'espèces (log10)",
    title = "Courbes d'accumulation d'espèces (Hill q=0)",
    subtitle = "Par station (méthodes segmentées pour Pitfall)"
  ) +
  facet_wrap(station ~ ., scales = "free") +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  scale_color_manual(values = present_colors) +
  # Largeur de trait par 'estimator' (si plusieurs estimateurs disponibles)
  scale_discrete_manual("linewidth", values = c("Rarefaction" = 1.5, "Extrapolation" = .8, "Asymptote" = .5),
                        na.translate = TRUE)

ggsave(file.path(out_dir, "accumulation_curves_by_station.png"),
       p_acc, width = 12, height = 8, dpi = 200)

# --------- 2) Résumé "Sample" (distribution par méthode) ----------
acc_sample <- acc %>%
  dplyr::filter(.data$estimator == "Sample")

p_sample <- ggplot(acc_sample, aes(x = methode, y = diversity, fill = methode)) +
  ggdist::stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA) +
  ggdist::stat_slab(linewidth = 2) +
  labs(
    x = NULL,
    y = "Richesse observée (q=0)",
    title = "Distribution des richesses observées par méthode (Sample)"
  ) +
  theme_bw() +
  guides(fill = "none") +
  scale_fill_manual(values = present_colors) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(out_dir, "accumulation_sample_summary.png"),
       p_sample, width = 10, height = 6, dpi = 200)

# --------- 3) Export des données ----------
readr::write_csv(acc, file.path(out_dir, "accumulation_hill_q0_all.csv"))

message("=== FIN 11_accumulation_curves.R ===  Résultats : ", normalizePath(out_dir))
