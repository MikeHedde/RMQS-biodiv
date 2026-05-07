# ============================================================
# PARAMÈTRES GLOBAUX
# ============================================================
in_path        <- "data/raw-data/1.faune/carabidae.csv"
out_dir_base   <- "figures/detectability/carabidae/sensitivity"
target_project <- "RMQS_2025"
taxo_level     <- c("ES", "S")
min_sites_vec  <- 5:10   # analyse de sensibilité
min_sites_ref  <- 7
require_contrast <- TRUE

# Méthodes analysées
methods_use <- c("Pitfall", "GPD")   # parmi "Pitfall","GPD","DVAC","TRI"

# Réplicats Pitfall
pit_reps <- c(10, 8, 6, 4, 2)

# Constantes effort
pit_diam_m         <- 0.05
pit_circ_m         <- pi * pit_diam_m
gpd_open_circ_m    <- pi * pit_diam_m
gpd_cross_m        <- 1.0
gpd_unit_intensity <- gpd_open_circ_m + gpd_cross_m
tm_unit            <- 0.0625

dir.create(out_dir_base, showWarnings = FALSE, recursive = TRUE)

# Sensitivity
# ============================================================
# SETUP DOSSIERS SENSITIVITY (min_sites 5..10)
# - Crée l’arborescence complète (comme le pipeline main)
# - Ne liste rien dans sensitivity (supposé vide)
# - Prêt à être appelé avant ta boucle de sensibilité
# ============================================================

# Sous-dossiers attendus 
subdirs <- c(
  "occ_model_output",
  "sp_protocol_favored",
  "DARMHA",
  "species_selection",
  "trait",
  "diversity"
)

# 0) Crée le dossier racine sensitivity
dir.create(out_dir_base, recursive = TRUE, showWarnings = FALSE)

# 1) Crée min_sites_X + sous-dossiers
for (k in min_sites_vec) {
  base_k <- file.path(out_dir_base, paste0("min_sites_", k))
  dir.create(base_k, recursive = TRUE, showWarnings = FALSE)
  
  for (sd in subdirs) {
    dir.create(file.path(base_k, sd), recursive = TRUE, showWarnings = FALSE)
  }
  
  # 2) (optionnel mais utile) README de traçabilité
  readme_path <- file.path(base_k, "README.txt")
  if (!file.exists(readme_path)) {
    writeLines(
      c(
        "Sensitivity run folder",
        paste0("min_sites = ", k),
        paste0("created   = ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "",
        "Expected contents:",
        paste0("- ", subdirs, collapse = "\n")
      ),
      con = readme_path
    )
  }
}

## Liste pour stocker les résultats
res_sens <- vector("list", length(min_sites_vec))
names(res_sens) <- paste0("k", min_sites_vec)

## Boucle principale
for (k in min_sites_vec) {
  
  k_min  <- k                 # <-- lock
  slot   <- paste0("k", k_min)
  
  min_sites <- k_min
  out_dir   <- file.path(out_dir_base, paste0("min_sites_", k_min))
  
  source("analyses/Occupancy/04_detection_effort.R")
  source("analyses/Occupancy/05_occupancy_fit.R")
  source("analyses/Occupancy/06_occupancy_figs.R")
  source("analyses/Occupancy/07_occupancy_stat.R")
  source("analyses/Occupancy/09_traits_detectability.R")
  
  res_sens[[slot]] <- list(
    k         = k_min,
    n_species = dplyr::n_distinct(tab_p$species),
    tab_p     = tab_p2,
    tab_beta  = tab_beta,
    tab_BL    = m_tr_int,
    traits_trends = as.data.frame(trends),
    traits_pairs  = as.data.frame(pairs(trends))
  )
}


# Analyses de stabilité à produire
## Stabilité du classement des protocoles
med_p <- lapply(res_sens, function(x)
  x$tab_p %>%
    group_by(method) %>%
    summarise(p_med = median(p_hat_clip, na.rm = TRUE))
)
med_ref <- med_p[[paste0("k", min_sites_ref)]]

stab_rank <- lapply(names(med_p), function(nm) {
  tibble(
    k = as.integer(sub("k", "", nm)),
    rho = cor(
      med_p[[nm]]$p_med,
      med_ref$p_med,
      method = "spearman"
    )
  )
}) %>% bind_rows()

##Stabilité de l’effet effort
eff_slope <- purrr::imap(res_sens, \(x, nm)
                         x$tab_beta %>%
                           dplyr::filter(component == "det", term == "eff_z") %>%
                           dplyr::transmute(estimate, se,
                                            lcl = estimate - 1.96*se,
                                            ucl = estimate + 1.96*se,
                                            k = as.integer(sub("k","", nm)))
)

eff_slope_df <- dplyr::bind_rows(eff_slope)

##Stabilité des effets de traits
library(broom.mixed)

trait_eff <- purrr::imap_dfr(res_sens, \(x, nm) {
  broom.mixed::tidy(x$tab_BL, effects = "fixed") %>%
    mutate(k = as.integer(sub("k","", nm)))
})

# Ex: ne garder que la pente "log_BL" (effet taille pour la méthode de référence Pitfall4)
trait_eff_logBL <- trait_eff %>%
  filter(term == "log_BL")

trait_eff_slopes <- trait_eff %>%
  filter(term %in% c("log_BL",
                     "methodPitfall2:log_BL",
                     "methodPitfall6:log_BL",
                     "methodPitfall8:log_BL",
                     "methodGPD:log_BL"))

trait_eff_method_slopes <- trait_eff_slopes %>%
  select(k, term, estimate, std.error) %>%
  tidyr::pivot_wider(names_from = term, values_from = c(estimate, std.error))


######Figures
p_A41 <- ggplot(stab_rank, aes(x = k, y = rho)) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "Minimum number of occupied sites (min_sites)",
    y = "Spearman rank correlation (ρ)",
    title = "Stability of protocol ranking across min_sites thresholds"
  ) +
  theme_minimal()
ggsave(file.path(out_dir_base, "Annex_A4_Stability_protocol_ranking.png"),
       p_A41, width = 8, height = 9, dpi = 300)

p_A42 <- ggplot(eff_slope_df, aes(x = estimate, y = factor(k))) +
  geom_point() +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Effect of sampling effort on detectability (β_eff_z)",
    y = "min_sites threshold",
    title = "Stability of the effort effect across minimum sites thresholds"
  ) +
  theme_minimal()
ggsave(file.path(out_dir_base, "Annex_A4_Stability_effort_effect.png"),
       p_A42, width = 8, height = 9, dpi = 300)

trait_logBL_stab <- trait_eff_logBL %>%
  mutate(
    lcl = estimate - 1.96 * std.error,
    ucl = estimate + 1.96 * std.error
  )

ggplot(trait_logBL_stab, aes(x = k, y = estimate)) +
  geom_line() +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  labs(
    x = "min_sites threshold",
    y = "Effect of body size on detectability (β_logBL)",
    title = "Stability of trait–detectability relationship"
  ) +
  theme_minimal()


# 1) garder les 3 termes d’intérêt
terms_A43 <- c("log_BL", "methodPitfall8:log_BL", "methodGPD:log_BL")

a4_df <- trait_eff %>%
  filter(term %in% terms_A4) %>%
  mutate(
    lcl = estimate - 1.96 * std.error,
    ucl = estimate + 1.96 * std.error,
    term_lab = dplyr::recode(
      term,
      "log_BL" = "Body size effect (reference method: Pitfall4)",
      "methodPitfall8:log_BL" = "Δ slope vs Pitfall4 (Pitfall8 × body size)",
      "methodGPD:log_BL" = "Δ slope vs Pitfall4 (GPD × body size)"
    )
  )

# 2) Figure annexe A4 : stabilité des pentes (coefficients)
p_A43 <- ggplot(a4_df, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ term_lab, ncol = 1, scales = "free_y") +
  labs(
    x = "Minimum number of occupied sites",
    y = "Coefficient on logit(detectability) scale (β)",
    title = "Sensitivity of trait–detectability slopes to minimum sites"
  ) +
  theme_minimal()
ggsave(file.path(out_dir_base, "Annex_A4_trait_method_interactions_sensitivity.png"),
       p_A43, width = 8, height = 9, dpi = 300)

