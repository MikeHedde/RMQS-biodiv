# ============================================================
# PUBLICATION-READY FIGURES
# Adapté au script actuel : stability_long, stability_by_iter,
# alpha_driver_by_iter, scenario_definitions
# ============================================================

library(tidyverse)
library(forcats)
library(patchwork)

PUB_DIR <- file.path(OUT_DIR, "figures_publication_ready")
dir.create(PUB_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 0. Recharge les objets nécessaires si absents de la session
# ------------------------------------------------------------

if (!exists("OUT_DIR")) {
  OUT_DIR <- "outputs_taxonomic_uncertainty_v9"
}

if (!exists("scenario_definitions")) {
  path_scen <- file.path(OUT_DIR, "scenario_definitions.csv")
  if (file.exists(path_scen)) {
    scenario_definitions <- readr::read_csv(path_scen, show_col_types = FALSE)
  } else if (exists("stability_by_iter")) {
    scenario_definitions <- stability_by_iter %>%
      distinct(scenario_family, scenario, baseline_scenario, unit_type)
  } else {
    stop("scenario_definitions introuvable : relance le script principal ou vérifie OUT_DIR.")
  }
}

if (!exists("stability_long")) {
  path_stab_long <- file.path(OUT_DIR, "stability_long.csv")
  if (file.exists(path_stab_long)) {
    stability_long <- readr::read_csv(path_stab_long, show_col_types = FALSE)
  } else {
    stop("stability_long introuvable : relance le script principal ou vérifie OUT_DIR.")
  }
}

if (!exists("stability_by_iter")) {
  path_stab_iter <- file.path(OUT_DIR, "stability_by_iter.csv")
  if (file.exists(path_stab_iter)) {
    stability_by_iter <- readr::read_csv(path_stab_iter, show_col_types = FALSE)
  } else {
    stop("stability_by_iter introuvable : relance le script principal ou vérifie OUT_DIR.")
  }
}

if (!exists("alpha_driver_by_iter")) {
  path_alpha <- file.path(OUT_DIR, "alpha_driver_by_iter.csv")
  if (file.exists(path_alpha)) {
    alpha_driver_by_iter <- readr::read_csv(path_alpha, show_col_types = FALSE)
  } else {
    warning("alpha_driver_by_iter introuvable : la figure drivers alpha ne pourra pas être produite.")
    alpha_driver_by_iter <- tibble()
  }
}

# Recrée aussi scenario_order si absent
if (!exists("scenario_order")) {
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
}

# ------------------------------------------------------------
# 1. Baselines à retirer des figures principales
# ------------------------------------------------------------

baseline_scenarios <- scenario_definitions %>%
  filter(scenario == baseline_scenario) %>%
  pull(scenario) %>%
  unique()

# ------------------------------------------------------------
# 2. Scénarios à garder dans les figures principales
# Le code ignore automatiquement les scénarios absents.
# ------------------------------------------------------------

scenario_key_main <- tribble(
  ~scenario, ~scenario_label, ~scenario_family_label, ~scenario_order,
  
  # Species-level error
  "species_same_genus_observed_10pct",
  "Diffuse congeneric\nerror 10%",
  "Species-level error",
  1,
  
  "species_same_genus_observed_rare_weighted_mean_10pct",
  "Rare-weighted\nerror 10%",
  "Species-level error",
  2,
  
  "species_same_genus_taxref_mainland_10pct",
  "TAXREF congeneric\nerror 10%",
  "Species-level error",
  3,
  
  "species_difficult_genera_observed_20pct",
  "Difficult genera\nerror 20%",
  "Species-level error",
  4,
  
  "species_difficult_genera_taxref_mainland_20pct",
  "Difficult genera\nTAXREF 20%",
  "Species-level error",
  5,
  
  "species_expert_complexes_observed_20pct",
  "Expert complexes\nerror 20%",
  "Species-level error",
  6,
  
  # Workflow decisions
  "rtu_expert_difficult_genera_to_genus",
  "Difficult genera\nreported as genus",
  "Workflow decision",
  7,
  
  "rtu_rare_species_to_genus",
  "Rare species\nreported as genus",
  "Workflow decision",
  8,
  
  "rtu_species_only_drop_unresolved",
  "Drop unresolved\nrecords",
  "Workflow decision",
  9,
  
  # Coarsening
  "rtu_genus_level",
  "Genus-level\nresolution",
  "Taxonomic coarsening",
  10,
  
  "rtu_family_level",
  "Family-level\nresolution",
  "Taxonomic coarsening",
  11
) %>%
  filter(scenario %in% scenario_definitions$scenario)

# Vérification utile
message("Scénarios retenus pour les figures principales :")
print(scenario_key_main)

# ------------------------------------------------------------
# 3. Labels propres des métriques
# ------------------------------------------------------------

metric_key <- tribble(
  ~inference_label, ~metric_label, ~metric_order,
  
  "q0 richness / number of units",
  "q0 richness",
  1,
  
  "Hill q1 / effective units",
  "Hill q1",
  2,
  
  "Hill q2 / dominant units",
  "Hill q2",
  3,
  
  "sample coverage Chao",
  "Sample coverage",
  4,
  
  "Bray-Curtis abundance distance",
  "Bray-Curtis",
  5,
  
  "Jaccard presence-absence distance",
  "Jaccard",
  6,
  
  "Sørensen presence-absence distance",
  "Sørensen",
  7,
  
  "ordination Procrustes R2",
  "Ordination",
  8
)

# Je retire volontairement "total abundance" des figures principales :
# c'est souvent trivialement stable et peu informatif.

# ============================================================
# FIGURE 1 — Taxonomic uncertainty structure in the dataset
# ============================================================

if (!exists("PUB_DIR")) {
  PUB_DIR <- file.path(OUT_DIR, "figures_publication_ready")
  dir.create(PUB_DIR, showWarnings = FALSE, recursive = TRUE)
}

if (!exists("taxo_audit")) {
  path_taxo <- file.path(OUT_DIR, "taxonomic_resolution_audit.csv")
  if (file.exists(path_taxo)) {
    taxo_audit <- readr::read_csv(path_taxo, show_col_types = FALSE)
  } else {
    stop("taxo_audit introuvable.")
  }
}

if (!exists("rare_contribution")) {
  path_rare <- file.path(OUT_DIR, "rare_species_contribution.csv")
  if (file.exists(path_rare)) {
    rare_contribution <- readr::read_csv(path_rare, show_col_types = FALSE)
  } else {
    stop("rare_contribution introuvable.")
  }
}

taxo_level_order <- c("species", "genus", "family", "order", "unusable")

p_fig1a <- taxo_audit %>%
  mutate(
    taxo_level = factor(taxo_level, levels = taxo_level_order)
  ) %>%
  select(taxo_level, n_rows, total_abundance) %>%
  pivot_longer(
    cols = c(n_rows, total_abundance),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  mutate(percent = 100 * value / sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    metric = recode(
      metric,
      n_rows = "Records",
      total_abundance = "Individuals"
    )
  ) %>%
  ggplot(aes(x = taxo_level, y = percent)) +
  geom_col(width = 0.75) +
  facet_wrap(~ metric, nrow = 1) +
  labs(
    title = "A. Initial taxonomic resolution",
    x = NULL,
    y = "Share of dataset (%)"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 11)
  )

p_fig1b <- rare_contribution %>%
  mutate(
    rare_status = recode(
      rare_status,
      rare_species = "Rare species",
      non_rare_species = "Non-rare species",
      non_species_RTU = "Non-species RTUs"
    ),
    rare_status = factor(
      rare_status,
      levels = c("Non-rare species", "Rare species", "Non-species RTUs")
    )
  ) %>%
  select(rare_status, abundance, n_rows, n_lb_nom) %>%
  pivot_longer(
    cols = c(abundance, n_rows, n_lb_nom),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  mutate(percent = 100 * value / sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    metric = recode(
      metric,
      abundance = "Individuals",
      n_rows = "Records",
      n_lb_nom = "Reported names"
    )
  ) %>%
  ggplot(aes(x = rare_status, y = percent)) +
  geom_col(width = 0.75) +
  facet_wrap(~ metric, nrow = 1) +
  labs(
    title = "B. Rare species and non-species RTUs",
    x = NULL,
    y = "Share of dataset (%)"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 11)
  )

# Optional panel C: expert targets if col_used exists in the session
expert_target_summary <- tibble()

if (exists("col_used")) {
  difficult_genera <- c("Folsomia", "Mesaphorura", "Isotoma")
  
  expert_complex_sets <- list(
    "Lepidocyrtus–Entomobrya complex" = c(
      "Lepidocyrtus lanuginosus",
      "Lepidocyrtus cyaneus",
      "Lepidocyrtus lignorum",
      "Entomobrya lanuginosa"
    ),
    "Ceratophysella complex" = c(
      "Ceratophysella bengtssoni",
      "Ceratophysella denticulata"
    ),
    "Sminthurides complex" = c(
      "Sminthurides malmgreni",
      "Sminthurides parvulus",
      "Sminthurides schoetti",
      "Sminthurides signatus"
    )
  )
  
  complex_lookup <- purrr::imap_dfr(
    expert_complex_sets,
    ~ tibble(species_label = .x, expert_group = .y)
  )
  
  expert_target_summary <- col_used %>%
    mutate(
      species_label_clean = coalesce(species_label, binomial),
      expert_group = case_when(
        genus %in% difficult_genera ~ paste0(genus, " difficult genus"),
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(complex_lookup, by = c("species_label_clean" = "species_label"), suffix = c("", "_complex")) %>%
    mutate(
      expert_group = coalesce(expert_group, expert_group_complex)
    ) %>%
    filter(!is.na(expert_group)) %>%
    group_by(expert_group) %>%
    summarise(
      abundance = sum(abundance, na.rm = TRUE),
      n_records = n(),
      n_stations = n_distinct(station),
      .groups = "drop"
    ) %>%
    mutate(
      percent_abundance = 100 * abundance / sum(col_used$abundance, na.rm = TRUE),
      expert_group = forcats::fct_reorder(expert_group, percent_abundance)
    )
}

if (nrow(expert_target_summary) > 0) {
  
  p_fig1c <- expert_target_summary %>%
    ggplot(aes(x = expert_group, y = percent_abundance)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(
      title = "C. Expert-informed uncertainty targets",
      x = NULL,
      y = "Share of total abundance (%)"
    ) +
    theme_classic(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold", size = 11)
    )
  
  p_fig1 <- (p_fig1a / p_fig1b / p_fig1c) +
    patchwork::plot_annotation(
      title = "Taxonomic structure of the Collembola dataset",
      subtitle = "Initial taxonomic resolution, rare taxa and expert-informed uncertainty targets",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11)
      )
    )
  
  ggsave(
    file.path(PUB_DIR, "Fig1_taxonomic_uncertainty_structure.png"),
    p_fig1,
    width = 9.5,
    height = 8.2,
    dpi = 600
  )
  
  ggsave(
    file.path(PUB_DIR, "Fig1_taxonomic_uncertainty_structure.pdf"),
    p_fig1,
    width = 9.5,
    height = 8.2
  )
  
} else {
  
  p_fig1 <- (p_fig1a / p_fig1b) +
    patchwork::plot_annotation(
      title = "Taxonomic structure of the Collembola dataset",
      subtitle = "Initial taxonomic resolution and contribution of rare species and non-species RTUs",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11)
      )
    )
  
  ggsave(
    file.path(PUB_DIR, "Fig1_taxonomic_uncertainty_structure.png"),
    p_fig1,
    width = 9.5,
    height = 6.2,
    dpi = 600
  )
  
  ggsave(
    file.path(PUB_DIR, "Fig1_taxonomic_uncertainty_structure.pdf"),
    p_fig1,
    width = 9.5,
    height = 6.2
  )
}

# ------------------------------------------------------------
# 4. FIGURE 2 — Main robustness heatmap
# ------------------------------------------------------------

robustness_pub <- stability_long %>%
  filter(!scenario %in% baseline_scenarios) %>%
  inner_join(scenario_key_main, by = "scenario") %>%
  inner_join(metric_key, by = "inference_label") %>%
  group_by(
    scenario_family,
    scenario,
    scenario_label,
    scenario_family_label,
    scenario_order,
    metric_label,
    metric_order
  ) %>%
  summarise(
    stability_mean = mean(stability, na.rm = TRUE),
    stability_p10 = quantile(stability, 0.10, na.rm = TRUE, names = FALSE),
    stability_p90 = quantile(stability, 0.90, na.rm = TRUE, names = FALSE),
    n_iter = n(),
    .groups = "drop"
  ) %>%
  mutate(
    stability_clipped = pmin(1, pmax(0, stability_mean)),
    loss_of_stability = 1 - stability_clipped,
    metric_label = factor(metric_label, levels = rev(metric_key$metric_label)),
    scenario_label = factor(
      scenario_label,
      levels = scenario_key_main %>%
        arrange(scenario_order) %>%
        pull(scenario_label) %>%
        unique()
    ),
    scenario_family_label = factor(
      scenario_family_label,
      levels = c("Species-level error", "Workflow decision", "Taxonomic coarsening")
    )
  )

p_robustness_main <- robustness_pub %>%
  ggplot(aes(x = scenario_label, y = metric_label, fill = loss_of_stability)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(. ~ scenario_family_label, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    na.value = "grey90",
    name = "Loss of\nstability"
  ) +
  labs(
    title = "Sensitivity of ecological inference to taxonomic uncertainty",
    subtitle = "Values are expressed relative to the appropriate baseline",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PUB_DIR, "Fig2_main_robustness_heatmap.png"),
  p_robustness_main,
  width = 9.5,
  height = 4.8,
  dpi = 600
)

ggsave(
  file.path(PUB_DIR, "Fig2_main_robustness_heatmap.pdf"),
  p_robustness_main,
  width = 9.5,
  height = 4.8
)

# ------------------------------------------------------------
# 5. FIGURE 3 — Artificial homogenisation
# ------------------------------------------------------------

homogenization_pub <- stability_by_iter %>%
  filter(!scenario %in% baseline_scenarios) %>%
  inner_join(scenario_key_main, by = "scenario") %>%
  select(
    scenario_family,
    scenario,
    scenario_label,
    scenario_family_label,
    scenario_order,
    iter,
    bray_mean_ratio,
    jaccard_mean_ratio,
    sorensen_mean_ratio
  ) %>%
  pivot_longer(
    cols = c(bray_mean_ratio, jaccard_mean_ratio, sorensen_mean_ratio),
    names_to = "distance_metric",
    values_to = "distance_ratio"
  ) %>%
  mutate(
    distance_label = recode(
      distance_metric,
      bray_mean_ratio = "Bray-Curtis",
      jaccard_mean_ratio = "Jaccard",
      sorensen_mean_ratio = "Sørensen"
    ),
    percent_change = 100 * (distance_ratio - 1)
  ) %>%
  group_by(
    scenario,
    scenario_label,
    scenario_family_label,
    scenario_order,
    distance_label
  ) %>%
  summarise(
    mean_change = mean(percent_change, na.rm = TRUE),
    p10 = quantile(percent_change, 0.10, na.rm = TRUE, names = FALSE),
    p90 = quantile(percent_change, 0.90, na.rm = TRUE, names = FALSE),
    n_iter = n(),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_key_main %>%
        arrange(scenario_order) %>%
        pull(scenario_label) %>%
        unique()
    ),
    scenario_family_label = factor(
      scenario_family_label,
      levels = c("Species-level error", "Workflow decision", "Taxonomic coarsening")
    ),
    distance_label = factor(distance_label, levels = c("Bray-Curtis", "Jaccard", "Sørensen"))
  )

p_homogenization_main <- homogenization_pub %>%
  ggplot(aes(x = scenario_label, y = mean_change, shape = distance_label)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed") +
  geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.12, na.rm = TRUE) +
  geom_point(size = 2.4, na.rm = TRUE) +
  facet_grid(. ~ scenario_family_label, scales = "free_x", space = "free_x") +
  labs(
    title = "Taxonomic uncertainty can artificially homogenise communities",
    subtitle = "Change in mean pairwise dissimilarity relative to the appropriate baseline",
    x = NULL,
    y = "Change in mean pairwise dissimilarity (%)",
    shape = "Dissimilarity"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PUB_DIR, "Fig3_main_homogenisation.png"),
  p_homogenization_main,
  width = 9.5,
  height = 4.2,
  dpi = 600
)

ggsave(
  file.path(PUB_DIR, "Fig3_main_homogenisation.pdf"),
  p_homogenization_main,
  width = 9.5,
  height = 4.2
)

# ------------------------------------------------------------
# 6. FIGURE 4 — Change in alpha-diversity driver coefficients
# ------------------------------------------------------------

driver_label_key <- c(
  "t360_mean" = "Temperature",
  "mos" = "Soil organic matter",
  "p_h" = "pH",
  "ph" = "pH"
)

alpha_metric_key <- tribble(
  ~alpha_metric, ~alpha_metric_label, ~alpha_metric_order,
  "q0", "q0 richness", 1,
  "q1", "Hill q1", 2,
  "q2", "Hill q2", 3,
  "coverage_chao", "Sample coverage", 4
)

alpha_driver_delta_pub <- alpha_driver_by_iter %>%
  filter(alpha_metric %in% alpha_metric_key$alpha_metric) %>%
  mutate(
    driver_label = recode(term, !!!driver_label_key),
    driver_label = if_else(is.na(driver_label), term, driver_label)
  )

baseline_alpha_coeff <- alpha_driver_delta_pub %>%
  filter(scenario == baseline_scenario) %>%
  select(
    scenario_family,
    baseline_scenario,
    alpha_metric,
    term,
    baseline_estimate = estimate
  )

alpha_driver_delta_pub <- alpha_driver_delta_pub %>%
  filter(!scenario %in% baseline_scenarios) %>%
  inner_join(
    baseline_alpha_coeff,
    by = c("scenario_family", "baseline_scenario", "alpha_metric", "term")
  ) %>%
  inner_join(scenario_key_main, by = "scenario") %>%
  inner_join(alpha_metric_key, by = "alpha_metric") %>%
  mutate(
    delta_estimate = estimate - baseline_estimate,
    scenario_label = factor(
      scenario_label,
      levels = scenario_key_main %>%
        arrange(scenario_order) %>%
        pull(scenario_label) %>%
        unique()
    ),
    scenario_family_label = factor(
      scenario_family_label,
      levels = c("Species-level error", "Workflow decision", "Taxonomic coarsening")
    ),
    alpha_metric_label = factor(
      alpha_metric_label,
      levels = alpha_metric_key %>%
        arrange(alpha_metric_order) %>%
        pull(alpha_metric_label)
    ),
    driver_label = factor(
      driver_label,
      levels = c("Temperature", "Soil organic matter", "pH")
    )
  ) %>%
  group_by(
    scenario,
    scenario_label,
    scenario_family_label,
    scenario_order,
    alpha_metric_label,
    driver_label
  ) %>%
  summarise(
    delta_mean = mean(delta_estimate, na.rm = TRUE),
    delta_p10 = quantile(delta_estimate, 0.10, na.rm = TRUE, names = FALSE),
    delta_p90 = quantile(delta_estimate, 0.90, na.rm = TRUE, names = FALSE),
    n_iter = n(),
    .groups = "drop"
  )

p_alpha_drivers_main <- alpha_driver_delta_pub %>%
  ggplot(aes(x = scenario_label, y = delta_mean, shape = driver_label)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed") +
  geom_errorbar(aes(ymin = delta_p10, ymax = delta_p90), width = 0.12, na.rm = TRUE) +
  geom_point(size = 2.2, na.rm = TRUE) +
  facet_grid(alpha_metric_label ~ scenario_family_label, scales = "free_x", space = "free_x") +
  labs(
    title = "Taxonomic uncertainty can alter environment–diversity relationships",
    subtitle = "Difference in standardised coefficients relative to the appropriate baseline",
    x = NULL,
    y = expression(Delta * " standardised coefficient"),
    shape = "Driver"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PUB_DIR, "Fig4_main_alpha_driver_delta_coefficients.png"),
  p_alpha_drivers_main,
  width = 10.5,
  height = 6.8,
  dpi = 600
)

ggsave(
  file.path(PUB_DIR, "Fig4_main_alpha_driver_delta_coefficients.pdf"),
  p_alpha_drivers_main,
  width = 10.5,
  height = 6.8
)

# ------------------------------------------------------------
# 7. Supplementary full heatmap
# ------------------------------------------------------------

p_robustness_supp <- stability_long %>%
  filter(!scenario %in% baseline_scenarios) %>%
  inner_join(metric_key, by = "inference_label") %>%
  group_by(scenario_family, scenario, metric_label, metric_order) %>%
  summarise(
    stability_mean = mean(stability, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    stability_clipped = pmin(1, pmax(0, stability_mean)),
    loss_of_stability = 1 - stability_clipped,
    scenario = factor(scenario, levels = scenario_order),
    metric_label = factor(metric_label, levels = rev(metric_key$metric_label))
  ) %>%
  ggplot(aes(x = scenario, y = metric_label, fill = loss_of_stability)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(. ~ scenario_family, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    na.value = "grey90",
    name = "Loss of\nstability"
  ) +
  labs(
    title = "Full sensitivity analysis across all taxonomic uncertainty scenarios",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "right"
  )

ggsave(
  file.path(PUB_DIR, "FigS_full_robustness_heatmap.png"),
  p_robustness_supp,
  width = 14,
  height = 5.8,
  dpi = 600
)

ggsave(
  file.path(PUB_DIR, "FigS_full_robustness_heatmap.pdf"),
  p_robustness_supp,
  width = 14,
  height = 5.8
)


# ============================================================
# FIGURE 5 — Consequences for taxonomic inventories
# Gamma diversity and species-list alteration
# ============================================================
# -----------------------------
# 5b. INVENTAIRES TAXONOMIQUES PAR SCENARIO
# -----------------------------

scenario_inventory_by_iter <- scenario_data %>%
  distinct(
    scenario_family,
    scenario,
    baseline_scenario,
    unit_type,
    iter,
    taxon_unit
  ) %>%
  arrange(scenario_family, scenario, iter, taxon_unit)

readr::write_csv(
  scenario_inventory_by_iter,
  file.path(OUT_DIR, "scenario_inventory_by_iter.csv")
)

# Changements d'inventaire par rapport à la baseline propre à chaque famille.
# Pour chaque scénario, on mesure :
# - unités retenues depuis la baseline
# - unités perdues par rapport à la baseline
# - unités ajoutées par rapport à la baseline

baseline_inventory <- scenario_inventory_by_iter %>%
  filter(scenario == baseline_scenario) %>%
  group_by(scenario_family, baseline_scenario) %>%
  summarise(
    baseline_units = list(unique(taxon_unit)),
    gamma_baseline = n_distinct(taxon_unit),
    .groups = "drop"
  )

scenario_inventory_change <- scenario_inventory_by_iter %>%
  filter(scenario != baseline_scenario) %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, iter) %>%
  summarise(
    scenario_units = list(unique(taxon_unit)),
    gamma_scenario = n_distinct(taxon_unit),
    .groups = "drop"
  ) %>%
  left_join(
    baseline_inventory,
    by = c("scenario_family", "baseline_scenario")
  ) %>%
  mutate(
    n_retained = purrr::map2_int(
      scenario_units,
      baseline_units,
      ~ length(intersect(.x, .y))
    ),
    n_lost = purrr::map2_int(
      baseline_units,
      scenario_units,
      ~ length(setdiff(.x, .y))
    ),
    n_gained = purrr::map2_int(
      scenario_units,
      baseline_units,
      ~ length(setdiff(.x, .y))
    ),
    gamma_change = gamma_scenario - gamma_baseline,
    gamma_change_pct = 100 * gamma_change / gamma_baseline,
    lost_pct = 100 * n_lost / gamma_baseline,
    gained_pct = 100 * n_gained / gamma_baseline
  ) %>%
  select(
    scenario_family,
    scenario,
    baseline_scenario,
    unit_type,
    iter,
    gamma_baseline,
    gamma_scenario,
    gamma_change,
    gamma_change_pct,
    n_retained,
    n_lost,
    n_gained,
    lost_pct,
    gained_pct
  )

readr::write_csv(
  scenario_inventory_change,
  file.path(OUT_DIR, "scenario_inventory_change.csv")
)

message("Inventaire par scénario écrit : scenario_inventory_by_iter.csv")
message("Changements d'inventaire écrits : scenario_inventory_change.csv")

message("Publication-ready figures saved in: ", PUB_DIR)
