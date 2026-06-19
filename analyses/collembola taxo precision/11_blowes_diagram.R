# ============================================================
# 11_blowes_diagram — diagramme alpha–gamma–occupation type Blowes
# Scénarios RMQS : changement apparent induit par l'incertitude taxonomique
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)

WORKFLOW_DIR <- if (!is.na(.workflow_file)) {
  dirname(.workflow_file)
} else {
  getOption("taxo.workflow_dir", "analyses/collembola taxo precision")
}

if (!file.exists(file.path(WORKFLOW_DIR, "00_config.R"))) {
  WORKFLOW_DIR <- "analyses/collembola taxo precision"
}

source(file.path(WORKFLOW_DIR, "00_config.R"))
source(file.path(WORKFLOW_DIR, "functions_taxonomic_uncertainty.R"))

load_required_step("05_compute_alpha_metrics")

STEP_ID <- "11_blowes_diagram"

if (exists("message_header", mode = "function")) {
  message_header("11_blowes_diagram — alpha–gamma–occupation scenario space")
} else {
  message("11_blowes_diagram — alpha–gamma–occupation scenario space")
}

# ------------------------------------------------------------
# 1. Choix des scénarios à afficher dans la figure principale
# ------------------------------------------------------------

BLOWES_MAIN_SCENARIOS <- c(
  "species_same_genus_observed_10pct",
  "species_same_genus_observed_rare_weighted_mean_10pct",
  "species_same_genus_taxref_mainland_10pct",
  "rtu_rare_species_to_genus",
  "rtu_species_only_drop_unresolved",
  "rtu_genus_level",
  "rtu_family_level"
)

scenario_blowes_labels <- tibble::tribble(
  ~scenario, ~scenario_label, ~scenario_group_label,
  "species_same_genus_observed_10pct", "Diffuse congeneric\nerror 10%", "Species-level error",
  "species_same_genus_observed_rare_weighted_mean_10pct", "Rare-weighted\nerror 10%", "Species-level error",
  "species_same_genus_taxref_mainland_10pct", "TAXREF congeneric\nerror 10%", "Species-level error",
  "rtu_rare_species_to_genus", "Rare species\nreported as genus", "Workflow decision",
  "rtu_species_only_drop_unresolved", "Drop unresolved\nrecords", "Workflow decision",
  "rtu_genus_level", "Genus-level\nresolution", "Taxonomic coarsening",
  "rtu_family_level", "Family-level\nresolution", "Taxonomic coarsening"
)

scenario_group_levels <- c("Species-level error", "Workflow decision", "Taxonomic coarsening")

# ------------------------------------------------------------
# 2. Calcul des métriques alpha, gamma et occupation moyenne
# ------------------------------------------------------------
# alpha = richesse locale moyenne q0
# gamma = nombre total d'unités taxonomiques dans le jeu de données
# occupation moyenne = alpha / gamma
#
# Dans le cadre de Blowes, la droite 1:1 correspond à une occupation
# moyenne inchangée. Si alpha augmente plus que gamma, ou si gamma baisse
# plus que alpha, l'occupation moyenne augmente : homogénéisation apparente.
# Si gamma augmente plus que alpha, l'occupation moyenne baisse :
# différenciation apparente.

blowes_metrics_by_iter <- scenario_summary_by_iter %>%
  mutate(
    mean_occupancy = mean_local_q0 / gamma_taxon_units
  )

baseline_metrics <- blowes_metrics_by_iter %>%
  filter(scenario == baseline_scenario) %>%
  transmute(
    scenario_family,
    baseline_scenario = scenario,
    alpha_baseline = mean_local_q0,
    gamma_baseline = gamma_taxon_units,
    occupancy_baseline = mean_occupancy
  ) %>%
  distinct()

blowes_delta_by_iter <- blowes_metrics_by_iter %>%
  left_join(baseline_metrics, by = c("scenario_family", "baseline_scenario")) %>%
  mutate(
    delta_alpha = mean_local_q0 - alpha_baseline,
    delta_gamma = gamma_taxon_units - gamma_baseline,
    delta_occupancy = mean_occupancy - occupancy_baseline,
    rel_delta_alpha_pct = 100 * delta_alpha / alpha_baseline,
    rel_delta_gamma_pct = 100 * delta_gamma / gamma_baseline,
    rel_delta_occupancy_pct = 100 * delta_occupancy / occupancy_baseline,
    beta_signature = case_when(
      is.na(rel_delta_alpha_pct) | is.na(rel_delta_gamma_pct) ~ NA_character_,
      abs(rel_delta_alpha_pct) < 0.5 & abs(rel_delta_gamma_pct) < 0.5 ~ "No meaningful change",
      rel_delta_occupancy_pct > 0 ~ "Apparent homogenisation",
      rel_delta_occupancy_pct < 0 ~ "Apparent differentiation",
      TRUE ~ "No meaningful change"
    )
  ) %>%
  left_join(scenario_blowes_labels, by = "scenario") %>%
  mutate(
    scenario_label = coalesce(scenario_label, scenario),
    scenario_group_label = coalesce(scenario_group_label, scenario_family),
    scenario_group_label = factor(scenario_group_label, levels = scenario_group_levels),
    beta_signature = factor(
      beta_signature,
      levels = c("Apparent differentiation", "No meaningful change", "Apparent homogenisation")
    )
  )

safe_q <- function(x, p) {
  if (all(is.na(x))) NA_real_ else as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

blowes_delta_summary <- blowes_delta_by_iter %>%
  group_by(scenario_family, scenario, baseline_scenario, unit_type, scenario_label, scenario_group_label) %>%
  summarise(
    n_iter = n(),
    rel_delta_alpha_med = median(rel_delta_alpha_pct, na.rm = TRUE),
    rel_delta_alpha_p10 = safe_q(rel_delta_alpha_pct, 0.10),
    rel_delta_alpha_p90 = safe_q(rel_delta_alpha_pct, 0.90),
    rel_delta_gamma_med = median(rel_delta_gamma_pct, na.rm = TRUE),
    rel_delta_gamma_p10 = safe_q(rel_delta_gamma_pct, 0.10),
    rel_delta_gamma_p90 = safe_q(rel_delta_gamma_pct, 0.90),
    rel_delta_occupancy_med = median(rel_delta_occupancy_pct, na.rm = TRUE),
    rel_delta_occupancy_p10 = safe_q(rel_delta_occupancy_pct, 0.10),
    rel_delta_occupancy_p90 = safe_q(rel_delta_occupancy_pct, 0.90),
    alpha_baseline = first(alpha_baseline),
    gamma_baseline = first(gamma_baseline),
    occupancy_baseline = first(occupancy_baseline),
    beta_signature = case_when(
      is.na(rel_delta_alpha_med) | is.na(rel_delta_gamma_med) ~ NA_character_,
      abs(rel_delta_alpha_med) < 0.5 & abs(rel_delta_gamma_med) < 0.5 ~ "No meaningful change",
      rel_delta_occupancy_med > 0 ~ "Apparent homogenisation",
      rel_delta_occupancy_med < 0 ~ "Apparent differentiation",
      TRUE ~ "No meaningful change"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    beta_signature = factor(
      beta_signature,
      levels = c("Apparent differentiation", "No meaningful change", "Apparent homogenisation")
    )
  )

blowes_delta_summary_main <- blowes_delta_summary %>%
  filter(scenario %in% BLOWES_MAIN_SCENARIOS) %>%
  mutate(
    scenario = factor(scenario, levels = BLOWES_MAIN_SCENARIOS),
    scenario_label = factor(
      scenario_label,
      levels = scenario_blowes_labels$scenario_label[match(BLOWES_MAIN_SCENARIOS, scenario_blowes_labels$scenario)]
    )
  ) %>%
  arrange(scenario)

readr::write_csv(blowes_delta_by_iter, file.path(OUT_DIR, "blowes_delta_by_iter.csv"))
readr::write_csv(blowes_delta_summary, file.path(OUT_DIR, "blowes_delta_summary_all.csv"))
readr::write_csv(blowes_delta_summary_main, file.path(OUT_DIR, "blowes_delta_summary_main.csv"))

# ------------------------------------------------------------
# 3. Figure principale
# ------------------------------------------------------------

range_vals <- range(
  c(
    0,
    blowes_delta_summary_main$rel_delta_alpha_p10,
    blowes_delta_summary_main$rel_delta_alpha_p90,
    blowes_delta_summary_main$rel_delta_gamma_p10,
    blowes_delta_summary_main$rel_delta_gamma_p90
  ),
  na.rm = TRUE
)

if (!all(is.finite(range_vals))) range_vals <- c(-10, 10)

pad <- diff(range_vals) * 0.15
if (!is.finite(pad) || pad == 0) pad <- 2
plot_lims <- c(floor(range_vals[1] - pad), ceiling(range_vals[2] + pad))

fig_blowes_alpha_gamma <- ggplot(
  blowes_delta_summary_main,
  aes(x = rel_delta_alpha_med, y = rel_delta_gamma_med)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_segment(
    aes(
      x = rel_delta_alpha_p10,
      xend = rel_delta_alpha_p90,
      y = rel_delta_gamma_med,
      yend = rel_delta_gamma_med
    ),
    linewidth = 0.45,
    alpha = 0.55
  ) +
  geom_segment(
    aes(
      x = rel_delta_alpha_med,
      xend = rel_delta_alpha_med,
      y = rel_delta_gamma_p10,
      yend = rel_delta_gamma_p90
    ),
    linewidth = 0.45,
    alpha = 0.55
  ) +
  geom_segment(
    aes(
      x = 0, y = 0,
      xend = rel_delta_alpha_med,
      yend = rel_delta_gamma_med,
      colour = beta_signature
    ),
    arrow = grid::arrow(length = grid::unit(0.12, "cm")),
    linewidth = 0.45,
    alpha = 0.55
  ) +
  geom_point(
    aes(colour = beta_signature, shape = scenario_group_label),
    size = 3.2
  ) +
  annotate(
    "label",
    x = plot_lims[1] + 0.05 * diff(plot_lims),
    y = plot_lims[2] - 0.08 * diff(plot_lims),
    label = "Apparent differentiation\nmean occupancy decreases",
    hjust = 0,
    vjust = 1,
    label.size = 0.15,
    size = 3
  ) +
  annotate(
    "label",
    x = plot_lims[2] - 0.05 * diff(plot_lims),
    y = plot_lims[1] + 0.08 * diff(plot_lims),
    label = "Apparent homogenisation\nmean occupancy increases",
    hjust = 1,
    vjust = 0,
    label.size = 0.15,
    size = 3
  ) +
  coord_equal(xlim = plot_lims, ylim = plot_lims, clip = "off") +
  scale_x_continuous(labels = scales::label_number(suffix = "%")) +
  scale_y_continuous(labels = scales::label_number(suffix = "%")) +
  labs(
    title = "Taxonomic uncertainty scenarios in alpha–gamma space",
    subtitle = "Changes are relative to the appropriate scenario-specific baseline; dashed line = unchanged mean occupancy",
    x = expression(paste("Change in mean local richness, ", Delta, alpha)),
    y = expression(paste("Change in regional richness, ", Delta, gamma)),
    colour = "Apparent beta-diversity signature",
    shape = "Scenario family"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank()
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  fig_blowes_alpha_gamma <- fig_blowes_alpha_gamma +
    ggrepel::geom_text_repel(
      aes(label = scenario_label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.alpha = 0.35,
      show.legend = FALSE
    )
} else {
  fig_blowes_alpha_gamma <- fig_blowes_alpha_gamma +
    geom_text(
      aes(label = scenario_label),
      size = 3,
      nudge_x = 0.6,
      nudge_y = 0.6,
      check_overlap = TRUE,
      show.legend = FALSE
    )
}

ggsave(
  filename = file.path(PUB_DIR, "FigX_blowes_alpha_gamma_scenario_space.pdf"),
  plot = fig_blowes_alpha_gamma,
  width = 9,
  height = 6.8
)

ggsave(
  filename = file.path(PUB_DIR, "FigX_blowes_alpha_gamma_scenario_space.png"),
  plot = fig_blowes_alpha_gamma,
  width = 9,
  height = 6.8,
  dpi = 300
)

# ------------------------------------------------------------
# 4. Version facettée par famille de scénarios, utile en supplément
# ------------------------------------------------------------

fig_blowes_alpha_gamma_facets <- fig_blowes_alpha_gamma +
  facet_wrap(~ scenario_group_label, nrow = 1) +
  labs(
    title = "Taxonomic uncertainty scenarios in alpha–gamma space",
    subtitle = "Facetted by scenario family"
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(PUB_DIR, "FigSX_blowes_alpha_gamma_scenario_space_facets.pdf"),
  plot = fig_blowes_alpha_gamma_facets,
  width = 11,
  height = 5.8
)

ggsave(
  filename = file.path(PUB_DIR, "FigSX_blowes_alpha_gamma_scenario_space_facets.png"),
  plot = fig_blowes_alpha_gamma_facets,
  width = 11,
  height = 5.8,
  dpi = 300
)

# ------------------------------------------------------------
# 5. Sauvegarde cache
# ------------------------------------------------------------

save_step(STEP_ID, c(
  "BLOWES_MAIN_SCENARIOS",
  "scenario_blowes_labels",
  "blowes_metrics_by_iter",
  "baseline_metrics",
  "blowes_delta_by_iter",
  "blowes_delta_summary",
  "blowes_delta_summary_main",
  "fig_blowes_alpha_gamma",
  "fig_blowes_alpha_gamma_facets"
))

message("Figures Blowes écrites dans : ", PUB_DIR)
