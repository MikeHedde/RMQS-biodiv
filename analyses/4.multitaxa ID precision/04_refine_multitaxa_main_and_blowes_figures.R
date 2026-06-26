# =============================================================================
# 04_refine_multitaxa_main_and_blowes_figures.R
# Refined publication figures for the RMQS 2024 taxonomic-uncertainty paper.
#
# Replaces the first visual drafts of:
#   Fig1_taxonomic_context
#   Fig3_directional_workflow_effects
#   Fig4_GDM_ecological_inference
#   Fig5_multitaxa_Blowes_synthesis
# and creates two improved supplementary GDM figures.
#
# It reads existing outputs only. No simulations or GDMs are rerun.
# =============================================================================

# ---- packages ----------------------------------------------------------------
pkgs <- c("dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork", "stringr")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install first: install.packages(c(", paste(sprintf('"%s"', missing), collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(stringr)
})

# ---- paths -------------------------------------------------------------------
AUDIT_DIR  <- "outputs_multitaxa_2024"
RESULT_DIR <- file.path(AUDIT_DIR, "uncertainty_results")
FIG_DIR    <- file.path(RESULT_DIR, "publication_figures")
MAIN_DIR   <- file.path(FIG_DIR, "main_text")
SUPP_DIR   <- file.path(FIG_DIR, "supplement")
dir.create(MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

MAIN_ASSEMBLAGES <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "carabidae_pitfall",
  "formicidae_pitfall",
  "isopoda_pitfall"
)

ERROR_SCENARIOS <- c(
  "species_observed_congeneric_10pct",
  "species_observed_rare_weighted_10pct"
)

WORKFLOW_SCENARIOS <- c(
  "rtu_rare_species_to_genus",
  "rtu_species_only_drop_unresolved",
  "rtu_genus_level",
  "rtu_family_level"
)

MAIN_SCENARIOS <- c(ERROR_SCENARIOS, WORKFLOW_SCENARIOS)

# ---- labels and publication style -------------------------------------------
assemblage_labels <- c(
  collembola_soil_core = "Collembola\nsoil cores",
  araneae_pitfall = "Araneae\npitfall traps",
  carabidae_pitfall = "Carabidae\npitfall traps",
  formicidae_pitfall = "Formicidae\npitfall traps",
  isopoda_pitfall = "Isopoda\npitfall traps",
  isopoda_hand_sorting = "Isopoda\nhand sorting",
  diplopoda_pitfall = "Diplopoda\npitfall traps",
  diplopoda_hand_sorting = "Diplopoda\nhand sorting"
)

# Colour-blind aware Okabe-Ito-derived palette; names match *display labels*.
assemblage_colours_id <- c(
  collembola_soil_core = "#009E73",
  araneae_pitfall = "#CC79A7",
  carabidae_pitfall = "#D55E00",
  formicidae_pitfall = "#0072B2",
  isopoda_pitfall = "#56B4E9",
  isopoda_hand_sorting = "#E69F00",
  diplopoda_pitfall = "#7A7A7A",
  diplopoda_hand_sorting = "#4D4D4D"
)
assemblage_colours <- setNames(
  unname(assemblage_colours_id),
  unname(assemblage_labels[names(assemblage_colours_id)])
)

scenario_info <- tibble::tribble(
  ~scenario, ~scenario_label, ~scenario_short, ~scenario_class,
  "species_observed_congeneric_10pct", "Observed-pool\nerror 10%", "Observed-pool\nerror 10%", "Identification error",
  "species_observed_rare_weighted_10pct", "Rare-weighted\nerror 10%", "Rare-weighted\nerror 10%", "Identification error",
  "rtu_rare_species_to_genus", "Rare taxa\nreported at genus", "Rare taxa at\ngenus", "Workflow decision",
  "rtu_species_only_drop_unresolved", "Drop unresolved\nrecords", "Drop unresolved\nrecords", "Workflow decision",
  "rtu_genus_level", "Genus-level\nresolution", "Genus-level\nresolution", "Taxonomic coarsening",
  "rtu_family_level", "Family-level\nresolution", "Family-level\nresolution", "Taxonomic coarsening"
)

scenario_levels <- scenario_info$scenario_label[
  match(MAIN_SCENARIOS, scenario_info$scenario)
]

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = rel(1.15)),
      plot.subtitle = element_text(hjust = 0, size = rel(0.95)),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", size = rel(0.92)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(7, 8, 7, 8)
    )
}

save_plot <- function(plot, filename, width, height, directory = MAIN_DIR) {
  ggsave(
    filename = file.path(directory, paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "mm",
    device = grDevices::pdf
  )
  ggsave(
    filename = file.path(directory, paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "mm",
    dpi = 500
  )
}

summarise_iter <- function(data, value_col, grouping) {
  data %>%
    filter(is.finite(.data[[value_col]])) %>%
    group_by(across(all_of(grouping))) %>%
    summarise(
      median = median(.data[[value_col]], na.rm = TRUE),
      p10 = quantile(.data[[value_col]], 0.10, na.rm = TRUE, names = FALSE),
      p90 = quantile(.data[[value_col]], 0.90, na.rm = TRUE, names = FALSE),
      n = n(),
      .groups = "drop"
    )
}

# ---- input tables ------------------------------------------------------------
need <- c(
  file.path(AUDIT_DIR, "assemblage_manifest.csv"),
  file.path(AUDIT_DIR, "taxonomic_resolution_audit_all.csv"),
  file.path(RESULT_DIR, "results_by_iter_all.csv")
)
if (any(!file.exists(need))) stop("Missing expected input file(s):\n", paste(need[!file.exists(need)], collapse = "\n"))

manifest <- read_csv(file.path(AUDIT_DIR, "assemblage_manifest.csv"), show_col_types = FALSE)
tax_audit <- read_csv(file.path(AUDIT_DIR, "taxonomic_resolution_audit_all.csv"), show_col_types = FALSE)
raw <- read_csv(file.path(RESULT_DIR, "results_by_iter_all.csv"), show_col_types = FALSE)

metadata <- manifest %>%
  transmute(
    assemblage_id,
    group_label,
    method,
    n_sampled_stations,
    n_positive_stations,
    n_species,
    n_genera,
    n_families,
    share_species_abundance,
    share_genus_abundance,
    share_family_or_coarser_abundance,
    juvenile_share,
    assemblage_label = recode(assemblage_id, !!!assemblage_labels, .default = assemblage_id)
  )

raw_contrasts <- raw %>%
  filter(!is.na(baseline_scenario), scenario != baseline_scenario) %>%
  left_join(metadata, by = "assemblage_id") %>%
  left_join(scenario_info, by = "scenario") %>%
  mutate(
    scenario_label = factor(scenario_label, levels = scenario_levels),
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])
    )
  )

# =============================================================================
# Figure 1 — revised taxonomic / ontogenetic context
# =============================================================================
resolution_pal <- c(
  "Species" = "#0072B2",
  "Genus" = "#E69F00",
  "Family" = "#009E73",
  "Order" = "#CC79A7",
  "Other / unusable" = "#999999"
)

meta_main <- metadata %>%
  filter(assemblage_id %in% MAIN_ASSEMBLAGES) %>%
  mutate(
    assemblage_label = factor(
      assemblage_label,
      levels = rev(unname(assemblage_labels[MAIN_ASSEMBLAGES]))
    )
  )

p1a <- tax_audit %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    resolution %in% c("species", "genus", "family", "order", "other_or_unusable")
  ) %>%
  left_join(meta_main %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
  mutate(
    resolution = recode(
      resolution,
      species = "Species",
      genus = "Genus",
      family = "Family",
      order = "Order",
      other_or_unusable = "Other / unusable"
    ),
    resolution = factor(resolution, levels = names(resolution_pal))
  ) %>%
  ggplot(aes(x = 100 * abundance_share, y = assemblage_label, fill = resolution)) +
  geom_col(width = 0.74, colour = "white", linewidth = 0.28) +
  scale_fill_manual(values = resolution_pal, name = "Reported resolution") +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    labels = \(x) paste0(x, "%"), expand = c(0, 0)
  ) +
  labs(
    title = "Reported taxonomic resolution",
    x = "Share of total abundance",
    y = NULL
  ) +
  theme_paper() +
  theme(panel.grid.major.y = element_blank())

context_main <- meta_main %>%
  transmute(
    assemblage_label,
    `Not resolved to species` = 100 * (share_genus_abundance + share_family_or_coarser_abundance),
    `Explicit non-adults` = 100 * juvenile_share
  ) %>%
  pivot_longer(-assemblage_label, names_to = "source", values_to = "share_pct") %>%
  mutate(source = factor(source, levels = c("Not resolved to species", "Explicit non-adults")))

p1b <- ggplot(context_main, aes(x = share_pct, y = assemblage_label)) +
  geom_segment(aes(x = 0, xend = share_pct, yend = assemblage_label), colour = "grey78", linewidth = 0.45) +
  geom_point(size = 2.6, colour = "#333333") +
  geom_text(
    aes(label = sprintf("%.0f%%", share_pct)),
    hjust = -0.25, size = 2.7, colour = "#333333"
  ) +
  facet_wrap(~ source, ncol = 1, scales = "free_x") +
  scale_x_continuous(
    labels = \(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(title = "Sources of taxonomic uncertainty", x = "Share of total abundance", y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

p_fig1 <- (p1a + p1b) +
  plot_layout(widths = c(1.55, 1)) +
  plot_annotation(
    title = "Taxonomic and ontogenetic context of the focal RMQS assemblages",
    subtitle = "Each group × protocol dataset is analysed independently.",
    tag_levels = "A"
  )

save_plot(p_fig1, "Fig1_taxonomic_context_refined", 180, 105)

# =============================================================================
# Figure 3 — revised directional workflow effects
# =============================================================================
workflow_long <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% WORKFLOW_SCENARIOS
  ) %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label,
    `Regional richness (γ)` = gamma_change_pct,
    `Mean Bray–Curtis` = bray_mean_change_pct,
    `Mean Sørensen` = sorensen_mean_change_pct,
    `Explained environmental turnover (GDM)` = gdm_bray_delta_explained
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label, scenario, scenario_label),
    names_to = "response",
    values_to = "change"
  ) %>%
  summarise_iter(
    value_col = "change",
    grouping = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "response")
  ) %>%
  mutate(
    response = factor(
      response,
      levels = c(
        "Regional richness (γ)",
        "Mean Bray–Curtis",
        "Mean Sørensen",
        "Explained environmental turnover (GDM)"
      )
    ),
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])
    ),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label[match(WORKFLOW_SCENARIOS, scenario_info$scenario)])
  )

p_fig3 <- ggplot(
  workflow_long,
  aes(x = scenario_label, y = median, colour = assemblage_label)
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_errorbar(
    aes(ymin = p10, ymax = p90),
    width = 0, linewidth = 0.42,
    position = position_dodge(width = 0.45),
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.1,
    position = position_dodge(width = 0.45),
    na.rm = TRUE
  ) +
  facet_wrap(~ response, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
  labs(
    title = "Workflow choices generate directional changes in inventories, beta diversity and environmental turnover",
    subtitle = "Points are medians; intervals indicate the 10th–90th percentiles across stochastic iterations.",
    x = NULL,
    y = "Change relative to the appropriate baseline (%)\n(explained turnover: percentage points)"
  ) +
  theme_paper(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.spacing = grid::unit(5, "mm")
  )

save_plot(p_fig3, "Fig3_directional_workflow_effects_refined", 190, 145)

# =============================================================================
# Figure 4 — revised ecological-inference heatmaps
# =============================================================================
gdm_long <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% MAIN_SCENARIOS
  ) %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label,
    `Bray–Curtis GDM` = gdm_bray_delta_explained,
    `Sørensen GDM` = gdm_sorensen_delta_explained,
    `Bray–Curtis prediction loss` = 1 - gdm_bray_predicted_stability,
    `Sørensen prediction loss` = 1 - gdm_sorensen_predicted_stability
  )

gdm_fit <- gdm_long %>%
  pivot_longer(
    cols = c(`Bray–Curtis GDM`, `Sørensen GDM`),
    names_to = "distance", values_to = "delta_explained"
  ) %>%
  summarise_iter(
    value_col = "delta_explained",
    grouping = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "distance")
  ) %>%
  mutate(
    assemblage_label = factor(assemblage_label, levels = rev(unname(assemblage_labels[MAIN_ASSEMBLAGES]))),
    scenario_label = factor(scenario_label, levels = scenario_levels)
  )

gdm_pred <- gdm_long %>%
  pivot_longer(
    cols = c(`Bray–Curtis prediction loss`, `Sørensen prediction loss`),
    names_to = "distance", values_to = "prediction_loss"
  ) %>%
  mutate(distance = str_remove(distance, " prediction loss")) %>%
  summarise_iter(
    value_col = "prediction_loss",
    grouping = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "distance")
  ) %>%
  mutate(
    assemblage_label = factor(assemblage_label, levels = rev(unname(assemblage_labels[MAIN_ASSEMBLAGES]))),
    scenario_label = factor(scenario_label, levels = scenario_levels)
  )

p4a <- ggplot(gdm_fit, aes(x = scenario_label, y = assemblage_label, fill = median)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = ifelse(is.na(median), "", sprintf("%+.1f", median))), size = 2.6) +
  facet_wrap(~ distance, ncol = 1) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, na.value = "grey86",
    name = "Δ explained deviance\n(percentage points)"
  ) +
  labs(title = "GDM fit", x = NULL, y = NULL) +
  theme_paper(base_size = 8.2) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

p4b <- ggplot(gdm_pred, aes(x = scenario_label, y = assemblage_label, fill = median)) +
  geom_tile(colour = "white", linewidth = 0.38) +
  geom_text(aes(label = ifelse(is.na(median), "", sprintf("%.2f", median))), size = 2.6) +
  facet_wrap(~ distance, ncol = 1) +
  scale_fill_gradientn(
    colours = c("#FFF7EC", "#FEC44F", "#F03B20", "#7F0000"),
    limits = c(0, 1), oob = squish,
    na.value = "grey86",
    name = "Loss of predicted-\ndissimilarity stability"
  ) +
  labs(title = "GDM prediction stability", x = NULL, y = NULL) +
  theme_paper(base_size = 8.2) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

p_fig4 <- (p4a + p4b) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Taxonomic workflows can modify environmental turnover inference",
    subtitle = "Each scenario is compared with the GDM fitted to the same assemblage and sites.",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_plot(p_fig4, "Fig4_GDM_ecological_inference_refined", 190, 150)

# =============================================================================
# Figure 5 — revised Blowes alpha-gamma-occupancy synthesis
# =============================================================================
summarise_blowes <- function(data) {
  data %>%
    filter(is.finite(mean_q0_change_pct), is.finite(gamma_change_pct)) %>%
    mutate(
      alpha_ratio = 1 + mean_q0_change_pct / 100,
      gamma_ratio = 1 + gamma_change_pct / 100
    ) %>%
    filter(alpha_ratio > 0, gamma_ratio > 0) %>%
    mutate(occupancy_ratio = alpha_ratio / gamma_ratio) %>%
    group_by(assemblage_id, assemblage_label, scenario, scenario_label) %>%
    summarise(
      alpha = 100 * (median(alpha_ratio) - 1),
      gamma = 100 * (median(gamma_ratio) - 1),
      alpha_p10 = 100 * (quantile(alpha_ratio, .10) - 1),
      alpha_p90 = 100 * (quantile(alpha_ratio, .90) - 1),
      gamma_p10 = 100 * (quantile(gamma_ratio, .10) - 1),
      gamma_p90 = 100 * (quantile(gamma_ratio, .90) - 1),
      occupancy_change = 100 * (median(occupancy_ratio) - 1),
      .groups = "drop"
    )
}

blowes_workflow <- raw_contrasts %>%
  filter(assemblage_id %in% MAIN_ASSEMBLAGES, scenario %in% WORKFLOW_SCENARIOS) %>%
  summarise_blowes() %>%
  mutate(
    assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label[match(WORKFLOW_SCENARIOS, scenario_info$scenario)])
  )

blowes_error <- raw_contrasts %>%
  filter(assemblage_id %in% MAIN_ASSEMBLAGES, scenario %in% ERROR_SCENARIOS) %>%
  summarise_blowes() %>%
  mutate(
    assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label[match(ERROR_SCENARIOS, scenario_info$scenario)])
  )

limits_from <- function(x, minimum_span) {
  vals <- c(x$alpha_p10, x$alpha_p90, x$gamma_p10, x$gamma_p90, 0)
  vals <- vals[is.finite(vals)]
  lo <- min(vals); hi <- max(vals)
  span <- max(hi - lo, minimum_span)
  c(floor((lo - .08 * span) / 10) * 10, ceiling((hi + .08 * span) / 10) * 10)
}

plot_blowes <- function(data, title, subtitle, ncol, lims) {
  ggplot(data, aes(x = alpha, y = gamma, colour = assemblage_label)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.32) +
    geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.32) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_segment(aes(x = alpha_p10, xend = alpha_p90, y = gamma, yend = gamma), linewidth = 0.45, na.rm = TRUE) +
    geom_segment(aes(x = alpha, xend = alpha, y = gamma_p10, yend = gamma_p90), linewidth = 0.45, na.rm = TRUE) +
    geom_point(size = 2.8) +
    facet_wrap(~ scenario_label, ncol = ncol) +
    coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
    scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Change in mean local richness, Δα (%)",
      y = "Change in regional richness, Δγ (%)"
    ) +
    theme_paper(base_size = 8.5) +
    theme(panel.spacing = grid::unit(4, "mm"))
}

p5a <- plot_blowes(
  blowes_workflow,
  "Workflow decisions and taxonomic coarsening",
  "Dashed line: unchanged mean occupancy. Points below the line indicate apparent homogenisation.",
  ncol = 2,
  lims = limits_from(blowes_workflow, 40)
)

p5b <- plot_blowes(
  blowes_error,
  "Moderate identification error within the observed pool",
  "Axes are rescaled because these scenarios cause smaller shifts than coarsening workflows.",
  ncol = 2,
  lims = limits_from(blowes_error, 30)
)

p_fig5 <- (p5a / p5b) +
  plot_layout(heights = c(1.35, 1)) +
  plot_annotation(
    title = "Taxonomic workflows shift apparent occupancy patterns in alpha–gamma space",
    subtitle = "Each point is a scenario median; crosshairs show 10th–90th percentiles across stochastic iterations.",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_plot(p_fig5, "Fig5_multitaxa_Blowes_synthesis_refined", 190, 225)

# =============================================================================
# Supplementary GDM figures
# =============================================================================

# S6: Baseline explanatory power, needed to contextualise absolute changes.
baseline_gdm <- raw %>%
  filter(scenario == baseline_scenario, assemblage_id %in% MAIN_ASSEMBLAGES) %>%
  left_join(metadata %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
  transmute(
    assemblage_label,
    `Bray–Curtis` = gdm_bray_explained_baseline,
    `Sørensen` = gdm_sorensen_explained_baseline
  ) %>%
  pivot_longer(-assemblage_label, names_to = "distance", values_to = "explained") %>%
  summarise_iter(
    value_col = "explained",
    grouping = c("assemblage_label", "distance")
  ) %>%
  mutate(assemblage_label = factor(assemblage_label, levels = rev(unname(assemblage_labels[MAIN_ASSEMBLAGES]))))

p_s6 <- ggplot(baseline_gdm, aes(x = median, y = assemblage_label, colour = distance)) +
  geom_errorbarh(aes(xmin = p10, xmax = p90), height = 0, linewidth = 0.45, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Bray–Curtis" = "#0072B2", "Sørensen" = "#D55E00"), name = "GDM") +
  scale_x_continuous(labels = \(x) paste0(round(x, 1), "%")) +
  labs(
    title = "Baseline explanatory power of environmental-turnover models",
    subtitle = "Explained deviance of GDMs fitted to baseline communities.",
    x = "Explained deviance",
    y = NULL
  ) +
  theme_paper()

save_plot(p_s6, "FigS6_baseline_GDM_explained_deviance", 170, 105, directory = SUPP_DIR)

# S7: Predictor composition shifts, only if terms were extracted.
terms_path <- file.path(RESULT_DIR, "gdm_terms_summary.csv")
if (file.exists(terms_path)) {
  gdm_terms <- read_csv(terms_path, show_col_types = FALSE)
  required_term_cols <- c(
    "assemblage_id", "scenario", "distance", "predictor",
    "baseline_ispline_weight_median", "scenario_ispline_weight_median"
  )
  
  if (all(required_term_cols %in% names(gdm_terms))) {
    predictor_shift <- gdm_terms %>%
      filter(assemblage_id %in% MAIN_ASSEMBLAGES, scenario %in% MAIN_SCENARIOS) %>%
      group_by(assemblage_id, scenario, distance) %>%
      mutate(
        baseline_share = baseline_ispline_weight_median / sum(baseline_ispline_weight_median, na.rm = TRUE),
        scenario_share = scenario_ispline_weight_median / sum(scenario_ispline_weight_median, na.rm = TRUE),
        delta_share_pp = 100 * (scenario_share - baseline_share)
      ) %>%
      ungroup() %>%
      left_join(metadata %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
      left_join(scenario_info %>% select(scenario, scenario_label), by = "scenario") %>%
      mutate(
        scenario_label = factor(scenario_label, levels = scenario_levels),
        assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES]))
      )
    
    p_s7 <- ggplot(predictor_shift, aes(x = predictor, y = delta_share_pp, colour = assemblage_label)) +
      geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.3) +
      geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.8, alpha = 0.9) +
      facet_grid(distance ~ scenario_label, scales = "free_y") +
      scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
      labs(
        title = "Scenario-induced shifts in the relative contribution of GDM predictors",
        subtitle = "Changes in each predictor's share of the summed I-spline coefficients.",
        x = NULL,
        y = "Change in relative contribution (percentage points)"
      ) +
      theme_paper(base_size = 7.4) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    
    save_plot(p_s7, "FigS7_GDM_predictor_composition_shifts", 245, 160, directory = SUPP_DIR)
  }
}

# export exact coordinate tables used for Figure 5
write_csv(blowes_workflow, file.path(FIG_DIR, "Blowes_workflow_summary_refined.csv"))
write_csv(blowes_error, file.path(FIG_DIR, "Blowes_error_summary_refined.csv"))

message("\nRefined figures written to:")
message(MAIN_DIR)
message(SUPP_DIR)
