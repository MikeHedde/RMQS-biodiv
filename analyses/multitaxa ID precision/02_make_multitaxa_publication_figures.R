# =============================================================================
# 02_make_multitaxa_publication_figures.R
# Publication-ready figures for RMQS 2024 multi-taxon taxonomic-uncertainty
# analysis.
#
# This script only reads output tables already produced by:
#   00_prepare_multitaxa_inputs_2024.R
#   01_run_multitaxa_uncertainty_2024_v1_5_gdm_dataframe_fix.R
#
# It does NOT rerun simulations or GDMs.
# =============================================================================

# ---- 0. Packages -------------------------------------------------------------
required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "forcats",
  "ggplot2", "scales", "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    "Install missing packages first:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ---- 1. Paths and user choices ----------------------------------------------
# Change these two paths only if your project structure differs.
AUDIT_DIR  <- "outputs_multitaxa_2024"
RESULT_DIR <- file.path(AUDIT_DIR, "uncertainty_results")
FIG_DIR    <- file.path(RESULT_DIR, "publication_figures")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "main_text"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "supplement"), recursive = TRUE, showWarnings = FALSE)

# Main-text assemblages: high-coverage datasets that carry the main narrative.
MAIN_ASSEMBLAGES <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "carabidae_pitfall",
  "formicidae_pitfall",
  "isopoda_pitfall"
)

# Secondary assemblages are retained for robustness / stage-related appendices.
SECONDARY_ASSEMBLAGES <- c(
  "diplopoda_hand_sorting",
  "diplopoda_pitfall",
  "isopoda_hand_sorting"
)

# Main-text scenarios. Family coarsening remains displayed as NA / grey where
# it is not meaningful (Carabidae and Formicidae are already single-family sets).
MAIN_SCENARIOS <- c(
  "species_observed_congeneric_10pct",
  "species_observed_rare_weighted_10pct",
  "rtu_rare_species_to_genus",
  "rtu_species_only_drop_unresolved",
  "rtu_genus_level",
  "rtu_family_level"
)

WORKFLOW_SCENARIOS <- c(
  "rtu_rare_species_to_genus",
  "rtu_species_only_drop_unresolved",
  "rtu_genus_level",
  "rtu_family_level"
)

ALL_SCENARIOS_ORDERED <- c(
  "adult_only_strict",
  "species_observed_congeneric_5pct",
  "species_observed_congeneric_10pct",
  "species_observed_rare_weighted_10pct",
  "rtu_rare_species_to_genus",
  "rtu_species_only_drop_unresolved",
  "rtu_genus_level",
  "rtu_family_level"
)

# ---- 2. Labels, colours, theme ----------------------------------------------
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

assemblage_colours <- c(
  collembola_soil_core = "#009E73",
  araneae_pitfall = "#CC79A7",
  carabidae_pitfall = "#D55E00",
  formicidae_pitfall = "#0072B2",
  isopoda_pitfall = "#56B4E9",
  isopoda_hand_sorting = "#E69F00",
  diplopoda_pitfall = "#7A7A7A",
  diplopoda_hand_sorting = "#4D4D4D"
)

scenario_info <- tibble::tribble(
  ~scenario, ~scenario_label, ~scenario_short, ~scenario_class,
  "adult_only_strict", "Adults only", "Adults only", "Ontogenetic filtering",
  "species_observed_congeneric_5pct", "Observed-pool\nerror 5%", "Observed-pool\nerror 5%", "Identification error",
  "species_observed_congeneric_10pct", "Observed-pool\nerror 10%", "Observed-pool\nerror 10%", "Identification error",
  "species_observed_rare_weighted_10pct", "Rare-weighted\nerror 10%", "Rare-weighted\nerror 10%", "Identification error",
  "rtu_rare_species_to_genus", "Rare taxa\nreported at genus", "Rare taxa at\ngenus", "Workflow decision",
  "rtu_species_only_drop_unresolved", "Drop unresolved\nrecords", "Drop unresolved\nrecords", "Workflow decision",
  "rtu_genus_level", "Genus-level\nresolution", "Genus-level\nresolution", "Taxonomic coarsening",
  "rtu_family_level", "Family-level\nresolution", "Family-level\nresolution", "Taxonomic coarsening"
)

resolution_colours <- c(
  species = "#0072B2",
  genus = "#E69F00",
  family = "#009E73",
  order = "#CC79A7",
  other_or_unusable = "#999999"
)

resolution_colours_labelled <- c(
  "Species" = "#0072B2",
  "Genus" = "#E69F00",
  "Family" = "#009E73",
  "Order" = "#CC79A7",
  "Other / unusable" = "#999999"
)

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.15), hjust = 0),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", size = rel(0.92)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_figure <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(FIG_DIR, "main_text", paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "mm",
    device = grDevices::pdf
  )
  ggsave(
    filename = file.path(FIG_DIR, "main_text", paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "mm",
    dpi = 450
  )
}

save_supplement <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(FIG_DIR, "supplement", paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "mm",
    device = grDevices::pdf
  )
  ggsave(
    filename = file.path(FIG_DIR, "supplement", paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "mm",
    dpi = 450
  )
}

summarise_distribution <- function(data, value_col, group_cols) {
  data %>%
    filter(is.finite(.data[[value_col]])) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      median = median(.data[[value_col]], na.rm = TRUE),
      p10 = quantile(.data[[value_col]], probs = 0.10, na.rm = TRUE, names = FALSE),
      p90 = quantile(.data[[value_col]], probs = 0.90, na.rm = TRUE, names = FALSE),
      n_iter = n(),
      .groups = "drop"
    )
}

# ---- 3. Read inputs ----------------------------------------------------------
required_files <- c(
  file.path(AUDIT_DIR, "assemblage_manifest.csv"),
  file.path(AUDIT_DIR, "taxonomic_resolution_audit_all.csv"),
  file.path(AUDIT_DIR, "stage_audit_all.csv"),
  file.path(RESULT_DIR, "results_by_iter_all.csv"),
  file.path(RESULT_DIR, "results_summary.csv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing required output files:\n",
    paste(" -", missing_files, collapse = "\n"),
    "\n\nRun the input-preparation and uncertainty scripts first."
  )
}

manifest <- readr::read_csv(file.path(AUDIT_DIR, "assemblage_manifest.csv"), show_col_types = FALSE)
tax_audit <- readr::read_csv(file.path(AUDIT_DIR, "taxonomic_resolution_audit_all.csv"), show_col_types = FALSE)
stage_audit <- readr::read_csv(file.path(AUDIT_DIR, "stage_audit_all.csv"), show_col_types = FALSE)
raw <- readr::read_csv(file.path(RESULT_DIR, "results_by_iter_all.csv"), show_col_types = FALSE)
results_summary <- readr::read_csv(file.path(RESULT_DIR, "results_summary.csv"), show_col_types = FALSE)

gdm_terms_path <- file.path(RESULT_DIR, "gdm_terms_summary.csv")
gdm_terms <- if (file.exists(gdm_terms_path)) {
  readr::read_csv(gdm_terms_path, show_col_types = FALSE)
} else {
  tibble()
}

# Derived labels / metadata used by all figures.
metadata <- manifest %>%
  transmute(
    assemblage_id,
    group_label,
    method,
    n_sampled_stations,
    n_positive_stations,
    total_abundance,
    n_species,
    n_genera,
    n_families,
    share_species_abundance,
    share_genus_abundance,
    share_family_or_coarser_abundance,
    juvenile_share,
    family_coarsening_applicable
  ) %>%
  mutate(
    assemblage_label = recode(assemblage_id, !!!assemblage_labels, .default = assemblage_id),
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[c(MAIN_ASSEMBLAGES, SECONDARY_ASSEMBLAGES)])
    )
  )

raw <- raw %>%
  left_join(metadata %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
  left_join(scenario_info, by = "scenario") %>%
  mutate(
    scenario_label = coalesce(scenario_label, scenario),
    scenario_short = coalesce(scenario_short, scenario),
    scenario_class = coalesce(scenario_class, scenario_family),
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[c(MAIN_ASSEMBLAGES, SECONDARY_ASSEMBLAGES)])
    )
  )

# Keep only genuine scenario-versus-baseline contrasts.
raw_contrasts <- raw %>%
  filter(!is.na(baseline_scenario), scenario != baseline_scenario)

# =============================================================================
# MAIN TEXT
# =============================================================================

# ---- Figure 1. Taxonomic context --------------------------------------------
main_meta <- metadata %>%
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
  left_join(main_meta %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
  mutate(
    resolution = factor(
      resolution,
      levels = c("species", "genus", "family", "order", "other_or_unusable"),
      labels = c("Species", "Genus", "Family", "Order", "Other / unusable")
    )
  ) %>%
  ggplot(aes(x = abundance_share * 100, y = assemblage_label, fill = resolution)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = resolution_colours_labelled, name = "Reported resolution") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Reported taxonomic resolution",
    x = "Share of total abundance",
    y = NULL
  ) +
  theme_paper() +
  theme(legend.position = "bottom")

context_long <- main_meta %>%
  transmute(
    assemblage_label,
    `Not resolved to species` = 100 * (share_genus_abundance + share_family_or_coarser_abundance),
    `Explicit non-adults` = 100 * juvenile_share
  ) %>%
  pivot_longer(-assemblage_label, names_to = "source", values_to = "share_pct") %>%
  mutate(
    source = factor(source, levels = c("Not resolved to species", "Explicit non-adults"))
  )

p1b <- ggplot(context_long, aes(x = share_pct, y = assemblage_label)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(size = 2.5, colour = "#333333") +
  facet_wrap(~ source, ncol = 1, scales = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title = "Sources of uncertainty",
    x = "Share of total abundance",
    y = NULL
  ) +
  theme_paper() +
  theme(legend.position = "none")

p_fig1 <- (p1a + p1b) +
  plot_layout(widths = c(1.55, 1)) +
  plot_annotation(
    title = "Taxonomic and ontogenetic context of the focal RMQS assemblages",
    subtitle = "Analyses are performed separately for each assemblage × protocol dataset.",
    tag_levels = "A"
  ) &
  theme(plot.title = element_text(face = "bold"))

save_figure(p_fig1, "Fig1_taxonomic_context", width = 180, height = 108)

# ---- Figure 2. Balanced alpha-beta-inference robustness heatmap -------------
metric_map_main <- tibble::tribble(
  ~metric_key, ~metric_label, ~domain,
  "q0_stability", "q0 richness", "Alpha diversity",
  "q1_stability", "Hill q1", "Alpha diversity",
  "q2_stability", "Hill q2", "Alpha diversity",
  "bray_stability", "Bray–Curtis", "Beta diversity",
  "sorensen_stability", "Sørensen", "Beta diversity",
  "ordination_procrustes_r2", "Ordination", "Beta diversity",
  "gdm_bray_predicted_stability", "GDM Bray", "Ecological inference",
  "gdm_sorensen_predicted_stability", "GDM Sørensen", "Ecological inference"
)

heatmap_main <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% MAIN_SCENARIOS
  ) %>%
  select(
    assemblage_id, assemblage_label, scenario, scenario_label, scenario_class,
    q0_stability, q1_stability, q2_stability,
    bray_stability, sorensen_stability, ordination_procrustes_r2,
    gdm_bray_predicted_stability, gdm_sorensen_predicted_stability
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label, scenario, scenario_label, scenario_class),
    names_to = "metric_key", values_to = "stability"
  ) %>%
  left_join(metric_map_main, by = "metric_key") %>%
  summarise_distribution(
    value_col = "stability",
    group_cols = c(
      "assemblage_id", "assemblage_label", "scenario", "scenario_label",
      "scenario_class", "metric_key", "metric_label", "domain"
    )
  )

# Explicit NA tiles make inapplicable scenarios visible rather than silently absent.
heatmap_grid <- tidyr::expand_grid(
  assemblage_id = MAIN_ASSEMBLAGES,
  scenario = MAIN_SCENARIOS,
  metric_key = metric_map_main$metric_key
) %>%
  left_join(metadata %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
  left_join(scenario_info %>% select(scenario, scenario_label, scenario_class), by = "scenario") %>%
  left_join(metric_map_main, by = "metric_key") %>%
  left_join(
    heatmap_main %>% select(assemblage_id, scenario, metric_key, median),
    by = c("assemblage_id", "scenario", "metric_key")
  ) %>%
  mutate(
    loss = pmin(1, pmax(0, 1 - median)),
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[match(MAIN_SCENARIOS, scenario_info$scenario)]
    ),
    metric_label = factor(
      metric_label,
      levels = rev(metric_map_main$metric_label)
    ),
    domain = factor(
      domain,
      levels = rev(c("Alpha diversity", "Beta diversity", "Ecological inference"))
    ),
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])
    )
  )

p_fig2 <- ggplot(heatmap_grid, aes(x = scenario_label, y = metric_label, fill = loss)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_vline(xintercept = c(2.5, 4.5), colour = "grey45", linewidth = 0.35) +
  facet_grid(domain ~ assemblage_label, scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = c("#FFF7EC", "#FEC44F", "#F03B20", "#7F0000"),
    limits = c(0, 1),
    oob = scales::squish,
    na.value = "grey86",
    name = "Loss of\nstability"
  ) +
  labs(
    title = "Taxonomic uncertainty affects alpha diversity, beta diversity and ecological inference differently",
    subtitle = "Loss of stability relative to the scenario-specific baseline; grey cells indicate non-applicable or non-estimable comparisons.",
    x = NULL,
    y = NULL
  ) +
  theme_paper(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.spacing = grid::unit(1.5, "mm")
  )

save_figure(p_fig2, "Fig2_cross_taxon_robustness", width = 190, height = 150)

# ---- Figure 3. Directional workflow effects, balanced across alpha/beta/GDM -
workflow_metrics <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% WORKFLOW_SCENARIOS
  ) %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label,
    gamma_change_pct,
    bray_mean_change_pct,
    sorensen_mean_change_pct,
    gdm_bray_relative_change_pct =
      if_else(
        is.finite(gdm_bray_explained_baseline) & abs(gdm_bray_explained_baseline) > 1e-8,
        100 * gdm_bray_delta_explained / gdm_bray_explained_baseline,
        NA_real_
      )
  ) %>%
  pivot_longer(
    cols = c(
      gamma_change_pct, bray_mean_change_pct,
      sorensen_mean_change_pct, gdm_bray_relative_change_pct
    ),
    names_to = "response_key", values_to = "change_pct"
  ) %>%
  mutate(
    response = recode(
      response_key,
      gamma_change_pct = "Δ regional richness (γ)",
      bray_mean_change_pct = "Δ mean Bray-Curtis",
      sorensen_mean_change_pct = "Δ mean Sørensen",
      gdm_bray_relative_change_pct = "Δ explained environmental turnover (GDM)"
    )
  ) %>%
  summarise_distribution(
    value_col = "change_pct",
    group_cols = c(
      "assemblage_id", "assemblage_label", "scenario", "scenario_label",
      "response_key", "response"
    )
  ) %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[match(WORKFLOW_SCENARIOS, scenario_info$scenario)]
    ),
    response = factor(
      as.character(response),
      levels = c(
        "Δ regional richness (γ)",
        "Δ mean Bray-Curtis",
        "Δ mean Sørensen",
        "Δ explained environmental turnover (GDM)"
      )
    ),
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])
    )
  )

p_fig3 <- ggplot(
  workflow_metrics,
  aes(x = scenario_label, y = median, colour = assemblage_label, group = assemblage_label)
) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey35") +
  geom_errorbar(
    aes(ymin = p10, ymax = p90),
    width = 0,
    linewidth = 0.35,
    position = position_dodge(width = 0.62),
    na.rm = TRUE
  ) +
  geom_point(
    size = 1.9,
    position = position_dodge(width = 0.62),
    na.rm = TRUE
  ) +
  facet_wrap(~ response, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
  labs(
    title = "Workflow choices produce directional changes in inventories, beta diversity and environmental turnover inference",
    subtitle = "Points are medians and bars are 10th–90th percentiles across stochastic iterations.",
    x = NULL,
    y = "Change relative to the appropriate baseline (%)"
  ) +
  theme_paper(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.spacing = grid::unit(4, "mm")
  )

save_figure(p_fig3, "Fig3_directional_workflow_effects", width = 190, height = 145)

# ---- Figure 4. GDM sensitivity detail ---------------------------------------
gdm_detail <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% MAIN_SCENARIOS
  ) %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label,
    bray_relative_explained_change = if_else(
      is.finite(gdm_bray_explained_baseline) & abs(gdm_bray_explained_baseline) > 1e-8,
      100 * gdm_bray_delta_explained / gdm_bray_explained_baseline,
      NA_real_
    ),
    sorensen_relative_explained_change = if_else(
      is.finite(gdm_sorensen_explained_baseline) & abs(gdm_sorensen_explained_baseline) > 1e-8,
      100 * gdm_sorensen_delta_explained / gdm_sorensen_explained_baseline,
      NA_real_
    ),
    bray_prediction_loss = 1 - gdm_bray_predicted_stability,
    sorensen_prediction_loss = 1 - gdm_sorensen_predicted_stability
  )

gdm_explained <- gdm_detail %>%
  pivot_longer(
    cols = c(bray_relative_explained_change, sorensen_relative_explained_change),
    names_to = "distance", values_to = "relative_change"
  ) %>%
  mutate(
    distance = recode(
      distance,
      bray_relative_explained_change = "Bray-Curtis GDM",
      sorensen_relative_explained_change = "Sørensen GDM"
    )
  ) %>%
  summarise_distribution(
    value_col = "relative_change",
    group_cols = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "distance")
  )

gdm_prediction <- gdm_detail %>%
  pivot_longer(
    cols = c(bray_prediction_loss, sorensen_prediction_loss),
    names_to = "distance", values_to = "prediction_loss"
  ) %>%
  mutate(
    distance = recode(
      distance,
      bray_prediction_loss = "Bray-Curtis GDM",
      sorensen_prediction_loss = "Sørensen GDM"
    )
  ) %>%
  summarise_distribution(
    value_col = "prediction_loss",
    group_cols = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "distance")
  )

gdm_explained <- gdm_explained %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[match(MAIN_SCENARIOS, scenario_info$scenario)]
    ),
    assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES]))
  )

gdm_prediction <- gdm_prediction %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[match(MAIN_SCENARIOS, scenario_info$scenario)]
    ),
    assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES]))
  )

p4a <- ggplot(gdm_explained, aes(x = scenario_label, y = assemblage_label, fill = median)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(is.na(median), "", sprintf("%+.0f", median)), colour = !is.na(median) & abs(median) > 12), size = 2.4) +
  facet_wrap(~ distance, ncol = 1) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, na.value = "grey86",
    name = "Relative change in\nexplained deviance (%)"
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  labs(title = "GDM fit", x = NULL, y = NULL) +
  theme_paper(base_size = 8) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

p4b <- ggplot(gdm_prediction, aes(x = scenario_label, y = assemblage_label, fill = median)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = ifelse(is.na(median), "", sprintf("%.2f", median)), colour = !is.na(median) & median > 0.25), size = 2.4) +
  facet_wrap(~ distance, ncol = 1) +
  scale_fill_gradientn(
    colours = c("#FFF7EC", "#FEC44F", "#F03B20", "#7F0000"),
    limits = c(0, 1), oob = scales::squish,
    na.value = "grey86", name = "Loss of GDM\nprediction stability"
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  labs(title = "GDM prediction stability", x = NULL, y = NULL) +
  theme_paper(base_size = 8) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

p_fig4 <- (p4a + p4b) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Taxonomic uncertainty can weaken environmental turnover inference",
    subtitle = "Effects are expressed relative to the baseline GDM fitted to the same assemblage and sites.",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_figure(p_fig4, "Fig4_GDM_ecological_inference", width = 190, height = 150)

# =============================================================================
# SUPPLEMENT
# =============================================================================

# ---- Figure S1. Full taxonomic + ontogenetic audit --------------------------
all_meta <- metadata %>%
  mutate(
    assemblage_label = factor(
      assemblage_label,
      levels = rev(unname(assemblage_labels[c(MAIN_ASSEMBLAGES, SECONDARY_ASSEMBLAGES)]))
    )
  )

ps1a <- tax_audit %>%
  filter(resolution %in% c("species", "genus", "family", "order", "other_or_unusable")) %>%
  left_join(all_meta %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
  mutate(
    resolution = factor(
      resolution,
      levels = c("species", "genus", "family", "order", "other_or_unusable"),
      labels = c("Species", "Genus", "Family", "Order", "Other / unusable")
    )
  ) %>%
  ggplot(aes(x = abundance_share * 100, y = assemblage_label, fill = resolution)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = resolution_colours_labelled, name = "Reported resolution") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  labs(title = "All assemblages", x = "Share of total abundance", y = NULL) +
  theme_paper()

context_all <- all_meta %>%
  transmute(
    assemblage_label,
    `Not resolved to species` = 100 * (share_genus_abundance + share_family_or_coarser_abundance),
    `Explicit non-adults` = 100 * juvenile_share,
    `Prevalence` = 100 * n_positive_stations / n_sampled_stations
  ) %>%
  pivot_longer(-assemblage_label, names_to = "descriptor", values_to = "value_pct")

ps1b <- ggplot(context_all, aes(x = value_pct, y = assemblage_label)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(size = 2.1) +
  facet_wrap(~ descriptor, ncol = 1, scales = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Taxonomic, stage and occupancy context", x = NULL, y = NULL) +
  theme_paper() +
  theme(legend.position = "none")

p_s1 <- (ps1a + ps1b) +
  plot_layout(widths = c(1.55, 1)) +
  plot_annotation(tag_levels = "A")

save_supplement(p_s1, "FigS1_full_taxonomic_ontogenetic_audit", width = 190, height = 150)

# ---- Figure S2. Full robustness heatmap, including secondary assemblies ----
metric_map_full <- tibble::tribble(
  ~metric_key, ~metric_label, ~domain,
  "q0_stability", "q0 richness", "Alpha diversity",
  "q1_stability", "Hill q1", "Alpha diversity",
  "q2_stability", "Hill q2", "Alpha diversity",
  "bray_stability", "Bray–Curtis", "Beta diversity",
  "sorensen_stability", "Sørensen", "Beta diversity",
  "jaccard_stability", "Jaccard", "Beta diversity",
  "ordination_procrustes_r2", "Ordination", "Beta diversity",
  "gdm_bray_predicted_stability", "GDM Bray", "Ecological inference",
  "gdm_sorensen_predicted_stability", "GDM Sørensen", "Ecological inference"
)

heatmap_full_data <- raw_contrasts %>%
  filter(scenario %in% ALL_SCENARIOS_ORDERED) %>%
  select(
    assemblage_id, assemblage_label, scenario, scenario_label,
    q0_stability, q1_stability, q2_stability,
    bray_stability, sorensen_stability, jaccard_stability,
    ordination_procrustes_r2,
    gdm_bray_predicted_stability, gdm_sorensen_predicted_stability
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label, scenario, scenario_label),
    names_to = "metric_key", values_to = "stability"
  ) %>%
  left_join(metric_map_full, by = "metric_key") %>%
  summarise_distribution(
    value_col = "stability",
    group_cols = c("assemblage_id", "assemblage_label", "scenario", "scenario_label",
                   "metric_key", "metric_label", "domain")
  ) %>%
  mutate(
    loss = pmin(1, pmax(0, 1 - median)),
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[match(ALL_SCENARIOS_ORDERED, scenario_info$scenario)]
    ),
    metric_label = factor(metric_label, levels = rev(metric_map_full$metric_label)),
    domain = factor(domain, levels = rev(c("Alpha diversity", "Beta diversity", "Ecological inference")))
  )

p_s2 <- ggplot(heatmap_full_data, aes(x = scenario_label, y = metric_label, fill = loss)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(domain ~ assemblage_label, scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = c("#FFF7EC", "#FEC44F", "#F03B20", "#7F0000"),
    limits = c(0, 1), oob = scales::squish,
    na.value = "grey86", name = "Loss of\nstability"
  ) +
  labs(
    title = "Full cross-assemblage robustness analysis",
    subtitle = "Includes secondary assemblages, 5% error, Jaccard and ontogenetic filtering where applicable.",
    x = NULL, y = NULL
  ) +
  theme_paper(base_size = 7.3) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.spacing = grid::unit(1.1, "mm")
  )

save_supplement(p_s2, "FigS2_full_robustness_heatmap", width = 250, height = 190)

# ---- Figure S3. Full directional alpha and beta responses ------------------
directional_full <- raw_contrasts %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label,
    gamma_change_pct,
    mean_q0_change_pct,
    mean_q1_change_pct,
    mean_q2_change_pct,
    bray_mean_change_pct,
    sorensen_mean_change_pct,
    jaccard_mean_change_pct
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label, scenario, scenario_label),
    names_to = "response_key", values_to = "change_pct"
  ) %>%
  mutate(
    response = recode(
      response_key,
      gamma_change_pct = "Gamma richness",
      mean_q0_change_pct = "Mean q0 richness",
      mean_q1_change_pct = "Mean Hill q1",
      mean_q2_change_pct = "Mean Hill q2",
      bray_mean_change_pct = "Mean Bray-Curtis",
      sorensen_mean_change_pct = "Mean Sørensen",
      jaccard_mean_change_pct = "Mean Jaccard"
    )
  ) %>%
  summarise_distribution(
    value_col = "change_pct",
    group_cols = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "response")
  ) %>%
  filter(!is.na(median)) %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[match(ALL_SCENARIOS_ORDERED, scenario_info$scenario)]
    )
  )

p_s3 <- ggplot(
  directional_full,
  aes(x = scenario_label, y = median, colour = assemblage_label, group = assemblage_label)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
  geom_errorbar(
    aes(ymin = p10, ymax = p90),
    width = 0, linewidth = 0.28,
    position = position_dodge(width = 0.62),
    na.rm = TRUE
  ) +
  geom_point(size = 1.5, position = position_dodge(width = 0.62), na.rm = TRUE) +
  facet_wrap(~ response, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
  labs(
    title = "Full directional effects of taxonomic uncertainty",
    x = NULL, y = "Change relative to the appropriate baseline (%)"
  ) +
  theme_paper(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

save_supplement(p_s3, "FigS3_full_alpha_beta_directional_effects", width = 220, height = 200)

# ---- Figure S4. Ontogenetic filtering in relevant assemblages ---------------
ontogenetic <- raw_contrasts %>%
  filter(scenario == "adult_only_strict") %>%
  transmute(
    assemblage_id, assemblage_label,
    gamma_change_pct,
    mean_q0_change_pct,
    bray_mean_change_pct,
    sorensen_mean_change_pct,
    gdm_bray_relative_change_pct =
      if_else(
        is.finite(gdm_bray_explained_baseline) & abs(gdm_bray_explained_baseline) > 1e-8,
        100 * gdm_bray_delta_explained / gdm_bray_explained_baseline,
        NA_real_
      )
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label),
    names_to = "response_key", values_to = "change_pct"
  ) %>%
  mutate(
    response = recode(
      response_key,
      gamma_change_pct = "Gamma richness",
      mean_q0_change_pct = "Mean q0 richness",
      bray_mean_change_pct = "Mean Bray-Curtis",
      sorensen_mean_change_pct = "Mean Sørensen",
      gdm_bray_relative_change_pct = "GDM Bray explained deviance"
    )
  ) %>%
  summarise_distribution(
    value_col = "change_pct",
    group_cols = c("assemblage_id", "assemblage_label", "response")
  )

if (nrow(ontogenetic) > 0L) {
  p_s4 <- ggplot(ontogenetic, aes(x = response, y = median, colour = assemblage_label)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey35") +
    geom_errorbar(aes(ymin = p10, ymax = p90), width = 0, linewidth = 0.4) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
    labs(
      title = "Sensitivity to excluding explicitly non-adult individuals",
      subtitle = "Applied only where a strict adult-only scenario was biologically meaningful.",
      x = NULL, y = "Change relative to the appropriate baseline (%)"
    ) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_supplement(p_s4, "FigS4_ontogenetic_filtering", width = 170, height = 100)
}

# ---- Figure S5. Detailed GDM sensitivity ------------------------------------
gdm_full <- raw_contrasts %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label,
    bray_relative_explained_change = if_else(
      is.finite(gdm_bray_explained_baseline) & abs(gdm_bray_explained_baseline) > 1e-8,
      100 * gdm_bray_delta_explained / gdm_bray_explained_baseline,
      NA_real_
    ),
    sorensen_relative_explained_change = if_else(
      is.finite(gdm_sorensen_explained_baseline) & abs(gdm_sorensen_explained_baseline) > 1e-8,
      100 * gdm_sorensen_delta_explained / gdm_sorensen_explained_baseline,
      NA_real_
    ),
    bray_prediction_loss = 1 - gdm_bray_predicted_stability,
    sorensen_prediction_loss = 1 - gdm_sorensen_predicted_stability
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label, scenario, scenario_label),
    names_to = "response_key", values_to = "value"
  ) %>%
  mutate(
    response_type = case_when(
      str_detect(response_key, "relative_explained_change") ~ "Relative change in explained deviance (%)",
      TRUE ~ "Loss of predicted-dissimilarity stability"
    ),
    distance = case_when(
      str_detect(response_key, "^bray") ~ "Bray-Curtis",
      TRUE ~ "Sørensen"
    )
  ) %>%
  summarise_distribution(
    value_col = "value",
    group_cols = c(
      "assemblage_id", "assemblage_label", "scenario", "scenario_label",
      "response_type", "distance"
    )
  ) %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[match(ALL_SCENARIOS_ORDERED, scenario_info$scenario)]
    )
  )

p_s5 <- ggplot(gdm_full, aes(x = scenario_label, y = median, colour = assemblage_label)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
  geom_errorbar(aes(ymin = p10, ymax = p90), width = 0, linewidth = 0.25, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  facet_grid(response_type ~ distance, scales = "free_y") +
  scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
  labs(
    title = "Full GDM sensitivity analysis",
    x = NULL, y = NULL
  ) +
  theme_paper(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

save_supplement(p_s5, "FigS5_full_GDM_sensitivity", width = 220, height = 170)

# ---- Optional Figure S6. Predictor-weight changes, only if extraction worked
if (
  nrow(gdm_terms) > 0L &&
  all(c("distance", "predictor", "baseline_ispline_weight_median",
        "scenario_ispline_weight_median") %in% names(gdm_terms))
) {
  gdm_weights <- gdm_terms %>%
    filter(
      assemblage_id %in% MAIN_ASSEMBLAGES,
      scenario %in% MAIN_SCENARIOS
    ) %>%
    group_by(assemblage_id, scenario, distance) %>%
    mutate(
      baseline_weight_share = baseline_ispline_weight_median / sum(baseline_ispline_weight_median, na.rm = TRUE),
      scenario_weight_share = scenario_ispline_weight_median / sum(scenario_ispline_weight_median, na.rm = TRUE),
      delta_weight_share = 100 * (scenario_weight_share - baseline_weight_share)
    ) %>%
    ungroup() %>%
    left_join(scenario_info %>% select(scenario, scenario_label), by = "scenario") %>%
    left_join(metadata %>% select(assemblage_id, assemblage_label), by = "assemblage_id")
  
  p_s6 <- ggplot(
    gdm_weights,
    aes(x = predictor, y = delta_weight_share, colour = assemblage_label)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.6, alpha = 0.8) +
    facet_grid(distance ~ scenario_label, scales = "free_y") +
    scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
    labs(
      title = "Changes in relative GDM predictor contributions",
      subtitle = "Differences in I-spline coefficient shares between each scenario and its baseline.",
      x = NULL, y = "Change in relative predictor contribution (percentage points)"
    ) +
    theme_paper(base_size = 7) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_supplement(p_s6, "FigS6_GDM_predictor_weight_changes", width = 250, height = 150)
} else {
  message(
    "No usable GDM term-weight table found: FigS6 was not produced. ",
    "This does not affect the other figures."
  )
}

# ---- 4. Export compact data tables used by the figures ----------------------
readr::write_csv(heatmap_grid, file.path(FIG_DIR, "figure_data_Fig2.csv"))
readr::write_csv(workflow_metrics, file.path(FIG_DIR, "figure_data_Fig3.csv"))
readr::write_csv(gdm_explained, file.path(FIG_DIR, "figure_data_Fig4_explained_deviance.csv"))
readr::write_csv(gdm_prediction, file.path(FIG_DIR, "figure_data_Fig4_prediction_stability.csv"))

message("\nPublication figures completed.")
message("Main-text figures: ", file.path(FIG_DIR, "main_text"))
message("Supplementary figures: ", file.path(FIG_DIR, "supplement"))
message("\nSuggested main-text sequence:")
message("  Fig. 1  Taxonomic and ontogenetic context")
message("  Fig. 2  Cross-taxon robustness heatmap")
message("  Fig. 3  Directional workflow effects")
message("  Fig. 4  GDM ecological inference")
message("\nA Blowes alpha-gamma-occupancy synthesis can be added later as a final conceptual figure.")

