# =============================================================================
# 03_make_multitaxa_blowes_figures.R
# Empirical alpha-gamma-occupancy ("Blowes") synthesis for RMQS 2024
# multi-taxon taxonomic-uncertainty scenarios.
#
# It reads results_by_iter_all.csv produced by:
#   01_run_multitaxa_uncertainty_2024_v1_5_gdm_dataframe_fix.R
#
# It does NOT rerun any scenarios or models.
#
# Key identity used here:
#   mean occupancy = mean local richness (alpha) / regional richness (gamma)
# when scenarios are evaluated on the same set of stations.
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

# ---- 1. Paths and choices ---------------------------------------------------
AUDIT_DIR  <- "outputs_multitaxa_2024"
RESULT_DIR <- file.path(AUDIT_DIR, "uncertainty_results")
FIG_DIR    <- file.path(RESULT_DIR, "publication_figures")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "main_text"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "supplement"), recursive = TRUE, showWarnings = FALSE)

MAIN_ASSEMBLAGES <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "carabidae_pitfall",
  "formicidae_pitfall",
  "isopoda_pitfall"
)

# These scenarios form the central empirical Blowes synthesis.
ERROR_SCENARIOS <- c(
  "species_observed_congeneric_5pct",
  "species_observed_congeneric_10pct",
  "species_observed_rare_weighted_10pct"
)

WORKFLOW_SCENARIOS <- c(
  "rtu_rare_species_to_genus",
  "rtu_species_only_drop_unresolved",
  "rtu_genus_level",
  "rtu_family_level"
)

STAGE_SCENARIO <- "adult_only_strict"

# ---- 2. Display labels and style --------------------------------------------
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
  ~scenario, ~scenario_label, ~scenario_class,
  "adult_only_strict", "Adults only", "Ontogenetic filtering",
  "species_observed_congeneric_5pct", "Observed-pool\nerror 5%", "Identification error",
  "species_observed_congeneric_10pct", "Observed-pool\nerror 10%", "Identification error",
  "species_observed_rare_weighted_10pct", "Rare-weighted\nerror 10%", "Identification error",
  "rtu_rare_species_to_genus", "Rare taxa\nreported at genus", "Workflow decision",
  "rtu_species_only_drop_unresolved", "Drop unresolved\nrecords", "Workflow decision",
  "rtu_genus_level", "Genus-level\nresolution", "Taxonomic coarsening",
  "rtu_family_level", "Family-level\nresolution", "Taxonomic coarsening"
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
      panel.grid.major = element_line(colour = "grey91", linewidth = 0.28),
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_main <- function(plot, filename, width, height) {
  ggsave(
    file.path(FIG_DIR, "main_text", paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "mm",
    device = grDevices::pdf
  )
  ggsave(
    file.path(FIG_DIR, "main_text", paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "mm",
    dpi = 450
  )
}

save_supplement <- function(plot, filename, width, height) {
  ggsave(
    file.path(FIG_DIR, "supplement", paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "mm",
    device = grDevices::pdf
  )
  ggsave(
    file.path(FIG_DIR, "supplement", paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "mm",
    dpi = 450
  )
}

# Exact coordinate transformation:
# alpha ratio = 1 + change in mean q0 / 100
# gamma ratio = 1 + change in gamma / 100
# mean occupancy ratio = alpha ratio / gamma ratio
#
# Calculating medians on ratios, then converting to % change, preserves the
# alpha-gamma relationship for the plotted point.
summarise_blowes <- function(data) {
  data %>%
    filter(
      is.finite(mean_q0_change_pct),
      is.finite(gamma_change_pct)
    ) %>%
    mutate(
      alpha_ratio = 1 + mean_q0_change_pct / 100,
      gamma_ratio = 1 + gamma_change_pct / 100
    ) %>%
    filter(alpha_ratio > 0, gamma_ratio > 0) %>%
    mutate(
      occupancy_ratio = alpha_ratio / gamma_ratio,
      delta_log_occupancy = log(occupancy_ratio)
    ) %>%
    group_by(
      assemblage_id, assemblage_label, scenario, scenario_label,
      scenario_class
    ) %>%
    summarise(
      alpha_ratio_median = median(alpha_ratio, na.rm = TRUE),
      alpha_ratio_p10 = quantile(alpha_ratio, probs = 0.10, na.rm = TRUE, names = FALSE),
      alpha_ratio_p90 = quantile(alpha_ratio, probs = 0.90, na.rm = TRUE, names = FALSE),
      gamma_ratio_median = median(gamma_ratio, na.rm = TRUE),
      gamma_ratio_p10 = quantile(gamma_ratio, probs = 0.10, na.rm = TRUE, names = FALSE),
      gamma_ratio_p90 = quantile(gamma_ratio, probs = 0.90, na.rm = TRUE, names = FALSE),
      occupancy_ratio_median = median(occupancy_ratio, na.rm = TRUE),
      occupancy_ratio_p10 = quantile(occupancy_ratio, probs = 0.10, na.rm = TRUE, names = FALSE),
      occupancy_ratio_p90 = quantile(occupancy_ratio, probs = 0.90, na.rm = TRUE, names = FALSE),
      delta_log_occupancy_median = median(delta_log_occupancy, na.rm = TRUE),
      n_iter = n(),
      .groups = "drop"
    ) %>%
    mutate(
      delta_alpha_pct = 100 * (alpha_ratio_median - 1),
      delta_alpha_p10 = 100 * (alpha_ratio_p10 - 1),
      delta_alpha_p90 = 100 * (alpha_ratio_p90 - 1),
      delta_gamma_pct = 100 * (gamma_ratio_median - 1),
      delta_gamma_p10 = 100 * (gamma_ratio_p10 - 1),
      delta_gamma_p90 = 100 * (gamma_ratio_p90 - 1),
      delta_occupancy_pct = 100 * (occupancy_ratio_median - 1),
      delta_occupancy_p10 = 100 * (occupancy_ratio_p10 - 1),
      delta_occupancy_p90 = 100 * (occupancy_ratio_p90 - 1),
      apparent_signature = case_when(
        occupancy_ratio_median > 1 + 1e-8 ~ "Apparent homogenisation",
        occupancy_ratio_median < 1 - 1e-8 ~ "Apparent differentiation",
        TRUE ~ "No occupancy change"
      )
    )
}

square_limits <- function(data, padding = 0.08, minimum_span = 20) {
  vals <- c(
    data$delta_alpha_p10, data$delta_alpha_p90,
    data$delta_gamma_p10, data$delta_gamma_p90,
    0
  )
  vals <- vals[is.finite(vals)]
  
  lo <- min(vals)
  hi <- max(vals)
  span <- max(hi - lo, minimum_span)
  
  lo <- lo - padding * span
  hi <- hi + padding * span
  
  # Round outwards to readable multiples of 10.
  c(
    floor(lo / 10) * 10,
    ceiling(hi / 10) * 10
  )
}

plot_blowes_space <- function(data, title, subtitle, ncol, limits) {
  ggplot(
    data,
    aes(
      x = delta_alpha_pct,
      y = delta_gamma_pct,
      colour = assemblage_label
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.32, colour = "grey45") +
    geom_vline(xintercept = 0, linewidth = 0.32, colour = "grey45") +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", linewidth = 0.55, colour = "black"
    ) +
    # Horizontal and vertical 10th–90th percentile ranges:
    geom_segment(
      aes(
        x = delta_alpha_p10, xend = delta_alpha_p90,
        y = delta_gamma_pct, yend = delta_gamma_pct
      ),
      linewidth = 0.38, alpha = 0.75, na.rm = TRUE
    ) +
    geom_segment(
      aes(
        x = delta_alpha_pct, xend = delta_alpha_pct,
        y = delta_gamma_p10, yend = delta_gamma_p90
      ),
      linewidth = 0.38, alpha = 0.75, na.rm = TRUE
    ) +
    geom_point(size = 2.5, na.rm = TRUE) +
    facet_wrap(~ scenario_label, ncol = ncol) +
    coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    scale_colour_manual(
      values = assemblage_colours,
      name = "Assemblage"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression("Change in mean local richness, " * Delta * alpha * " (%)"),
      y = expression("Change in regional richness, " * Delta * gamma * " (%)")
    ) +
    theme_paper(base_size = 8.5) +
    theme(
      legend.position = "bottom",
      panel.spacing = grid::unit(4, "mm")
    )
}

# ---- 3. Read outputs ---------------------------------------------------------
raw_path <- file.path(RESULT_DIR, "results_by_iter_all.csv")
manifest_path <- file.path(AUDIT_DIR, "assemblage_manifest.csv")

if (!file.exists(raw_path)) {
  stop("Cannot find: ", raw_path)
}
if (!file.exists(manifest_path)) {
  stop("Cannot find: ", manifest_path)
}

raw <- readr::read_csv(raw_path, show_col_types = FALSE)
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

metadata <- manifest %>%
  transmute(
    assemblage_id,
    assemblage_label = recode(
      assemblage_id,
      !!!assemblage_labels,
      .default = assemblage_id
    )
  )

raw_contrasts <- raw %>%
  filter(
    !is.na(baseline_scenario),
    scenario != baseline_scenario
  ) %>%
  left_join(metadata, by = "assemblage_id") %>%
  left_join(scenario_info, by = "scenario") %>%
  mutate(
    scenario_label = coalesce(scenario_label, scenario),
    scenario_class = coalesce(scenario_class, scenario_family),
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels)
    )
  )

if (!all(c("mean_q0_change_pct", "gamma_change_pct") %in% names(raw_contrasts))) {
  stop(
    "The result table needs mean_q0_change_pct and gamma_change_pct ",
    "to construct the alpha-gamma-occupancy synthesis."
  )
}

# ---- 4. Main text Figure 5 ---------------------------------------------------
workflow_blowes <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% WORKFLOW_SCENARIOS
  ) %>%
  summarise_blowes() %>%
  mutate(
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])
    ),
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[
        match(WORKFLOW_SCENARIOS, scenario_info$scenario)
      ]
    )
  )

error_blowes <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% ERROR_SCENARIOS
  ) %>%
  summarise_blowes() %>%
  mutate(
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])
    ),
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[
        match(ERROR_SCENARIOS, scenario_info$scenario)
      ]
    )
  )

if (nrow(workflow_blowes) == 0L || nrow(error_blowes) == 0L) {
  stop(
    "No valid Blowes coordinates were created. Check that scenario IDs in ",
    "WORKFLOW_SCENARIOS and ERROR_SCENARIOS are present in results_by_iter_all.csv."
  )
}

workflow_limits <- square_limits(workflow_blowes, minimum_span = 40)
error_limits <- square_limits(error_blowes, minimum_span = 35)

p_workflow <- plot_blowes_space(
  workflow_blowes,
  title = "Workflow decisions and taxonomic coarsening",
  subtitle = "Dashed line: unchanged mean occupancy. Above = apparent differentiation; below = apparent homogenisation.",
  ncol = 2,
  limits = workflow_limits
)

p_errors <- plot_blowes_space(
  error_blowes,
  title = "Moderate identification error within the observed species pool",
  subtitle = "Axes are rescaled to resolve the smaller changes produced by identification-error scenarios.",
  ncol = 3,
  limits = error_limits
)

p_fig5 <- (p_workflow / p_errors) +
  plot_layout(heights = c(1.28, 1)) +
  plot_annotation(
    title = "Taxonomic uncertainty can mimic opposing beta-diversity signatures",
    subtitle = "Each point is the scenario median for one assemblage; crosshairs show the 10th–90th percentile range across stochastic iterations.",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_main(p_fig5, "Fig5_multitaxa_Blowes_synthesis", width = 190, height = 225)

# ---- 5. Supplementary Figure S6: derived occupancy changes -----------------
all_blowes <- raw_contrasts %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    scenario %in% c(ERROR_SCENARIOS, WORKFLOW_SCENARIOS)
  ) %>%
  summarise_blowes() %>%
  mutate(
    assemblage_label = factor(
      assemblage_label,
      levels = rev(unname(assemblage_labels[MAIN_ASSEMBLAGES]))
    ),
    scenario_label = factor(
      scenario_label,
      levels = scenario_info$scenario_label[
        match(c(ERROR_SCENARIOS, WORKFLOW_SCENARIOS), scenario_info$scenario)
      ]
    )
  )

p_s6 <- ggplot(
  all_blowes,
  aes(x = scenario_label, y = assemblage_label, fill = delta_occupancy_pct)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%+.0f", delta_occupancy_pct)),
    size = 2.5
  ) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    name = expression(Delta * " mean occupancy (%)"),
    na.value = "grey86"
  ) +
  labs(
    title = "Derived changes in mean taxon occupancy",
    subtitle = "Positive values indicate apparent homogenisation; negative values indicate apparent differentiation.",
    x = NULL,
    y = NULL
  ) +
  theme_paper(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid = element_blank()
  )

save_supplement(p_s6, "FigS6_multitaxa_Blowes_occupancy_heatmap", width = 185, height = 112)

# ---- 6. Supplementary Figure S7: adult-only / ontogenetic filtering --------
stage_blowes <- raw_contrasts %>%
  filter(scenario == STAGE_SCENARIO) %>%
  summarise_blowes() %>%
  mutate(
    assemblage_label = factor(
      assemblage_label,
      levels = unname(assemblage_labels)
    )
  )

if (nrow(stage_blowes) > 0L) {
  stage_limits <- square_limits(stage_blowes, minimum_span = 30)
  
  p_s7 <- ggplot(
    stage_blowes,
    aes(
      x = delta_alpha_pct,
      y = delta_gamma_pct,
      colour = assemblage_label
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.32, colour = "grey45") +
    geom_vline(xintercept = 0, linewidth = 0.32, colour = "grey45") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.55) +
    geom_segment(
      aes(
        x = delta_alpha_p10, xend = delta_alpha_p90,
        y = delta_gamma_pct, yend = delta_gamma_pct
      ),
      linewidth = 0.4, na.rm = TRUE
    ) +
    geom_segment(
      aes(
        x = delta_alpha_pct, xend = delta_alpha_pct,
        y = delta_gamma_p10, yend = delta_gamma_p90
      ),
      linewidth = 0.4, na.rm = TRUE
    ) +
    geom_point(size = 2.8) +
    coord_equal(xlim = stage_limits, ylim = stage_limits, expand = FALSE) +
    scale_colour_manual(values = assemblage_colours, name = "Assemblage") +
    labs(
      title = "Ontogenetic filtering in alpha-gamma space",
      subtitle = "Adults-only scenario, shown only for assemblages where explicitly non-adult records were present.",
      x = expression("Change in mean local richness, " * Delta * alpha * " (%)"),
      y = expression("Change in regional richness, " * Delta * gamma * " (%)")
    ) +
    theme_paper()
  
  save_supplement(p_s7, "FigS7_Blowes_ontogenetic_filtering", width = 140, height = 115)
}

# ---- 7. Export data tables ---------------------------------------------------
readr::write_csv(workflow_blowes, file.path(FIG_DIR, "Blowes_workflow_summary.csv"))
readr::write_csv(error_blowes, file.path(FIG_DIR, "Blowes_identification_error_summary.csv"))
readr::write_csv(all_blowes, file.path(FIG_DIR, "Blowes_all_main_scenarios_summary.csv"))
if (nrow(stage_blowes) > 0L) {
  readr::write_csv(stage_blowes, file.path(FIG_DIR, "Blowes_ontogenetic_filtering_summary.csv"))
}

message("\nBlowes figures completed.")
message("Main text: ", file.path(FIG_DIR, "main_text", "Fig5_multitaxa_Blowes_synthesis.pdf"))
message("Supplement: ", file.path(FIG_DIR, "supplement"))
message("\nImportant interpretation:")
message("  The alpha-gamma relationship is empirical and scenario-specific.")
message("  It does not represent temporal biodiversity change.")
message("  It visualises apparent homogenisation or differentiation generated")
message("  by taxonomic uncertainty in the same set of communities.")

