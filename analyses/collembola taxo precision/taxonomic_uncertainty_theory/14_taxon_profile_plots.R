# ============================================================
# 14_taxon_profile_plots.R
# Figures for empirical taxon profiles and cross-taxon simulations
# ============================================================

plot_taxon_profile_axes <- function(profile_summaries) {
  dat <- profile_axes_long(profile_summaries)

  ggplot(dat, aes(x = taxon, y = value)) +
    geom_point(size = 3) +
    facet_wrap(~ axis, scales = "free_y", nrow = 1) +
    labs(
      title = "Taxonomic and sampling profiles of empirical groups",
      subtitle = "These descriptors position groups in the theoretical uncertainty framework.",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "grey90", colour = NA),
      panel.grid.minor = element_blank()
    )
}

plot_cross_taxon_scenarios <- function(
  blowes_summary,
  response = c("delta_gamma", "delta_occupancy"),
  selected_scenarios = c(
    "id_observed_uniform_10pct",
    "id_regional_uniform_10pct",
    "id_observed_rare_weighted_10pct",
    "workflow_mixed_rtu_10pct",
    "workflow_drop_unresolved_10pct",
    "workflow_genus_level",
    "workflow_family_level"
  )
) {
  response <- match.arg(response)
  med_col <- paste0(response, "_med")
  p10_col <- paste0(response, "_p10")
  p90_col <- paste0(response, "_p90")

  dat <- blowes_summary %>%
    filter(scenario %in% selected_scenarios) %>%
    mutate(scenario_label = factor(scenario_label, levels = unique(scenario_label)))

  if (nrow(dat) == 0) stop("No selected scenarios found in `blowes_summary`.", call. = FALSE)

  y_lab <- if (response == "delta_gamma") {
    "Change in regional richness, Δγ (%)"
  } else {
    "Change in mean occupancy (%)"
  }

  ggplot(dat, aes(x = scenario_label, y = .data[[med_col]], colour = taxon, group = taxon)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_errorbar(aes(ymin = .data[[p10_col]], ymax = .data[[p90_col]]), width = 0.15) +
    geom_point(size = 2.5) +
    coord_flip() +
    labs(
      title = if (response == "delta_gamma") {
        "Taxon-specific theoretical vulnerability: regional richness"
      } else {
        "Taxon-specific theoretical vulnerability: mean occupancy"
      },
      subtitle = "Synthetic communities use each group's real regional taxonomic hierarchy and sampling descriptors.",
      x = NULL,
      y = y_lab,
      colour = "Taxon"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

plot_cross_taxon_blowes <- function(
  blowes_summary,
  selected_scenarios = c(
    "id_observed_uniform_10pct",
    "id_regional_uniform_10pct",
    "workflow_mixed_rtu_10pct",
    "workflow_drop_unresolved_10pct",
    "workflow_genus_level",
    "workflow_family_level"
  )
) {
  dat <- blowes_summary %>% filter(scenario %in% selected_scenarios)
  if (nrow(dat) == 0) stop("No selected scenarios found in `blowes_summary`.", call. = FALSE)

  ggplot(dat, aes(x = delta_alpha_med, y = delta_gamma_med, colour = taxon, shape = scenario_label)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 3) +
    labs(
      title = "Cross-taxon theoretical responses in alpha–gamma space",
      subtitle = "Dashed line: unchanged mean occupancy. Each point is a scenario median across simulated communities.",
      x = "Change in mean local richness, Δα (%)",
      y = "Change in regional richness, Δγ (%)",
      colour = "Taxon",
      shape = "Scenario"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}
