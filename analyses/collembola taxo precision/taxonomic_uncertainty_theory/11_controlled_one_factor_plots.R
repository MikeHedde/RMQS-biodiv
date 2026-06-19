# ============================================================
# 11_controlled_one_factor_plots.R
# Figures for the controlled OFAT experiment
# ============================================================

controlled_response_plot <- function(
  blowes_summary,
  scenario_domain = "identification_error",
  response = c("delta_gamma", "delta_occupancy"),
  title = NULL
) {
  response <- match.arg(response)

  names_needed <- c(
    "scenario_domain", "scenario_label", "factor_label", "factor_value",
    "intensity_label", "requested_intensity"
  )
  if (!all(names_needed %in% names(blowes_summary))) {
    stop("`blowes_summary` is missing required columns.", call. = FALSE)
  }

  response_cols <- switch(
    response,
    delta_gamma = c("delta_gamma_med", "delta_gamma_p10", "delta_gamma_p90"),
    delta_occupancy = c("delta_occupancy_med", "delta_occupancy_p10", "delta_occupancy_p90")
  )
  y_label <- switch(
    response,
    delta_gamma = expression(paste(Delta, gamma, " (%)")),
    delta_occupancy = expression(paste(Delta, " mean occupancy (%)"))
  )

  dat <- blowes_summary %>%
    filter(scenario_domain == !!scenario_domain) %>%
    mutate(intensity_label = factor(intensity_label, levels = unique(intensity_label[order(requested_intensity, na.last = TRUE)])))

  if (nrow(dat) == 0) stop("No rows available for the requested scenario domain.", call. = FALSE)

  if (is.null(title)) {
    title <- switch(
      scenario_domain,
      identification_error = "Controlled effects of identification-error mechanisms",
      reporting_workflow = "Controlled effects of reporting workflows",
      fixed_workflow = "Controlled effects of fixed taxonomic workflows",
      "Controlled theoretical responses"
    )
  }

  ggplot(dat, aes(x = factor_value, y = .data[[response_cols[1]]])) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey55") +
    geom_linerange(
      aes(ymin = .data[[response_cols[2]]], ymax = .data[[response_cols[3]]], colour = intensity_label),
      alpha = 0.55
    ) +
    geom_line(aes(colour = intensity_label, group = intensity_label), linewidth = 0.65) +
    geom_point(aes(colour = intensity_label), size = 1.8) +
    facet_grid(scenario_label ~ factor_label, scales = "free_x") +
    labs(
      x = NULL,
      y = y_label,
      colour = if (scenario_domain == "fixed_workflow") "Workflow" else "Requested rate",
      title = title,
      subtitle = "Lines show medians; vertical ranges show the 10th–90th percentiles across stochastic replicates."
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    )
}

controlled_fixed_workflow_plot <- function(
  blowes_summary,
  response = c("delta_gamma", "delta_occupancy"),
  title = NULL
) {
  response <- match.arg(response)
  response_cols <- switch(
    response,
    delta_gamma = c("delta_gamma_med", "delta_gamma_p10", "delta_gamma_p90"),
    delta_occupancy = c("delta_occupancy_med", "delta_occupancy_p10", "delta_occupancy_p90")
  )
  y_label <- switch(
    response,
    delta_gamma = expression(paste(Delta, gamma, " (%)")),
    delta_occupancy = expression(paste(Delta, " mean occupancy (%)"))
  )

  dat <- blowes_summary %>% filter(scenario_domain == "fixed_workflow")
  if (nrow(dat) == 0) stop("No fixed-workflow rows were found.", call. = FALSE)

  if (is.null(title)) title <- "Controlled effects of fixed taxonomic workflows"

  ggplot(dat, aes(x = factor_value, y = .data[[response_cols[1]]], colour = scenario_label)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey55") +
    geom_linerange(aes(ymin = .data[[response_cols[2]]], ymax = .data[[response_cols[3]]]), alpha = 0.55) +
    geom_line(aes(group = scenario_label), linewidth = 0.65) +
    geom_point(size = 1.8) +
    facet_wrap(~factor_label, scales = "free_x", nrow = 1) +
    labs(
      x = NULL,
      y = y_label,
      colour = "Workflow",
      title = title,
      subtitle = "Lines show medians; vertical ranges show the 10th–90th percentiles across stochastic replicates."
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
}

controlled_endpoint_heatmap <- function(endpoint_contrasts, response = "delta_gamma_med") {
  dat <- endpoint_contrasts %>%
    filter(response == !!response) %>%
    mutate(
      contrast_label = paste0(
        formatC(low_factor_value, format = "fg", digits = 3), " → ",
        formatC(high_factor_value, format = "fg", digits = 3)
      )
    )

  if (nrow(dat) == 0) stop("No endpoint contrasts found for the requested response.", call. = FALSE)

  ggplot(dat, aes(x = factor_label, y = scenario_label, fill = high_minus_low)) +
    geom_tile() +
    facet_wrap(~intensity_label, scales = "free_y") +
    labs(
      x = NULL,
      y = NULL,
      fill = "High − low response",
      title = paste("Endpoint contrasts:", response),
      subtitle = "Positive values mean that increasing the controlled parameter increases the response."
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}
