# ============================================================
# 08_plotting.R
# Theory and empirical overlays
# ============================================================

plot_theory_blowes <- function(
  theory_summary,
  scenario,
  empirical_points = NULL,
  point_alpha = 0.25
) {
  required <- c("scenario", "delta_alpha_med", "delta_gamma_med", "beta_signature")
  missing <- setdiff(required, names(theory_summary))
  if (length(missing) > 0) {
    stop("`theory_summary` lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  dat <- theory_summary %>% filter(scenario == !!scenario)
  if (nrow(dat) == 0) stop("No theoretical rows found for this scenario.", call. = FALSE)

  p <- ggplot(dat, aes(x = delta_alpha_med, y = delta_gamma_med)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = "dashed") +
    geom_point(aes(colour = beta_signature), alpha = point_alpha, size = 1.8) +
    coord_equal() +
    labs(
      x = expression(paste(Delta, alpha, " (%)")),
      y = expression(paste(Delta, gamma, " (%)")),
      colour = "Apparent beta signature",
      title = paste("Theoretical alpha-gamma space:", scenario),
      subtitle = "Dashed line: unchanged mean occupancy"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  if (!is.null(empirical_points)) {
    needed <- c("delta_alpha_pct", "delta_gamma_pct")
    if (!all(needed %in% names(empirical_points))) {
      stop("`empirical_points` must contain delta_alpha_pct and delta_gamma_pct.", call. = FALSE)
    }

    p <- p +
      geom_point(
        data = empirical_points,
        aes(x = delta_alpha_pct, y = delta_gamma_pct),
        inherit.aes = FALSE,
        shape = 21,
        fill = "white",
        colour = "black",
        stroke = 0.9,
        size = 3.4
      )
  }

  p
}

plot_effect_against_architecture <- function(
  theory_summary,
  scenario,
  x = "cv_species_per_genus",
  y = "delta_occupancy_med"
) {
  needed <- c("scenario", x, y)
  if (!all(needed %in% names(theory_summary))) {
    stop("Requested variables are absent from `theory_summary`.", call. = FALSE)
  }

  theory_summary %>%
    filter(scenario == !!scenario) %>%
    ggplot(aes(x = .data[[x]], y = .data[[y]])) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(alpha = 0.35) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(
      x = x,
      y = y,
      title = paste("Theoretical sensitivity under:", scenario)
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}
