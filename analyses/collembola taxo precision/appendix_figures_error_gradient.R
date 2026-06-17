# ============================================================
# appendix_figures_error_gradient — figures de l'annexe gradient d'erreur
# ============================================================

.workflow_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
WORKFLOW_DIR <- if (!is.na(.workflow_file)) dirname(.workflow_file) else getwd()
source(file.path(WORKFLOW_DIR, "R", "00_config.R"))
source(file.path(WORKFLOW_DIR, "R", "functions_taxonomic_uncertainty.R"))

load_required_step("09_appendix_error_gradient_analysis")

appendix_fig_dir <- file.path(PUB_DIR, "appendix_error_gradient")
dir.create(appendix_fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# FIGURE A1 — courbes de réponse des indices de biodiversité
# -----------------------------
plot_metrics <- c("mean_local_q0", "mean_local_q1", "mean_local_q2", "mean_coverage_chao", "gamma_taxon_units")

p_appendix_biodiv <- appendix_biodiversity_curve_summary %>%
  filter(metric %in% plot_metrics) %>%
  mutate(
    metric_label = factor(metric_label,
      levels = c("Alpha q0 richness", "Alpha Hill q1", "Alpha Hill q2", "Mean sample coverage", "Gamma richness")
    )
  ) %>%
  ggplot(aes(x = error_pct, y = mean_percent_change)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = c(-100 * ID_ERROR_APPENDIX_CRITICAL_REL_CHANGE,
                            100 * ID_ERROR_APPENDIX_CRITICAL_REL_CHANGE),
             linewidth = 0.25, linetype = "dotted") +
  geom_ribbon(aes(ymin = p10_percent_change, ymax = p90_percent_change), alpha = 0.18) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Response of biodiversity indices to increasing identification error",
    subtitle = paste0(
      "Same-genus reassignment pool: ", appendix_pool_name,
      " | Dotted lines = ±", round(100 * ID_ERROR_APPENDIX_CRITICAL_REL_CHANGE), "% relative change"
    ),
    x = "Identification error rate (%)",
    y = "Relative change vs species baseline (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(appendix_fig_dir, "FigA1_biodiversity_response_curves.png"), p_appendix_biodiv,
       width = 9.5, height = 7.2, dpi = 600)
ggsave(file.path(appendix_fig_dir, "FigA1_biodiversity_response_curves.pdf"), p_appendix_biodiv,
       width = 9.5, height = 7.2)

# -----------------------------
# FIGURE A2 — courbes de stabilité
# -----------------------------
plot_stability_levels <- c(
  "Alpha q0 richness", "Alpha Hill q1", "Alpha Hill q2", "Mean sample coverage",
  "Bray-Curtis", "Jaccard", "Sørensen", "Ordination (Procrustes R2)"
)

p_appendix_stability <- appendix_stability_curve_summary %>%
  mutate(metric_label = factor(metric_label, levels = plot_stability_levels)) %>%
  ggplot(aes(x = error_pct, y = mean_stability)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = ID_ERROR_APPENDIX_CRITICAL_STABILITY, linewidth = 0.25, linetype = "dotted") +
  geom_ribbon(aes(ymin = p10_stability, ymax = p90_stability), alpha = 0.18) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Stability of ecological inference along the identification-error gradient",
    subtitle = paste0(
      "Dotted line = critical stability threshold (", ID_ERROR_APPENDIX_CRITICAL_STABILITY, ")"
    ),
    x = "Identification error rate (%)",
    y = "Mean stability relative to species baseline"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(appendix_fig_dir, "FigA2_stability_response_curves.png"), p_appendix_stability,
       width = 9.5, height = 9.0, dpi = 600)
ggsave(file.path(appendix_fig_dir, "FigA2_stability_response_curves.pdf"), p_appendix_stability,
       width = 9.5, height = 9.0)

# -----------------------------
# FIGURE A3 — seuil critique par métrique
# -----------------------------
appendix_threshold_plot_data <- appendix_critical_thresholds %>%
  mutate(
    threshold_label = recode(
      threshold_type,
      relative_change = "Critical from relative change",
      stability = "Critical from stability"
    ),
    metric_label = forcats::fct_reorder(metric_label, critical_error_pct, .na_rm = FALSE, .desc = TRUE)
  )

p_appendix_thresholds <- appendix_threshold_plot_data %>%
  ggplot(aes(x = critical_error_pct, y = metric_label)) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ threshold_label, scales = "free_x") +
  labs(
    title = "Estimated critical error rate by metric",
    subtitle = "First error level at which the chosen criticality criterion is reached",
    x = "Critical error rate (%)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(appendix_fig_dir, "FigA3_critical_error_thresholds.png"), p_appendix_thresholds,
       width = 9.0, height = 5.0, dpi = 600)
ggsave(file.path(appendix_fig_dir, "FigA3_critical_error_thresholds.pdf"), p_appendix_thresholds,
       width = 9.0, height = 5.0)
