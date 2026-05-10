# ============================================================
# 14_sensitivity_to_min_sites.R
# Sensitivity of main conclusions to species inclusion threshold
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(broom.mixed)
  library(car)
  library(readr)
})

# ------------------------------------------------------------
# 0) Parameters
# ------------------------------------------------------------

min_sites_vec <- 5:10
min_sites_ref <- 7

out_dir_base <- file.path(
  "figures",
  "detectability",
  "carabidae",
  "sensitivity_min_sites"
)

dir.create(out_dir_base, recursive = TRUE, showWarnings = FALSE)

subdirs <- c(
  "occupancy",
  "detectability",
  "trait",
  "tables",
  "species_selection",
  "occ_model_output",
  "sp_protocol_favored"
)

for (k in min_sites_vec) {
  base_k <- file.path(out_dir_base, paste0("min_sites_", k))
  dir.create(base_k, recursive = TRUE, showWarnings = FALSE)
  
  for (sd in subdirs) {
    dir.create(file.path(base_k, sd), recursive = TRUE, showWarnings = FALSE)
  }
}

# ------------------------------------------------------------
# 1) Run pipeline for each threshold
# ------------------------------------------------------------

res_sens <- vector("list", length(min_sites_vec))
names(res_sens) <- paste0("k", min_sites_vec)

for (threshold in min_sites_vec) {
  
  message("Running sensitivity analysis for min_sites = ", threshold)
  
  env_k <- new.env(parent = globalenv())
  
  env_k$min_sites <- threshold
  env_k$res_sens <- res_sens
  env_k$out_dir <- file.path(out_dir_base, paste0("min_sites_", threshold))
  
  dir.create(env_k$out_dir, recursive = TRUE, showWarnings = FALSE)
  
  subdirs <- c(
    "species_selection",
    "occ_model_output",
    "occupancy",
    "detectability",
    "traits",
    "trait",
    "tables",
    "figures"
  )
  
  for (sd in subdirs) {
    dir.create(file.path(env_k$out_dir, sd), recursive = TRUE, showWarnings = FALSE)
  }
  
  source("analyses/Ecol Indic - Detection/04_detection_effort_reviewed.R", local = env_k)
  source("analyses/Ecol Indic - Detection/05_occupancy_fit_reviewed.R", local = env_k)
  source("analyses/Ecol Indic - Detection/06_occupancy_figs.R", local = env_k)
  source("analyses/Ecol Indic - Detection/07_occupancy_stat.R", local = env_k)
  source("analyses/Ecol Indic - Detection/09_traits_detectability.R", local = env_k)
  
  res_sens[[paste0("k", threshold)]] <- list(
    k = threshold,
    n_species = dplyr::n_distinct(env_k$tab_p$species),
    tab_p = env_k$tab_p,
    tab_beta = env_k$tab_beta,
    m_bl_add = env_k$m_bl_add,
    m_bl_int = env_k$m_bl_int,
    m_wing_add = env_k$m_wing_add,
    m_wing_int = env_k$m_wing_int
  )
}

saveRDS(
  res_sens,
  file.path(out_dir_base, "sensitivity_min_sites_results.rds")
)

# ------------------------------------------------------------
# 2) Number of species retained
# ------------------------------------------------------------

n_species_df <- imap_dfr(res_sens, \(x, nm) {
  tibble(
    k = x$k,
    n_species = x$n_species
  )
})

write_csv(
  n_species_df,
  file.path(out_dir_base, "tables/n_species_by_min_sites.csv")
)

p_nsp <- ggplot(n_species_df, aes(x = k, y = n_species)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    x = "Minimum number of occupied sites",
    y = "Number of species retained",
    title = "Species retained across inclusion thresholds"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir_base, "Annex_S5_n_species_by_min_sites.png"),
  p_nsp,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 3) Stability of protocol ranking
# ------------------------------------------------------------

method_order <- c(
  "Pitfall2", "Pitfall4", "Pitfall6",
  "Pitfall8", "Pitfall10", "GPD"
)

med_p <- imap(res_sens, \(x, nm) {
  
  x$tab_p %>%
    filter(method %in% method_order) %>%
    mutate(method = factor(method, levels = method_order)) %>%
    group_by(method) %>%
    summarise(
      p_med = median(p_hat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(method)
})

med_ref <- med_p[[paste0("k", min_sites_ref)]]

stab_rank <- imap_dfr(med_p, \(df, nm) {
  
  df2 <- df %>%
    left_join(
      med_ref %>% rename(p_med_ref = p_med),
      by = "method"
    )
  
  tibble(
    k = as.integer(sub("k", "", nm)),
    rho = cor(
      df2$p_med,
      df2$p_med_ref,
      method = "spearman",
      use = "complete.obs"
    )
  )
})

write_csv(
  stab_rank,
  file.path(out_dir_base, "tables/stability_protocol_ranking.csv")
)

p_rank <- ggplot(stab_rank, aes(x = k, y = rho)) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Minimum number of occupied sites",
    y = "Spearman rank correlation with k = 7",
    title = "Stability of protocol ranking"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir_base, "Annex_S5_protocol_ranking_stability.png"),
  p_rank,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 4) Stability of effort effect
# ------------------------------------------------------------

eff_slope_df <- imap_dfr(res_sens, \(x, nm) {
  
  x$tab_beta %>%
    filter(component == "det", term == "eff_z") %>%
    transmute(
      k = x$k,
      estimate = estimate,
      se = se,
      lcl = estimate - 1.96 * se,
      ucl = estimate + 1.96 * se
    )
})

write_csv(
  eff_slope_df,
  file.path(out_dir_base, "tables/stability_effort_effect.csv")
)

p_eff <- ggplot(eff_slope_df, aes(x = estimate, y = factor(k))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), height = 0.2) +
  geom_point(size = 2) +
  labs(
    x = "Effect of sampling effort on detectability (βeff_z)",
    y = "Minimum number of occupied sites",
    title = "Stability of effort effect"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir_base, "Annex_S5_effort_effect_stability.png"),
  p_eff,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 5) Stability of body-size effect
# ------------------------------------------------------------

bl_add_df <- imap_dfr(res_sens, \(x, nm) {
  
  broom.mixed::tidy(x$m_bl_add, effects = "fixed") %>%
    filter(term == "log_BL") %>%
    transmute(
      k = x$k,
      estimate,
      se = std.error,
      lcl = estimate - 1.96 * std.error,
      ucl = estimate + 1.96 * std.error
    )
})

write_csv(
  bl_add_df,
  file.path(out_dir_base, "tables/stability_body_size_additive_effect.csv")
)

p_bl <- ggplot(bl_add_df, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    x = "Minimum number of occupied sites",
    y = "Body-size effect on detectability (βlogBL)",
    title = "Stability of body-size effect"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir_base, "Annex_S5_body_size_effect_stability.png"),
  p_bl,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 6) Stability of method × body-size interaction
# ------------------------------------------------------------

bl_int_anova <- imap_dfr(res_sens, \(x, nm) {
  
  a <- as.data.frame(car::Anova(x$m_bl_int, type = 3))
  a$term <- rownames(a)
  
  a %>%
    filter(term == "method:log_BL") %>%
    transmute(
      k = x$k,
      chisq = Chisq,
      df = Df,
      p_value = `Pr(>Chisq)`
    )
})

write_csv(
  bl_int_anova,
  file.path(out_dir_base, "tables/stability_method_body_size_interaction.csv")
)

p_bl_int <- ggplot(bl_int_anova, aes(x = k, y = p_value)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_line() +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Minimum number of occupied sites",
    y = "P-value for method × body size",
    title = "Stability of method × body-size interaction"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir_base, "Annex_S5_method_body_size_interaction_stability.png"),
  p_bl_int,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 7) Stability of wing-development effects
# ------------------------------------------------------------

wing_add_df <- imap_dfr(res_sens, \(x, nm) {
  
  broom.mixed::tidy(x$m_wing_add, effects = "fixed") %>%
    filter(term == "full_wingnot_fully_winged") %>%
    transmute(
      k = x$k,
      estimate,
      se = std.error,
      lcl = estimate - 1.96 * std.error,
      ucl = estimate + 1.96 * std.error
    )
})

write_csv(
  wing_add_df,
  file.path(out_dir_base, "tables/stability_wing_additive_effect.csv")
)

p_wing_add <- ggplot(wing_add_df, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    x = "Minimum number of occupied sites",
    y = "Effect of not fully winged vs fully winged",
    title = "Stability of wing-development additive effect"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir_base, "Annex_S5_wing_additive_effect_stability.png"),
  p_wing_add,
  width = 7,
  height = 5,
  dpi = 300
)

wing_int_anova <- imap_dfr(res_sens, \(x, nm) {
  
  a <- as.data.frame(car::Anova(x$m_wing_int, type = 3))
  a$term <- rownames(a)
  
  a %>%
    filter(term == "method:full_wing") %>%
    transmute(
      k = x$k,
      chisq = Chisq,
      df = Df,
      p_value = `Pr(>Chisq)`
    )
})

write_csv(
  wing_int_anova,
  file.path(out_dir_base, "tables/stability_method_wing_interaction.csv")
)

p_wing_int <- ggplot(wing_int_anova, aes(x = k, y = p_value)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_line() +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Minimum number of occupied sites",
    y = "P-value for method × wing development",
    title = "Stability of method × wing-development interaction"
  ) +
  theme_minimal()

ggsave(
  file.path(out_dir_base, "Annex_S5_method_wing_interaction_stability.png"),
  p_wing_int,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 8) Summary table for manuscript / appendix
# ------------------------------------------------------------

summary_sensitivity <- n_species_df %>%
  left_join(stab_rank, by = "k") %>%
  left_join(
    eff_slope_df %>%
      select(k, effort_estimate = estimate, effort_lcl = lcl, effort_ucl = ucl),
    by = "k"
  ) %>%
  left_join(
    bl_add_df %>%
      select(k, body_size_estimate = estimate, body_size_lcl = lcl, body_size_ucl = ucl),
    by = "k"
  ) %>%
  left_join(
    bl_int_anova %>%
      select(k, p_method_body_size = p_value),
    by = "k"
  ) %>%
  left_join(
    wing_add_df %>%
      select(k, wing_estimate = estimate, wing_lcl = lcl, wing_ucl = ucl),
    by = "k"
  ) %>%
  left_join(
    wing_int_anova %>%
      select(k, p_method_wing = p_value),
    by = "k"
  )

write_csv(
  summary_sensitivity,
  file.path(out_dir_base, "tables/summary_sensitivity_min_sites.csv")
)

# ------------------------------------------------------------
# 9) Compact combined figure
# ------------------------------------------------------------

p_combined <- cowplot::plot_grid(
  p_nsp,
  p_rank,
  p_eff,
  p_bl,
  p_bl_int,
  p_wing_int,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave(
  file.path(out_dir_base, "Annex_S5_sensitivity_min_sites_combined.png"),
  p_combined,
  width = 12,
  height = 12,
  dpi = 300
)

###############

library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

species_focus <- c(
  "Abax parallelepipedus",
  "Carabus nemoralis",
  "Pseudoophonus rufipes",
  "Bembidion obtusum", 
  "Poecilus cupreus"
)

method_order <- c(
  "Pitfall2", "Pitfall4", "Pitfall6",
  "Pitfall8", "Pitfall10", "GPD"
)

p_hat_sens_species <- purrr::imap_dfr(res_sens, function(x, nm) {
  
  x$tab_p %>%
    filter(species %in% species_focus) %>%
    mutate(
      min_sites = x$k
    )
}) %>%
  mutate(
    method = factor(method, levels = method_order),
    species = factor(species, levels = species_focus)
  )

p_species_stability <- ggplot(
  p_hat_sens_species,
  aes(
    x = min_sites,
    y = p_hat,
    colour = method,
    group = method
  )
) +
  geom_ribbon(
    aes(ymin = lcl, ymax = ucl, fill = method),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ species, ncol = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Minimum number of occupied sites",
    y = "Detection probability p̂",
    colour = "Protocol",
    fill = "Protocol"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom"
  )

p_species_stability

ggsave(
  file.path(out_dir_base, "Annex_S5_sensitivity_min_sites_species.png"),
  p_species_stability,
  width = 12,
  height = 12,
  dpi = 300
)
