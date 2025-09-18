# ============================================================
# FIGURES & SYNTHÈSES p̂ et β (ψ ~ ALT, DOY)
# Requiert: tab_p, tab_beta (créés au script 05)
# ============================================================

tab_p   <- readr::read_csv(file.path(out_dir, "p_hat_by_method_unmarked_full.csv"), show_col_types = FALSE)
tab_beta<- readr::read_csv(file.path(out_dir, "psi_coef_ALT_DOY_by_species.csv"), show_col_types = FALSE)

if (nrow(tab_p) > 0) {
  comm_sum <- tab_p %>% group_by(method) %>%
    summarise(p_med = median(p_hat, na.rm=TRUE), n_sp = n_distinct(species), .groups="drop")
  
  tab_p$method <- factor(tab_p$method, levels = unique(tab_p$method))
  p_fig <- ggplot(tab_p, aes(x = method, y = p_hat)) +
    geom_violin(fill = "grey92") +
    geom_point(aes(colour = species), position = position_jitter(width = 0.08), alpha = 0.5) +
    geom_point(data = comm_sum, aes(y = p_med), color = "red", size = 3) +
    labs(y = "p̂ (détectabilité, occupancy; effort médian par méthode)",
         x = NULL, title = "Détectabilité par méthode — ajustée pour l’effort") +
    theme_minimal()
  ggsave(file.path(out_dir, "p_by_method_unmarked_effort.png"), p_fig, width = 9, height = 5, dpi = 200)
  
  readr::write_csv(comm_sum, file.path(out_dir, "community_detection_unmarked_effort_psi_cov.csv"))
}

tab_beta <- tab_beta %>% distinct(species, beta_alt, se_alt, beta_doy, se_doy)

g1 <- ggplot(tab_beta, aes(x = beta_alt)) + geom_density(na.rm=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Effect on occupancy (β_ALTITUDE_z)", y="Density", title="Community distribution of altitude effects on ψ") +
  theme_minimal()

g2 <- ggplot(tab_beta, aes(x = beta_doy)) + geom_density(na.rm=TRUE) +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Effect on occupancy (β_DOY_z)", y="Density", title="Community distribution of date (DOY) effects on ψ") +
  theme_minimal()

ggsave(file.path(out_dir, "psi_effects_altitude_density.png"), g1, width=6, height=4, dpi=200)
ggsave(file.path(out_dir, "psi_effects_doy_density.png"), g2, width=6, height=4, dpi=200)

# Violin simplifié ordonné
if (nrow(tab_p) > 0) {
  comm_sum <- tab_p %>% group_by(method) %>%
    summarise(p_med = median(p_hat, na.rm=TRUE), .groups="drop")
  tab_p$method <- factor(tab_p$method, levels = c("GPD", "Pitfall10", "Pitfall8", "Pitfall6", "Pitfall4", "Pitfall2", "DVAC", "TM"))
  
  p_fig2 <- ggplot(tab_p, aes(x = method, y = p_hat)) +
    geom_violin(fill = "grey92") +
    geom_jitter(color = "grey", size = 2, width = 0.3) +
    geom_point(data = comm_sum, aes(y = p_med), color = "red", size = 3) +
    labs(y = "p̂ (détectabilité, occupancy unmarked;\n effort médian par méthode)",
         x = NULL, title = "Détectabilité par méthode —\najustée pour l’effort (unmarked::occu)") +
    theme_minimal()
  ggsave(file.path(out_dir, "p_by_method_unmarked_effort_psi_cov.png"), p_fig2, width = 8, height = 5, dpi = 200)
}
