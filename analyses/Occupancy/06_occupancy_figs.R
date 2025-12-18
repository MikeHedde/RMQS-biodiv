# ============================================================
# FIGURES & SYNTHÈSES p̂ et β (ψ~ DO, HABY)
# Requiert: tab_p, tab_beta (créés au script 05)
# ============================================================

tab_p     <- readr::read_csv(file.path(out_dir, "p_hat_by_method_unmarked_full.csv"), show_col_types = FALSE)
tab_beta0 <- readr::read_csv(file.path(out_dir, "psi_coef_cov_by_species.csv"), show_col_types = FALSE)

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

tab_beta_long <- tab_beta0 %>%
  distinct(species,
           beta_doy, se_doy,
           beta_altitude, se_altitude,
           beta_ndvi_m, se_ndvi_m,
           beta_ndvi_sd, se_ndvi_sd) %>%
  pivot_longer(
    cols = c(beta_doy, beta_altitude, beta_ndvi_m, beta_ndvi_sd,
             se_doy,   se_altitude,   se_ndvi_m,   se_ndvi_sd),
    names_to = c(".value", "covar"),
    names_pattern = "(beta|se)_(.*)"
  ) %>%
  mutate(
    covar = recode(covar,
                   "altitude" = "Altitude",
                   "doy"      = "Day of year",
                   "ndvi_m"   = "mean NDVI",
                   "ndvi_sd"  = "sd NDVI"),
    covar = factor(covar, levels = c("Altitude", "Day of year", "mean NDVI", "sd NDVI"))
  )

p <- ggplot(tab_beta_long, aes(x = beta)) +
  geom_density(na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ covar, ncol = 2) +
  labs(x = "Effect on occupancy (β)", y = "Density") +
  theme_minimal()

ggsave(file.path(out_dir, "psi_effects_covar_density.png"), p, width=6, height=4, dpi=200)

# Violin simplifié ordonné
if (nrow(tab_p) > 0) {
  comm_sum <- tab_p %>% group_by(method) %>%
    summarise(p_med = median(p_hat, na.rm=TRUE), .groups="drop")
  tab_p$method <- factor(tab_p$method, levels = c("GPD", "Pitfall10", "Pitfall8", "Pitfall6", "Pitfall4", "Pitfall2", "DVAC", "TM"))
  
  p_fig2 <- ggplot(tab_p, aes(x = method, y = p_hat)) +
    geom_boxplot(fill = "grey92", outliers = F) +
    geom_jitter(color = "grey", size = 2, width = 0.3) +
    labs(y = "p̂ (detectability)",
         x = NULL, #title = "Détectabilité par méthode —\najustée pour l’effort (unmarked::occu)"
         ) +
    geom_point(data = comm_sum, aes(y = p_med), color = "red", size = 3) +
    theme_minimal()
  ggsave(file.path(out_dir, "p_by_method_unmarked_effort_psi_cov.png"), p_fig2, width = 8, height = 5, dpi = 200)

tab_p2 <- tab_p %>% 
  left_join(comm_sum)%>%
  mutate(rel_p = p_hat-p_med,
         genus = word(species,1))
  

  p_fig2bis <- ggplot(tab_p2, aes(x = method, y = rel_p)) +
    geom_jitter(aes(colour = species), width = 0.2) +
    geom_hline(yintercept = 0)+
    #geom_point(data = comm_sum, aes(y = p_med), color = "red", size = 1) +
    labs(y = "p̂ (relative detectability)\n (p_hat - p_median)",
         x = NULL, #title = "Différence de détectabilité (p_hat - p_median) par méthode —\najustée pour l’effort (unmarked::occu)"
         ) +
    facet_wrap(species~.)+
    coord_flip()+
    theme_minimal()+
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "p_by_method_unmarked_effort_psi_cov_sp.png"), p_fig2bis, width = 12, height = 12, dpi = 600)
  }
