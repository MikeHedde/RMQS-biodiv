# ============================================================
# 11_diversity_stats_clean.R
# - Coverage: glmmTMB beta logit (method x habitat) + trend ntrap
# - Diversity: lmer per metric with coverage + method x habitat
# - Residual plot: coverage-adjusted deviations by method, habitat, metric
# ============================================================

# ---------------------------
# 0) Helpers
# ---------------------------
clip01 <- function(x, eps = 1e-4) pmin(pmax(x, eps), 1 - eps)

as_ntrap <- function(method_chr) {
  dplyr::case_when(
    method_chr == "Pitfall2"  ~ 2,
    method_chr == "Pitfall4"  ~ 4,
    method_chr == "Pitfall6"  ~ 6,
    method_chr == "Pitfall8"  ~ 8,
    method_chr == "Pitfall10" ~ 10,
    method_chr == "GPD"       ~ 4,
    TRUE ~ NA_real_
  )
}

# ---------------------------
# 1) Habitat join (une seule fois, propre)
# ---------------------------
hab <- read.csv("data/derived-data/liste_habitat.csv", sep = ";", header = TRUE) %>%
  rename(id_station = STATION) %>%
  mutate(id_station = as.numeric(id_station)) %>%
  # garde seulement les colonnes utiles si besoin :
  # select(id_station, HABITAT_CODE, HABITAT_TYPE, ...)
  mutate(
    HABITAT_CODE = factor(
      HABITAT_CODE,
      levels = c("Forest", "Schrubland", "Grassland / Meadow", "Crop")
    )
  )

# ---------------------------
# 2) Préparer div_cov2 (coverage + diversité + habitat + ids)
# ---------------------------

div_cov <- read.csv("data/derived-data/Hill_div_results.csv")

div_cov2 <- div_cov %>%
  mutate(
    cov_clip = clip01(coverage),
    ntrap    = as_ntrap(as.character(method))
  ) 

# Assure l’ordre des méthodes (utile pour emmeans/plots)
div_cov2 <- div_cov2 %>%
  mutate(method = factor(method, levels = c("Pitfall2","Pitfall4","Pitfall6","Pitfall8","Pitfall10","GPD"))) %>%
  filter(!is.na(HABITAT_CODE))

cov_sum <- div_cov2 %>%
  group_by(method, HABITAT_CODE, ntrap) %>%
  summarise(
    cov_med = median(coverage, na.rm = TRUE),
    cov_q1  = quantile(coverage, 0.25, na.rm = TRUE),
    cov_q3  = quantile(coverage, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

cov_sum$HABITAT_CODE <- factor(cov_sum$HABITAT_CODE, levels = c("Crop", "Grassland / Meadow", "Schrubland", "Forest"))

cov_pit <- cov_sum %>%
  filter(method != "GPD")

cov_gpd <- cov_sum %>%
  filter(method == "GPD") %>%
  mutate(ntrap_gpd = 4.3)   # décalage visuel

cov_plot <- ggplot(cov_pit,
                   aes(x = ntrap, y = cov_med)) +
  ## Pitfall
  geom_point(size = 3, colour = "black") +
  geom_errorbar(aes(ymin = cov_q1, ymax = cov_q3),
                width = 0.15, linewidth = 0.6, colour = "black") +
  
  ## GPD (couche dédiée)
  geom_point(
    data = cov_gpd,
    aes(x = ntrap_gpd, y = cov_med),
    shape = 21, size = 3,
    fill = "white", colour = "black", stroke = 0.8
  ) +
  geom_errorbar(
    data = cov_gpd,
    aes(x = ntrap_gpd, ymin = cov_q1, ymax = cov_q3),
    width = 0.15, linewidth = 0.6, colour = "black"
  ) +
  
  facet_grid(~ HABITAT_CODE) +
  
  scale_x_continuous(
    breaks = c(2, 4, 6, 8, 10),
    labels = c("2", "4", "6", "8", "10"),
    name = "Number of pitfall traps\n(sampling effort)"
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    name = "Inventory coverage\n(Good’s coverage)"
  ) +
  
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 8))
  )

cov_plot
ggsave(file.path(out_dir, "diversity/coverage.png"),
       cov_plot, width = 12, height = 7, dpi = 300)

# Long
div_long <- div_cov2 %>%
  select(id_station, HABITAT_CODE, method, coverage,
         richness, Hill_q1, Hill_q2) %>%
  pivot_longer(
    cols = c(richness, Hill_q1, Hill_q2),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric,
                    levels = c("richness","Hill_q1","Hill_q2"),
                    labels = c("Richness (q = 0)","Hill diversity (q = 1)","Hill diversity (q = 2)")
    )
  )

# -------------------------
# 1) Choisir une couverture de référence (éviter extrapolation)
# -------------------------
cov_ref <- 0.95
# Optionnel : vérif rapide
# div_long %>% summarise(min=min(coverage,na.rm=TRUE), max=max(coverage,na.rm=TRUE))

# -------------------------
# 2) Fit + tests LRT + contrasts à coverage fixée
# -------------------------
fit_one_metric <- function(df) {
  
  df0 <- df %>%
    select(value, coverage, method, HABITAT_CODE, id_station) %>%
    drop_na() %>%
    droplevels()
  
  # id_station doit être un facteur pour lmer
  df0$id_station <- factor(df0$id_station)
  
  m0 <- lmer(value ~ coverage + (1 | id_station),
             data = df0, REML = FALSE)
  
  m1 <- lmer(value ~ coverage + method + HABITAT_CODE + (1 | id_station),
             data = df0, REML = FALSE)
  
  m2 <- lmer(value ~ coverage + method * HABITAT_CODE + (1 | id_station),
             data = df0, REML = FALSE)
  
  list(
    data = df0,
    m0 = m0,
    m1 = m1,
    m2 = m2,
    lrt_method = anova(m0, m1),
    lrt_inter  = anova(m1, m2),
    emm = emmeans(m2, ~ method | HABITAT_CODE,
                  at = list(coverage = cov_ref)),
    contr = contrast(
      emmeans(m2, ~ method | HABITAT_CODE,
              at = list(coverage = cov_ref)),
      method = "trt.vs.ctrl",
      ref = "Pitfall10",
      adjust = "holm"
    )
  )
}

fits <- div_long %>%
  group_by(metric) %>%
  group_split() %>%
  setNames(levels(div_long$metric)) %>%
  lapply(fit_one_metric)

# -------------------------
# 3) Export tables “propres” pour annexe
# -------------------------
lrt_tbl <- bind_rows(lapply(names(fits), function(nm){
  tibble::tibble(
    metric = nm,
    test_method = fits[[nm]]$lrt_method$`Pr(>Chisq)`[2],
    test_inter  = fits[[nm]]$lrt_inter$`Pr(>Chisq)`[2]
  )
}))

# Contrastes (Δ vs Pitfall10) à coverage=cov_ref
contr_tbl <- bind_rows(lapply(names(fits), function(nm){
  as.data.frame(summary(fits[[nm]]$contr)) %>%
    mutate(metric = nm, coverage_ref = cov_ref)
}))

write.csv(lrt_tbl,  file.path(out_dir, "diversity/LRT_method_interaction_by_metric.csv"), row.names = FALSE)
write.csv(contr_tbl,file.path(out_dir, "diversity/contrasts_vs_Pitfall10_at_covref.csv"), row.names = FALSE)

# -------------------------
# 4) Figure UNIQUE : effets méthode restants à couverture fixée
# -------------------------
plot_df <- contr_tbl %>%
  # garder colonnes standard emmeans
  mutate(
    lwr = estimate-SE,
    upr = estimate+SE
  ) %>%
  mutate(
    contrast = factor(contrast, levels = rev(unique(contrast))),
    metric   = factor(metric, levels = levels(div_long$metric))
  )

p <- ggplot(plot_df, aes(x = estimate, y = contrast)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = lwr, xmax = upr), height = 0.2) +
  geom_point(size = 2) +
  facet_grid(HABITAT_CODE ~ metric, scales = "free_x") +
  labs(
    x = paste0("Method effect on diversity at fixed coverage (Δ vs Pitfall10; coverage = ", cov_ref, ")"),
    y = NULL
  ) +
  theme_minimal()+
  theme(text = element_text(size = 18))

ggsave(file.path(out_dir, "diversity/Fig_method_effects_at_fixed_coverage.png"),
       p, width = 12, height = 10, dpi = 300)
