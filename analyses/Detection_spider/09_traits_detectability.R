# ============================================================
# TRAITS -> DETECTABILITY
# Clean modelling structure:
# 1. Body length: additive vs method interaction
# 2. Wing development: additive vs method interaction, controlling for body size
# 3. Associated figures
# ============================================================

# ------------------------------------------------------------
# 1) Build common trait-detectability dataset
# ------------------------------------------------------------

trait_gn <- trait_sp %>%
  mutate(genus = word(species, 1)) %>%
  group_by(genus) %>%
  summarise(
    pct_winged_gn = median(pct_winged, na.rm = TRUE),
    .groups = "drop"
  )

dat_traits <- tab_p %>%
  left_join(BL, by = "species") %>%
  left_join(trait_sp, by = "species") %>%
  mutate(genus = word(species, 1)) %>%
  left_join(trait_gn, by = "genus") %>%
  mutate(
    pct_winged = ifelse(is.na(pct_winged), pct_winged_gn, pct_winged),
    wing_imputed_genus = is.na(pct_winged) & !is.na(pct_winged_gn),
    
    full_wing = case_when(
      pct_winged == 1 ~ "fully_winged",
      pct_winged <  1 ~ "not_fully_winged",
      TRUE ~ NA_character_
    ),
    
    method = factor(
      method,
      levels = c("Pitfall2", "Pitfall4", "Pitfall6",
                 "Pitfall8", "Pitfall10", "GPD")
    ),
    
    log_BL = log10(mean_BL),
    
    p_clip = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    lcl_clip = pmin(pmax(lcl, 1e-5), 1 - 1e-5),
    ucl_clip = pmin(pmax(ucl, 1e-5), 1 - 1e-5),
    
    logit_p = qlogis(p_clip),
    se_logit = (qlogis(ucl_clip) - qlogis(lcl_clip)) / (2 * 1.96),
    weight = 1 / se_logit^2,
    
    species = factor(species)
  ) %>%
  filter(
    !is.na(method),
    !is.na(log_BL),
    !is.na(pct_winged),
    is.finite(logit_p),
    is.finite(se_logit),
    se_logit > 0,
    is.finite(weight),
    weight > 0
  )

write_csv(dat_traits, file.path(out_dir, "trait/traits_detectability_dataset.csv"))

# ------------------------------------------------------------
# 2) Body length models
# ------------------------------------------------------------

m_bl_add <- lmer(
  logit_p ~ method + log_BL + (1 | species),
  data = dat_traits,
  weights = weight,
  REML = TRUE
)

m_bl_int <- lmer(
  logit_p ~ method * log_BL + (1 | species),
  data = dat_traits,
  weights = weight,
  REML = TRUE
)

capture.output(
  summary(m_bl_add),
  car::Anova(m_bl_add, type = 3),
  file = file.path(out_dir, "trait/model_BL_additive_summary.txt")
)

capture.output(
  summary(m_bl_int),
  car::Anova(m_bl_int, type = 3),
  file = file.path(out_dir, "trait/model_BL_method_interaction_summary.txt")
)

bl_trends <- emtrends(m_bl_int, ~ method, var = "log_BL")
bl_trends_pairs <- pairs(bl_trends, adjust = "tukey")

capture.output(
  bl_trends,
  bl_trends_pairs,
  file = file.path(out_dir, "trait/model_BL_emtrends_by_method.txt")
)

# ------------------------------------------------------------
# 3) Wing development models
#    Additive effect vs method-dependent effect,
#    while controlling for body size
# ------------------------------------------------------------

m_wing_add <- lmer(
  logit_p ~ method + log_BL + full_wing + (1 | species),
  data = dat_traits,
  weights = weight,
  REML = TRUE
)

m_wing_int <- lmer(
  logit_p ~ method + log_BL + full_wing + method:full_wing + (1 | species),
  data = dat_traits,
  weights = weight,
  REML = TRUE
)

capture.output(
  summary(m_wing_add),
  car::Anova(m_wing_add, type = 3),
  file = file.path(out_dir, "trait/model_wing_additive_summary.txt")
)

capture.output(
  summary(m_wing_int),
  car::Anova(m_wing_int, type = 3),
  file = file.path(out_dir, "trait/model_wing_method_interaction_summary.txt")
)

wing_emm <- emmeans(
  m_wing_int,
  ~ full_wing | method
)

wing_contrasts <- pairs(
  wing_emm,
  adjust = "none"
)

wing_contrasts_by_method <- as.data.frame(wing_contrasts)

wing_effect_by_method <- contrast(
  wing_emm,
  method = "revpairwise",
  by = "method",
  adjust = "none"
) %>%
  as.data.frame()


capture.output(
  wing_emm,
  wing_contrasts,
  file = file.path(out_dir, "traits/model_wing_emmeans_by_method.txt")
)

# ------------------------------------------------------------
# 4) Figure A — Body length effect by method
# ------------------------------------------------------------

log_seq <- seq(
  quantile(dat_traits$log_BL, 0.05, na.rm = TRUE),
  quantile(dat_traits$log_BL, 0.95, na.rm = TRUE),
  length.out = 100
)

pred_bl <- emmeans(
  m_bl_int,
  ~ method | log_BL,
  at = list(log_BL = log_seq)
) %>%
  as.data.frame() %>%
  mutate(
    p_fit = plogis(emmean),
    p_lcl = plogis(lower.CL),
    p_ucl = plogis(upper.CL),
    method = factor(method, levels = levels(dat_traits$method))
  )%>%
  arrange(method, log_BL)

p_BL <- ggplot(dat_traits, aes(x = log_BL, y = p_clip, colour = method)) +
  geom_point(
    alpha = 0.35,
    size = 2,
    colour = "grey40"
  ) +
  geom_ribbon(
    data = pred_bl,
    aes(x = log_BL, ymin = p_lcl, ymax = p_ucl),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  geom_line(
    data = pred_bl %>% arrange(method, log_BL),
    aes(x = log_BL, y = p_fit, group = method),
    inherit.aes = FALSE,
    linewidth = 1.3,
    colour = "black"
  )+
  facet_wrap(~ method, ncol = 6) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Body length log10(mm)",
    y = "Detection probability p̂"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "trait/Fig_trait_BL_detectability_by_method.png"),
  p_BL,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 5) Figure B — Wing development effect by method
#    Predictions at median body length
# ------------------------------------------------------------

pred_wing_cat <- emmeans(
  m_wing_int,
  ~ method * full_wing,
  at = list(
    log_BL = median(dat_traits$log_BL, na.rm = TRUE)
  )
) %>%
  as.data.frame() %>%
  mutate(
    p_fit = plogis(emmean),
    p_lcl = plogis(lower.CL),
    p_ucl = plogis(upper.CL),
    method = factor(method, levels = levels(dat_traits$method)),
    full_wing = factor(
      full_wing,
      levels = c("not_fully_winged", "fully_winged"),
      labels = c("Not fully", "Full")
    )
  )

dat_wing_plot <- dat_traits %>%
  filter(!is.na(full_wing)) %>%
  mutate(
    full_wing = factor(
      full_wing,
      levels = c("not_fully_winged", "fully_winged"),
      labels = c("Not fully", "Full")
    )
  )

p_wing <- ggplot(
  dat_wing_plot,
  aes(x = full_wing, y = p_clip, colour = method)
) +
  geom_violin(
    fill = "grey90",
    colour = NA,
    alpha = 0.5
  )+
  geom_jitter(
    width = 0.12,
    alpha = 0.5,
    size = 1.7
  ) +
  geom_pointrange(
    data = pred_wing_cat,
    aes(
      x = full_wing,
      y = p_fit,
      ymin = p_lcl,
      ymax = p_ucl
    ),
    inherit.aes = FALSE,
    size = 0.8
  ) +
  facet_wrap(~ method, ncol = 6) +
  coord_cartesian(ylim = c(0, 1))+
  labs(
    x = "Wing development",
    y = "Detection probability p̂"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(out_dir, "trait/Fig_trait_wing_binary_detectability_by_method.png"),
  p_wing,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 6) Combined figure
# ------------------------------------------------------------

p_traits <- cowplot::plot_grid(
  p_BL,
  p_wing,
  ncol = 1,
  labels = c("A", "B"),
  align = "v"
)

ggsave(
  file.path(out_dir, "trait/Fig_traits_detectability_combined.png"),
  p_traits,
  width = 12,
  height = 10,
  dpi = 300
)
