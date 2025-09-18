# ============================================================
# TRAITS → DÉTECTABILITÉ (p̂) ~ méthode * log(masse)
# Sorties: modèles, figures
# ============================================================

tab_p <- readr::read_csv(file.path(out_dir, "p_hat_by_method_unmarked_full.csv"), show_col_types = FALSE)

trait <- raw %>%
  filter(PROJET == target_project) %>%
  select(LB_NOM, `MASSE (mg)`) %>%
  rename(species = LB_NOM) %>%
  group_by(species) %>%
  summarise(mass = median(`MASSE (mg)`, na.rm = TRUE), .groups="drop")

dat <- tab_p %>% left_join(trait, by = "species")

# Mise à l’échelle logit + poids méta depuis IC95%
dat2 <- dat %>%
  mutate(
    p_clip  = pmin(pmax(p_hat, 1e-5), 1 - 1e-5),
    logit_p = logit(p_clip),
    se_p    = ifelse(!is.na(lcl) & !is.na(ucl), (ucl - lcl) / (2 * 1.96), NA_real_),
    se_logit = ifelse(!is.na(se_p), se_p / (p_clip * (1 - p_clip)), NA_real_),
    weight   = ifelse(is.finite(1/(se_logit^2)), 1/(se_logit^2), 1),
    log_mass = log10(mass)
  ) %>%
  filter(!is.na(log_mass), is.finite(log_mass)) %>%
  mutate(species = factor(species))

m_mass_int <- lme4::lmer(logit_p ~ method * log_mass + (1 | species),
                         data = dat2, weights = weight, REML = TRUE)
s <- summary(m_mass_int); capture.output(s, file = file.path(out_dir, "model_traits_detectability_summary.txt"))

trends <- emmeans::emtrends(m_mass_int, ~ method, var = "log_mass")
capture.output(trends, file = file.path(out_dir, "model_traits_emtrends.txt"))
capture.output(pairs(trends), file = file.path(out_dir, "model_traits_emtrends_pairs.txt"))

# Prédictions lissées sur gamme de masses par méthode
newgrid <- dat2 %>%
  group_by(method) %>%
  summarise(m_min = quantile(log_mass, 0.05), m_max = quantile(log_mass, 0.95), .groups="drop") %>%
  rowwise() %>% do(tibble(method = .$method, log_mass = seq(.$m_min, .$m_max, length.out = 100))) %>% ungroup()

pred <- cbind(newgrid,
              predict(m_mass_int, newdata = newgrid, re.form = NA, se.fit = TRUE)) %>%
  as_tibble() %>%
  mutate(p_fit = plogis(fit), p_lcl = plogis(fit - 1.96*se.fit), p_ucl = plogis(fit + 1.96*se.fit))

ggplot(dat2, aes(x = log_mass, y = p_clip, colour = method)) +
  geom_point(alpha = 0.5) +
  geom_line(data = pred, aes(y = p_fit), linewidth = 1) +
  #geom_ribbon(data = pred, aes(ymin = p_lcl, ymax = p_ucl, fill = method), alpha = 0.15, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(method~.) +
  labs(x = "log10(body mass, mg)", y = "Detectability p̂ (back-transformed)",
       title = "Effect of body mass on detectability by method") +
  theme_minimal() +
  theme(legend.position = "none") ->
  p_mass

ggsave(file.path(out_dir, "detectability_vs_mass_by_method.png"), p_mass, width = 9, height = 6, dpi = 200)
