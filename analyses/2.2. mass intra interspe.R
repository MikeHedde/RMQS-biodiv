## Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, tidyverse, stringr,
                 car, broom, paletteer, cowplot, ggpubr,
                 moments, mgcv)

#Matrice de données environnementales
hab <- read.csv("data/derived-data/liste_habitat.csv", sep = ";", h = T) %>%
  select(STATION, HABITAT_TYPE) %>%
  mutate(STATION = as.factor(STATION)) %>%
  unique()

# Définition des mots à exclure dans les remarques sur les individus
mots_a_exclure <- c("endommagé", "a peser", "non", "abdomen creux", "corps coupé en 2", "abdomen absent", "repeser", "absence d'élytres", "vide")

# Chargement des données de masses des Carabidae et fusion avec les données environnementales
mass0 <- read.csv("data/raw-data/1.faune/carabidae.csv", h=TRUE, sep=";", dec=".", fileEncoding = "ISO-8859-1") %>%
  separate(ID_ECHANTILLON, "_", into = c("pj", "YEAR", "CITY", "STATION", "ECH")) %>%
  dplyr::filter(!grepl(paste(mots_a_exclure, collapse = "|"), 
                REMARQUE_INDIVIDU, ignore.case = TRUE),
                RANK %in% c("ES", "S"),
                ABONDANCE_TOTALE == 1,
                !is.na(MASSE)) %>%
  mutate(MASSE = as.numeric(MASSE), STATION = as.factor(STATION)) %>%
  select(STATION, NOM_VALIDE, MASSE, ABONDANCE_TOTALE) %>%
  left_join(hab)

# distribution des espèces
ab_prop <- mass0 %>%
  group_by(STATION, HABITAT_TYPE, NOM_VALIDE) %>%
  summarise(nb = n(), m_mass = mean(MASSE), 
            etendue = max(MASSE)-min(MASSE),
            etendue_norm = etendue/mean(MASSE),
            skewn = skewness(MASSE), kurt = kurtosis(MASSE))  %>%
  mutate(freq = nb/sum(nb)*100)

#En présence de P. cupreus
dom_prop <- ab_prop %>%
  #filter(hab_class %in% c("crop", "non-crop")) %>%
  group_by(STATION) %>%
  filter(any(NOM_VALIDE == "Poecilus cupreus (Linnaeus, 1758)")) %>%
  mutate(
    mass_poecilus = m_mass[NOM_VALIDE == "Poecilus cupreus (Linnaeus, 1758)"],
    ecart_mass = m_mass - mass_poecilus,
    ecart_mass_prc = ecart_mass/mass_poecilus,
    ecart_mass_ratio = m_mass / mass_poecilus,
    freq_poecilus = freq[NOM_VALIDE == "Poecilus cupreus (Linnaeus, 1758)"]
  ) %>%
  ungroup() %>%
  mutate(
    pc_mass_class = case_when(
      mass_poecilus < 80 ~ "low",
      mass_poecilus >= 80 ~ "high"),
    pc_freq_class = ntile(freq_poecilus, 4)) %>%
  group_by(STATION, HABITAT_TYPE, NOM_VALIDE, pc_mass_class)%>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  #filter(!NOM_VALIDE == "Poecilus cupreus (Linnaeus, 1758)") %>%
  mutate(ecart_mass_class = case_when(
      ecart_mass >= 0 ~ "2.heavier",
      ecart_mass < 0 ~ "1.lighter"
    ))
  
N = 7
dom_prop1 <- dom_prop %>%
  filter(nb > N,
         NOM_VALIDE != "Poecilus cupreus (Linnaeus, 1758)")

dom_prop2 <- dom_prop %>%
  filter(nb > N,
         NOM_VALIDE != "Poecilus cupreus (Linnaeus, 1758)") %>%
  mutate(pc_freq_class = case_when(pc_freq_class == 1 ~ 
                                     paste("1st quartile\n", round(min(freq_poecilus[pc_freq_class == 1]), 1), "-", round(max(freq_poecilus[pc_freq_class == 1]), 1), "%"),
                                   pc_freq_class == 2 ~ 
                                     paste("2nd quartile\n", round(min(freq_poecilus[pc_freq_class == 2]), 1), "-", round(max(freq_poecilus[pc_freq_class == 2]), 1), "%"),
                                   pc_freq_class == 3 ~ 
                                     paste("3rd quartile\n", round(min(freq_poecilus[pc_freq_class == 3]), 1), "-", round(max(freq_poecilus[pc_freq_class == 3]), 1), "%"),
                                   pc_freq_class == 4 ~ 
                                     paste("4th quartile\n", round(min(freq_poecilus[pc_freq_class == 4]), 1), "-", round(max(freq_poecilus[pc_freq_class == 4]), 1), "%")
  ))

## Calcul des moments de la distribution
breath <- ggplot(dom_prop2, aes(x=abs(log(ecart_mass_ratio)), 
                                y = etendue_norm))+
  geom_smooth(aes(colour = as.factor(pc_freq_class)),
              method = "lm", se = T, formula = y ~ poly(x, 1))+
  labs(colour = "P. cupreus \nfrequency quantile",
       x = "",
       y = "Niche breath \n(normalised width of mass distribution)")+
  geom_hline(yintercept = 0)+ geom_vline(xintercept = 0)+
  geom_point(aes(colour = as.factor(pc_freq_class)), alpha = 0.4)+
  scale_colour_paletteer_d("calecopal::sbchannel") +
  facet_grid(pc_freq_class~.)+
  theme_bw()+
  theme(legend.position = "none")

kurt <- ggplot(dom_prop2, aes(x=abs(log(ecart_mass_ratio)), y = kurt))+
  geom_smooth(aes(colour = as.factor(pc_freq_class)),
              method = "lm", se = T, formula = y ~ poly(x, 1))+
  labs(colour = "P. cupreus \nfrequency quantile",
       x = "Mass ratio \n(log of absolute species i / P. cupreus)",
       y = "Kurtosis of \nspecies mass distribution")+
  geom_hline(yintercept = 0)+ geom_vline(xintercept = 0)+
  geom_point(aes(colour = as.factor(pc_freq_class)), alpha = 0.4)+
  scale_colour_paletteer_d("calecopal::sbchannel") +
  facet_grid(pc_freq_class~.)+
  theme_bw()+
  theme(legend.position = "none")


skewn <- ggplot(dom_prop2, aes(x=abs(log(ecart_mass_ratio)), y = skewn))+
  geom_smooth(aes(colour = as.factor(pc_freq_class)),
              method = "lm", se = T, formula = y ~ poly(x, 1)
  )+
  labs(colour = "P. cupreus \nfrequency quantile",
       x = "",
       y = "Skewness of species \nmass distribution")+
  geom_hline(yintercept = 0)+ geom_vline(xintercept = 0)+
  geom_point(aes(colour = as.factor(pc_freq_class)), alpha = 0.4)+
  scale_colour_paletteer_d("calecopal::sbchannel") +
  facet_grid(pc_freq_class~.)+
  theme_bw()+
  theme(legend.position = "none")

p1 <- ggarrange(breath, kurt, skewn, common.legend = F, ncol = 3, 
          align = "h")

###############
summary(aov(skewn~pc_freq_class*abs(log(ecart_mass_ratio)), data=dom_prop1))
d = dom_prop1

##########
mod_skew1 <- lm(skewn~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==1,])
mod_skew2 <- lm(skewn~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==2,])
mod_skew3 <- lm(skewn~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==3,])
mod_skew4 <- lm(skewn~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==4,])

lwr_skew_slope_1 <- confint(mod_skew1)["abs(log(ecart_mass_ratio))", 1]
upr_skew_slope_1 <- confint(mod_skew1)["abs(log(ecart_mass_ratio))", 2]
lwr_skew_slope_2 <- confint(mod_skew2)["abs(log(ecart_mass_ratio))", 1]
upr_skew_slope_2 <- confint(mod_skew2)["abs(log(ecart_mass_ratio))", 2]
lwr_skew_slope_3 <- confint(mod_skew3)["abs(log(ecart_mass_ratio))", 1]
upr_skew_slope_3 <- confint(mod_skew3)["abs(log(ecart_mass_ratio))", 2]
lwr_skew_slope_4 <- confint(mod_skew4)["abs(log(ecart_mass_ratio))", 1]
upr_skew_slope_4 <- confint(mod_skew4)["abs(log(ecart_mass_ratio))", 2]

lwr_skew_interc_1 <- confint(mod_skew1)["(Intercept)", 1]
upr_skew_interc_1 <- confint(mod_skew1)["(Intercept)", 2]
lwr_skew_interc_2 <- confint(mod_skew2)["(Intercept)", 1]
upr_skew_interc_2 <- confint(mod_skew2)["(Intercept)", 2]
lwr_skew_interc_3 <- confint(mod_skew3)["(Intercept)", 1]
upr_skew_interc_3 <- confint(mod_skew3)["(Intercept)", 2]
lwr_skew_interc_4 <- confint(mod_skew4)["(Intercept)", 1]
upr_skew_interc_4 <- confint(mod_skew4)["(Intercept)", 2]
###########################################

mod_kurt1 <- lm(kurt~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==1,])
mod_kurt2 <- lm(kurt~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==2,])
mod_kurt3 <- lm(kurt~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==3,])
mod_kurt4 <- lm(kurt~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==4,])

lwr_kurt_slope_1 <- confint(mod_kurt1)["abs(log(ecart_mass_ratio))", 1]
upr_kurt_slope_1 <- confint(mod_kurt1)["abs(log(ecart_mass_ratio))", 2]
lwr_kurt_slope_2 <- confint(mod_kurt2)["abs(log(ecart_mass_ratio))", 1]
upr_kurt_slope_2 <- confint(mod_kurt2)["abs(log(ecart_mass_ratio))", 2]
lwr_kurt_slope_3 <- confint(mod_kurt3)["abs(log(ecart_mass_ratio))", 1]
upr_kurt_slope_3 <- confint(mod_kurt3)["abs(log(ecart_mass_ratio))", 2]
lwr_kurt_slope_4 <- confint(mod_kurt4)["abs(log(ecart_mass_ratio))", 1]
upr_kurt_slope_4 <- confint(mod_kurt4)["abs(log(ecart_mass_ratio))", 2]

lwr_kurt_interc_1 <- confint(mod_kurt1)["(Intercept)", 1]
upr_kurt_interc_1 <- confint(mod_kurt1)["(Intercept)", 2]
lwr_kurt_interc_2 <- confint(mod_kurt2)["(Intercept)", 1]
upr_kurt_interc_2 <- confint(mod_kurt2)["(Intercept)", 2]
lwr_kurt_interc_3 <- confint(mod_kurt3)["(Intercept)", 1]
upr_kurt_interc_3 <- confint(mod_kurt3)["(Intercept)", 2]
lwr_kurt_interc_4 <- confint(mod_kurt4)["(Intercept)", 1]
upr_kurt_interc_4 <- confint(mod_kurt4)["(Intercept)", 2]
#######################

mod_breath1 <- lm(etendue_norm~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==1,])
mod_breath2 <- lm(etendue_norm~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==2,])
mod_breath3 <- lm(etendue_norm~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==3,])
mod_breath4 <- lm(etendue_norm~abs(log(ecart_mass_ratio)), data=d[d$pc_freq_class ==4,])

lwr_breath_slope_1 <- confint(mod_breath1)["abs(log(ecart_mass_ratio))", 1]
upr_breath_slope_1 <- confint(mod_breath1)["abs(log(ecart_mass_ratio))", 2]
lwr_breath_slope_2 <- confint(mod_breath2)["abs(log(ecart_mass_ratio))", 1]
upr_breath_slope_2 <- confint(mod_breath2)["abs(log(ecart_mass_ratio))", 2]
lwr_breath_slope_3 <- confint(mod_breath3)["abs(log(ecart_mass_ratio))", 1]
upr_breath_slope_3 <- confint(mod_breath3)["abs(log(ecart_mass_ratio))", 2]
lwr_breath_slope_4 <- confint(mod_breath4)["abs(log(ecart_mass_ratio))", 1]
upr_breath_slope_4 <- confint(mod_breath4)["abs(log(ecart_mass_ratio))", 2]

lwr_breath_interc_1 <- confint(mod_breath1)["(Intercept)", 1]
upr_breath_interc_1 <- confint(mod_breath1)["(Intercept)", 2]
lwr_breath_interc_2 <- confint(mod_breath2)["(Intercept)", 1]
upr_breath_interc_2 <- confint(mod_breath2)["(Intercept)", 2]
lwr_breath_interc_3 <- confint(mod_breath3)["(Intercept)", 1]
upr_breath_interc_3 <- confint(mod_breath3)["(Intercept)", 2]
lwr_breath_interc_4 <- confint(mod_breath4)["(Intercept)", 1]
upr_breath_interc_4 <- confint(mod_breath4)["(Intercept)", 2]
#############################


conf_df <- data.frame(
  Quantile = rep(c("Q1", "Q2", "Q3", "Q4"), each = 6),
  Moment = rep(c("Breath", "Kurt.", "Skewn."), each = 2),
  Type = rep(c("Slope", "Intercept"), times = 2),
  Lwr = c(lwr_breath_slope_1, lwr_breath_interc_1,lwr_kurt_slope_1, lwr_kurt_interc_1, lwr_skew_slope_1, lwr_skew_interc_1,
          lwr_breath_slope_2, lwr_breath_interc_2,lwr_kurt_slope_2, lwr_kurt_interc_2, lwr_skew_slope_2, lwr_skew_interc_2,
          lwr_breath_slope_3, lwr_breath_interc_3,lwr_kurt_slope_3, lwr_kurt_interc_3, lwr_skew_slope_3, lwr_skew_interc_3,
          lwr_breath_slope_4, lwr_breath_interc_4,lwr_kurt_slope_4, lwr_kurt_interc_4, lwr_skew_slope_4, lwr_skew_interc_4),
  Upr = c(upr_breath_slope_1, upr_breath_interc_1,upr_kurt_slope_1, upr_kurt_interc_1, upr_skew_slope_1, upr_skew_interc_1,
          upr_breath_slope_2, upr_breath_interc_2,upr_kurt_slope_2, upr_kurt_interc_2, upr_skew_slope_2, upr_skew_interc_2,
          upr_breath_slope_3, upr_breath_interc_3,upr_kurt_slope_3, upr_kurt_interc_3, upr_skew_slope_3, upr_skew_interc_3,
          upr_breath_slope_4, upr_breath_interc_4,upr_kurt_slope_4, upr_kurt_interc_4, upr_skew_slope_4, upr_skew_interc_4),
  Estimate = c(coef(mod_breath1)[2], coef(mod_breath1)[1], coef(mod_kurt1)[2], coef(mod_kurt1)[1], coef(mod_skew1)[2], coef(mod_skew1)[1],
               coef(mod_breath2)[2], coef(mod_breath2)[1], coef(mod_kurt2)[2], coef(mod_kurt2)[1], coef(mod_skew2)[2], coef(mod_skew2)[1],
               coef(mod_breath3)[2], coef(mod_breath3)[1], coef(mod_kurt3)[2], coef(mod_kurt3)[1], coef(mod_skew3)[2], coef(mod_skew3)[1],
               coef(mod_breath4)[2], coef(mod_breath4)[1], coef(mod_kurt4)[2], coef(mod_kurt4)[1], coef(mod_skew4)[2], coef(mod_skew4)[1])
  )

# Graphique des intervalles de confiance
graph_conf <- 
  ggplot(conf_df, aes(x = Moment, y = Estimate, color = Quantile)) +
  geom_pointrange(aes(ymin = Lwr, ymax = Upr), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_colour_paletteer_d("calecopal::sbchannel") +
  labs(y = "Estimate", x = "Moment") +
  facet_grid(Type~., scales = "free")+
  theme_bw()
  theme(legend.position = "none")

  ggarrange(p1, graph_conf, common.legend = F, ncol = 2, 
            align = "h", widths = c(2.2,1))
  