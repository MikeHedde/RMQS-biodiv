## Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, tidyverse, stringr, car, broom, paletteer, cowplot, ggpubr)

# Définition des mots à exclure dans les remarques sur les individus
mots_a_exclure <- c("endommagé", "a peser", "non", "abdomen creux", "corps coupé en 2", "abdomen absent", "repeser", "absence d'élytres", "vide")

# Chargement des données environnementales
env0 <- read.csv("data/derived-data/all_env_variables.csv")

station_koppen <- read.csv("data/derived-data/station_koppen.csv", h = T, sep=",") %>%
  select(STATION, kopp) %>%
  mutate(STATION = as.factor(STATION))

remove_outliers <- function(data, column, threshold) {
  mean_val <- mean(data[[column]], na.rm = TRUE)
  sd_val <- sd(data[[column]], na.rm = TRUE)
  data %>% filter(abs(data[[column]] - mean_val) <= threshold * sd_val)
}

# Chargement des données de masses des Carabidae et fusion avec les données environnementales
mass0 <- read.csv("data/raw-data/1.faune/carabidae.csv", h=TRUE, sep=";", dec=".", fileEncoding = "ISO-8859-1") %>%
  separate(ID_ECHANTILLON, "_", into = c("pj", "YEAR", "CITY", "STATION", "ECH")) %>%
  filter(!grepl(paste(mots_a_exclure, collapse = "|"), REMARQUE_INDIVIDU, ignore.case = TRUE)) %>%
  mutate(MASSE = as.numeric(MASSE), STATION = as.factor(STATION)) %>%
  filter(STATION %in% env0$STATION) %>%
  select(STATION, NOM_VALIDE, male, female, MASSE)

env <- subset(env0, STATION %in% mass0$STATION) %>%
  mutate(hab_class = ifelse(HAB == 0, "Open", "Closed")) %>%
  select(STATION, hab_class) %>%
  mutate(STATION = as.factor(STATION)) %>%
  unique()
  
mass <- mass0 %>% 
  left_join(env)%>%
  filter(!is.na(MASSE))

# Transformation des colonnes `male` et `femelle` en une colonne `sex`
mass_long <- mass %>%
  filter(!is.na(MASSE)) %>%
  mutate(sex = case_when(
    male == 1 ~ "male", 
    female == 1 ~ "female"
  )) %>%
  filter(!is.na(sex)) %>%
  left_join(station_koppen)

mean_mass <- mass_long %>%
  group_by(NOM_VALIDE)%>%
  summarise(moy = mean(MASSE), nb = length(MASSE))



# Calcul des masses moyennes par station, habitat, espèce et sexe
rensch_wide <- mass_long %>%
  group_by(STATION, hab_class, NOM_VALIDE, sex) %>%
  summarise(mean_mass = mean(MASSE, na.rm = TRUE), .groups = 'drop') %>%
  filter(all(c("male", "female") %in% sex)) %>%  # Garder uniquement les localités avec les deux sexes présents
  pivot_wider(names_from = sex, values_from = mean_mass) %>%  # Transformation des données en format large
  mutate(log_female = log(female),
         log_male = log(male),
         ecart = (log_female - log_male) / log_female)

# Modèles de régression linéaire pour habitats ouverts et fermés
res.lm_open <- lm(log_female ~ log_male, data = subset(rensch_wide, hab_class == "Open"))
res.lm_closed <- lm(log_female ~ log_male, data = subset(rensch_wide, hab_class == "Closed"))

# Extraction des intervalles de confiance
confint_open <- confint(res.lm_open)["log_male", ]
confint_closed <- confint(res.lm_closed)["log_male", ]
intercept_open <- confint(res.lm_open)["(Intercept)", ]
intercept_closed <- confint(res.lm_closed)["(Intercept)", ]

# Création du dataframe pour la visualisation des intervalles de confiance
conf_df <- data.frame(
  Habitat = rep(c("Open", "Closed"), each = 2),
  Type = rep(c("Slope", "Intercept"), times = 2),
  Lwr = c(confint_open[1], intercept_open[1], confint_closed[1], intercept_closed[1]),
  Upr = c(confint_open[2], intercept_open[2], confint_closed[2], intercept_closed[2]),
  Estimate = c(coef(res.lm_open)["log_male"], coef(res.lm_open)["(Intercept)"],
               coef(res.lm_closed)["log_male"], coef(res.lm_closed)["(Intercept)"])
)

# Graphique des intervalles de confiance
graph_conf <- ggplot(conf_df, aes(x = Type, y = Estimate, color = Habitat)) +
  geom_pointrange(aes(ymin = Lwr, ymax = Upr), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_colour_paletteer_d("suffrager::classic") +
  labs(y = "Estimate", x = "Parameter") +
  theme_bw()+
  theme(legend.position = "none")

# Graphique de la distribution des masses par espèce avec couleur selon l'habitat
distri_esp <- ggplot(subset(mass, !is.na(MASSE)), 
                     aes(x = reorder(NOM_VALIDE, MASSE, median), y = MASSE, color = hab_class)) +
  geom_boxplot(outlier.shape = NA, colour = "black") +
  geom_jitter(alpha=0.1, width = 0.2) +
  labs(y = "Mass\n(mg)", x = "", color = "Habitat") +
  scale_y_continuous(trans='log10') +
  scale_colour_paletteer_d("suffrager::classic") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") 

# Graphique de la relation log(male) vs log(female) selon l'habitat
graph_relation <- ggplot(rensch_wide, 
                         aes(x = log_male, y = log_female, color = hab_class)) +
  geom_abline(slope=1, linetype = "dashed", color="black", size = 1.5) +
  geom_point(alpha = 0.1, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_paletteer_d("suffrager::classic") +
  labs(x = "Male mass \n(mg, log)",
       y = "Female mass \n(mg, log)",
       color = "Habitat") +
  theme_bw() +
  theme(legend.position = "none")

# Affichage des graphiques côte à côte avec ggarrange
graph_combined <- ggarrange(distri_esp, 
                            ggarrange(graph_relation, graph_conf, ncol = 2,
                                      labels = c("B", "C")), 
                            labels = c("A"),
                            ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

graph_combined



########################

# Calcul des masses moyennes par station, habitat, espèce et sexe
rensch_wide_kopp <- mass_long %>%
  group_by(kopp, STATION, hab_class, NOM_VALIDE, sex) %>%
  summarise(mean_mass = mean(MASSE, na.rm = TRUE), .groups = 'drop') %>%
  filter(all(c("male", "female") %in% sex)) %>%  # Garder uniquement les localités avec les deux sexes présents
  pivot_wider(names_from = sex, values_from = mean_mass) %>%  # Transformation des données en format large
  mutate(log_female = log(female),
         log_male = log(male),
         ecart = (log_female - log_male) / log_female)

# Modèles de régression linéaire pour habitats ouverts et fermés
res.lm_open_kopp <- lm(log_female ~ log_male, data = subset(rensch_wide_kopp, hab_class == "Open"))
res.lm_closed_kopp <- lm(log_female ~ log_male, data = subset(rensch_wide_kopp, hab_class == "Closed"))

# Extraction des intervalles de confiance
confint_open_kopp <- confint(res.lm_open_kopp)["log_male", ]
confint_closed_kopp <- confint(res.lm_closed_kopp)["log_male", ]
intercept_open_kopp <- confint(res.lm_open_kopp)["(Intercept)", ]
intercept_closed_kopp <- confint(res.lm_closed_kopp)["(Intercept)", ]

# Création du dataframe pour la visualisation des intervalles de confiance
conf_df_kopp <- data.frame(
  Habitat = rep(c("Open", "Closed"), each = 2),
  Type = rep(c("Slope", "Intercept"), times = 2),
  Lwr = c(confint_open_kopp[1], intercept_open_kopp[1], confint_closed_kopp[1], intercept_closed_kopp[1]),
  Upr = c(confint_open_kopp[2], intercept_open_kopp[2], confint_closed_kopp[2], intercept_closed_kopp[2]),
  Estimate = c(coef(res.lm_open)["log_male"], coef(res.lm_open)["(Intercept)"],
               coef(res.lm_closed)["log_male"], coef(res.lm_closed)["(Intercept)"])
)

# Graphique des intervalles de confiance
graph_conf_kopp <- ggplot(conf_df_kopp, aes(x = Type, y = Estimate, color = Habitat)) +
  geom_pointrange(aes(ymin = Lwr, ymax = Upr), position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_colour_paletteer_d("suffrager::classic") +
  labs(y = "Estimate", x = "Parameter") +
  theme_bw()+
  theme(legend.position = "none")


# Graphique de la relation log(male) vs log(female) selon l'habitat
graph_relation_kopp <- ggplot(rensch_wide_kopp, 
                         aes(x = log_male, y = log_female, color = hab_class)) +
  geom_abline(slope=1, linetype = "dashed", color="black", size = 1.5) +
  geom_point(alpha = 0.1, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_paletteer_d("suffrager::classic") +
  labs(x = "Male mass \n(mg, log)",
       y = "Female mass \n(mg, log)",
       color = "Habitat") +
  #facet_wrap(kopp~.)+
  theme_bw() +
  theme(legend.position = "none")

# Affichage des graphiques côte à côte avec ggarrange
graph_combined_kopp <- ggarrange(distri_esp, 
                            ggarrange(graph_relation_kopp, graph_conf_kopp, ncol = 2,
                                      labels = c("B", "C")), 
                            labels = c("A"),
                            ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

graph_combined_kopp
