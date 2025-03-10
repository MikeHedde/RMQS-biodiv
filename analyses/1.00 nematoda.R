# Analyse des indicateurs nématofauniques
# Auteur : Mickael Hedde
# Date : 10-03-2025
# Objectif : Décrire la valeurs des indicateurs basés sur la nématofaune 

# Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, tidyr, gdm, forcats, ggpubr, paletteer)

# Chargement des données de nematofaune
nem <- read.csv("data/raw-data/1.faune/nematoda.csv", h = T, sep = ";") %>%
  pivot_wider(id_cols = c(STATION, LAND_USE, ALT, X), names_from = indice,
              values_from = valeur) %>%
  mutate(LAND_USE = factor(LAND_USE, levels = c("Urban", "Annual crops", "Permanent crops", "Pastures and meadows", 
                                           "Forests & tree plantations")))

# réprésentation SI -EI
ggplot(nem, aes(x=SI, y = EI, colour = LAND_USE))+
  geom_point(size=5)+
  lims(y = c(0, 100), x=c(0,100))+
  geom_hline(yintercept = 50)+
  geom_vline(xintercept = 50)+
  scale_colour_paletteer_d("lisa::FridaKahlo") +
  labs(x="Structure Index", y = "Enrichment Index",
       colour = "Land use")+
  coord_fixed(ratio = 1)+
  theme_bw()
