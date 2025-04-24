# Analyse des indicateurs nématofauniques
# Auteur : Mickael Hedde
# Date : 17-03-2025
# Objectif : Décrire la valeurs des indicateurs basés sur la nématofaune 

# Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, tidyr, gdm, forcats, stringr, ggpubr, paletteer)

# Chargement des données de nematofaune
nem0 <- read.csv("data/raw-data/1.faune/nematoda.csv", h = T, sep = ";")  %>%
  mutate(LAND_USE = factor(LAND_USE, levels = c("Forests & tree plantations", "Pastures and meadows", 
                                                "Annual crops", "Permanent crops", "Urban" ),
                           ordered = T))%>%
  mutate(fact_col = case_when(LAND_USE ==  "Annual crops" ~ "#121510FF",
                              LAND_USE ==  "Pastures and meadows" ~ "#6D8325FF",
                              LAND_USE ==  "Urban" ~ "#D6CFB7FF",
                              LAND_USE ==  "Permanent crops" ~ "#E5AD4FFF", 
                              LAND_USE ==  "Forests & tree plantations" ~ "#BD5630FF"))
# Matrice d'indices nématologiques
nem_indice <- nem0 %>%
  pivot_wider(id_cols = c(STATION, LAND_USE, ALT, X, fact_col), names_from = indice,
              values_from = valeur) 

# Matrice ddes abondances
nem_ab<- nem0 %>%
  pivot_wider(id_cols = c(STATION, LAND_USE, fact_col), names_from = indice,
              values_from = valeur) %>%
  filter(!LAND_USE %in% c("Permanent crops", "Urban")) %>%
  mutate(prc_phyto = ab_phytofac/ab_tot*100) %>%
  select(STATION, LAND_USE, fact_col, ab_libres, ab_phytofac, ab_phytopar, prc_phyto) %>%
  pivot_longer(cols=4:7)


# Définition des couleurs pour chaque LAND_USE
land_use_colors <- c("Annual crops" = "#121510FF",
                     "Pastures and meadows" = "#6D8325FF",
                     "Urban" = "#D6CFB7FF",
                     "Permanent crops" = "#E5AD4FFF", 
                     "Forests & tree plantations" = "#BD5630FF")

indice_names <- as_labeller(c("ab_libres" = "Free-living nematode abundance \n(Individuals per 100 g dry soil)",
                              "ab_phytopar" = "Obligate plant-feeding nematode abundance \n(Individuals per 100 g dry soil)",
                              "ab_phytofac" = "Facultative plant-feeding  nematode abundance \n(Individuals per 100 g dry soil)",
                              "prc_phyto" = "Facultative plant-feeding  nematodes proportion \n(Percentage of the total community)"))  

# Réprésentation SI - EI
indices_plot <- ggplot(nem_indice, aes(x=SI, y=EI, colour=LAND_USE)) +
  geom_point(size=5) +
  scale_x_continuous(limits=c(0, 100), expand=c(0, 0)) +  
  scale_y_continuous(limits=c(0, 100), expand=c(0, 0)) +
  geom_hline(yintercept=50) +
  geom_vline(xintercept=50) +
  scale_color_manual(name="Land use", values=land_use_colors) +  # Utilisation du même code couleur
  labs(x="Structure Index", y="Enrichment Index") +
  coord_fixed(ratio=1) +
  theme_bw() +
  theme(legend.key.height=unit(1, "cm"))

# Réprésentation %phyto (abondance des groupes trophiques)
ab_plot <- ggplot(nem_ab, aes(x=LAND_USE, y=value, fill=LAND_USE)) +
  geom_boxplot() + 
  labs(x="", y="") +
  facet_wrap(name~., labeller=indice_names, dir="v", scales="free") +
  scale_fill_manual(name="Land use", values=land_use_colors) +  # Applique les mêmes couleurs que dans indices_plot
  theme_bw() +
  theme(legend.position ="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

# Assemblage des figures
assembled_nem <- ggarrange(indices_plot, ab_plot,
                           labels = c("A", "B"),
                           ncol = 2, nrow = 1)

# Sauvegarde de la figure finale
ggsave("figures/article ASE/nematodes.png", 
       width = 36, height = 12, units = "cm")
