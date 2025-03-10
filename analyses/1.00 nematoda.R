# Analyse des indicateurs nématofauniques
# Auteur : Mickael Hedde
# Date : 10-03-2025
# Objectif : Décrire la valeurs des indicateurs basés sur la nématofaune 

# Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, tidyr, gdm, forcats, stringr, ggpubr, paletteer)

# Chargement des données de nematofaune
nem0 <- read.csv("data/raw-data/1.faune/nematoda.csv", h = T, sep = ";")  %>%
  mutate(LAND_USE = factor(LAND_USE, levels = c("Forests & tree plantations", "Pastures and meadows", 
                                                "Annual crops", "Permanent crops", "Urban" ),
                           ordered = T))

# réprésentation SI -EI
nem_indice <- nem0 %>%
  pivot_wider(id_cols = c(STATION, LAND_USE, ALT, X), names_from = indice,
              values_from = valeur) 

indices_plot <- ggplot(nem_indice, aes(x=SI, y = EI, colour = str_wrap(LAND_USE, 15)))+
  geom_point(size=5)+
  lims(y = c(0, 100), x=c(0,100))+
  geom_hline(yintercept = 50)+
  geom_vline(xintercept = 50)+
  scale_colour_paletteer_d("lisa::FridaKahlo") +
  labs(x="Structure Index", y = "Enrichment Index",
       colour = "Land use")+
  coord_fixed(ratio = 1)+
  theme_bw()+
  theme(legend.key.height=unit(1, "cm")) 


# réprésentation %phyto
nem_ab<- nem0 %>%
  pivot_wider(id_cols = c(STATION, LAND_USE), names_from = indice,
              values_from = valeur) %>%
  
  mutate(prc_phyto = ab_phytofac/ab_tot*100) %>%
  select(STATION, LAND_USE, ab_libres, ab_phytofac, ab_phytopar, prc_phyto) %>%
  pivot_longer(cols=3:6)

indice_names <- as_labeller(c("ab_libres" = "Free-living nematode abundance \n(Individuals per 100 g dry soil)",
                     "ab_phytopar" = "Obligate plant-feeding nematode abundance \n(Individuals per 100 g dry soil)",
                     "ab_phytofac" = "Facultative plant-feeding  nematode abundance \n(Individuals per 100 g dry soil)",
                     "prc_phyto" = "Facultative plant-feeding  nematodes proportion \n(Percentage of the total community)"))  
  
ab_plot <- ggplot(nem_ab, aes(x=LAND_USE, y= value, fill = LAND_USE)) +
  geom_boxplot(outliers = FALSE) + 
  scale_fill_paletteer_d("lisa::FridaKahlo", direction = 1)+
  labs(x="", y="")+
  facet_wrap(name~., labeller = indice_names,
             dir = "v", scales = "free")+
  theme_bw()+
  theme(legend.position ="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

assembled_nem <- ggarrange(indices_plot, ab_plot,
                           labels = c("A", "B"),
                           ncol = 2, nrow = 1)    
  
ggsave("figures/article ASE/nematodes.png", 
       width = 21, height = 12, units = "cm")
