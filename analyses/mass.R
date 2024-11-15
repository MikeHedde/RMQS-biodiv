library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)

mass <- read.csv("C:/Users/heddemic/Documents/Codes R/RMQS/mass.csv", h=T, sep=";")

# mass par esp
mass_esp <- mass %>%
  select(NOM_VALIDE, MASSE) 

nb_mass <- mass_esp %>%
  group_by(NOM_VALIDE)%>%
  summarise(nb = n())%>%
  filter(nb>1)

mass_esp2 <- mass_esp%>%
  filter(NOM_VALIDE %in% nb_mass$NOM_VALIDE)

ggplot(mass_esp2, aes(x = reorder(NOM_VALIDE, MASSE, median), y = MASSE))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.5, width = 0.25, colour = "lightgreen")+
  labs(y = "Masse (mg)", x ="")+
  scale_y_continuous(trans='log10')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sp_list = c("Poecilus cupreus (Linnaeus, 1758)", 
            "Harpalus distinguendus (Duftschmid, 1812)",
            "Amara aenea (De Geer, 1774)", 
            "Carabus auratus Linnaeus, 1761",
            "Calathus fuscipes (Goeze, 1777)")

#masse en fonction de l'altitude
SP <- mass %>%
  filter(NOM_VALIDE %in% sp_list) %>%
  mutate(bin=cut_width(ALTITUDE, width=200, boundary=0))

ggplot(SP, aes(x = bin, y = MASSE))+
  geom_boxplot(outlier.shape = NA)+
  labs(y = "Masse (mg)", x ="Classe d'altitude (m)")+
  geom_jitter(alpha=0.5, width = 0.25, colour = "lightgreen")+
  facet_grid(NOM_VALIDE~., scales = "free_y")+
  theme_bw()


# abondance-mass scaling
mass_ab_esp <- mass %>%
  select(NOM_LOCALITE, NOM_VALIDE, MASSE, ABONDANCE_TOTALE) %>%
  group_by(NOM_LOCALITE, NOM_VALIDE) %>%
  summarise(mass= mean(MASSE), ab = sum(ABONDANCE_TOTALE))

ggplot(mass_ab_esp, aes(y = ab, x = mass))+
  geom_smooth(method = "lm", se = F, aes(colour = NOM_LOCALITE))+
  geom_density2d(color = "red") +
  geom_point(aes(color = NOM_LOCALITE), size = 3, alpha = 0.3)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
theme_bw()

# mass esp
ggplot(subset(mass_esp2, NOM_VALIDE %in% sp_list), aes(x = log10(MASSE), 
       fill = NOM_VALIDE))+
      geom_density(lwd = 1.2,
               linetype = 2,
               colour = 2, alpha = 0.1) +
  theme_bw()


# Rencsch rule: mass par esp par sex
rensch <- read.csv("C:/Users/heddemic/Documents/Codes R/RMQS/rensch.csv", h=T, sep=";")

rensch_wide <- rensch %>%
  group_by(station, hab_type, taxa, sex) %>%
  summarise(mean_mass = mean(mass)) %>% #calculer les masses moyennes par site, par esp et par sexe
  filter(!is.na(mean_mass)) %>%
  filter(all(c('male', 'female') %in% sex)) %>% # Vérifier que les deux sexes sont présents
  pivot_wider(names_from = sex, values_from = mean_mass)   # Élargir pour avoir une colonne pour chaque sexe

ggplot(rensch_wide, aes(x = log(male), y = log(female), group = hab_type)) +
  geom_point(aes(colour = hab_type)) +  # Ajouter des points
  geom_smooth(method = "glm", formula = y ~ poly(x, 2, raw=TRUE), se = FALSE, aes(colour = hab_type)) +  # Ajouter une régression linéaire
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Ajouter la bissectrice
  labs(x = "Log(Male Mass)", y = "Log(Female Mass)", title = "Rensch's Rule Representation") +
  #facet_wrap(~hab_type)+
  theme_minimal()

sd_wide <- rensch %>%
  mutate(genus = word(taxa, 1)) %>%
  group_by(station, hab_type, genus, taxa, sex) %>%
  summarise(sd_mass = sd(mass)) %>% #calculer les SD des masses par site, par esp et par sexe
  filter(!is.na(sd_mass)) %>%
  filter(all(c('male', 'female') %in% sex)) %>% # Vérifier que les deux sexes sont présents
  pivot_wider(names_from = sex, values_from = sd_mass) %>%  # Élargir pour avoir une colonne pour chaque sexe
  mutate(ratio = female/male-1) %>%
  group_by(hab_type, genus, taxa)%>%
  summarise(mean_sd = mean(ratio))

ggplot(sd_wide, aes(x=hab_type, y= reorder(taxa, mean_sd)))+
  geom_tile(aes(fill = mean_sd), color = "black")+
  labs(y = "Espèces", x ="Type d'habitat", fill = "Variabilité \nde la masse")+
  scale_fill_gradient2()+
  theme_bw()
