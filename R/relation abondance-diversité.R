# Chargement des packages nécessaires
install.packages("dplyr")
install.packages("vegan")
install.packages("ggplot2")
library(dplyr)
library(vegan)
library(ggplot2)
library(tibble)

# Chargement des données
sandbox <- read.csv("data/raw-data/sandbox.csv", h = T, sep = ";", encoding = "Latin-1")


mymat <- sandbox %>%
  filter(METHODE == "Pitfall") %>%
  filter(!ORDRE %in% c("VIDE", "MANQUANT")) %>%
  pivot_wider(id_cols = STATION, names_from = ORDRE, values_from = ABONDANCE_TOTALE, 
              values_fn = sum, values_fill = 0) 

ab <- rowSums(mymat[,-1])
div <- specnumber(mymat[,-1])
abdiv <- tibble(cbind(STATION = mymat$STATION, ab, div))

ggplot(abdiv, aes(x= ab, y = div))+
  geom_point()+
  labs(y = "Richesse spécifique", x = "Nombre d'individus collectés")+
  geom_smooth(method = lm, se = T)+
  theme_minimal()
