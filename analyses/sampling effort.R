# COURBES D ACCUMULATION D ESPECES

#Étapes :
  
# Transformation en matrice : Le code transforme d'abord les données en une matrice de présence/absence pour chaque site et chaque méthode.
# Boucle pour chaque site et méthode : Pour chaque combinacarabn, il calcule la courbe d'accumulation et ajuste un modèle asymptotique pour évaluer si un plateau est atteint.
# Boxplot : Un boxplot est généré pour montrer la proportion de sites ayant atteint le plateau pour chaque méthode.
    
    

# Chargement des packages nécessaires
install.packages("dplyr")
install.packages("tidyr")
install.packages("vegan")
install.packages("nls2")
library(dplyr)
library(tidyr)
library(vegan)
library(nls2)
library(ggplot2)

# Chargement des données
sandbox <- read.csv("data/raw-data/sandbox.csv", h = T, sep = ";", encoding = "Latin-1")
carab <- read.csv("data/raw-data/carabidae.csv", h = T, sep = ";", fileEncoding="latin1")


# 1. Calcul des courbes d'accumulation pour chaque méthode et chaque site
# Transformation en matrice de présence/absence pour chaque méthode
mymat <- carab %>%
  #filter(!ORDRE %in% c("VIDE", "MANQUANT")) %>%
  group_by(METHODE, STATION, REPETITION, LB_NOM) %>%
  summarise(ABONDANCE_TOTALE = sum(ABONDANCE_TOTALE), .groups = 'drop') %>%
  mutate(PRESENCE = ifelse(ABONDANCE_TOTALE > 0, 1, 0)) %>%
  select(-ABONDANCE_TOTALE) %>%
  pivot_wider(names_from = LB_NOM, values_from = PRESENCE, values_fill = 0) %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>%  # filtre les lignes (échantillons) à abondance nulle
  filter(STATION != c(136, 391,  392,  393, 633, 1180 ,1243, 1768, 1996, 2276 ))  #1 seul barber
  
# Calcul des courbes d'accumulation pour chaque méthode et pour tous les sites
site_accumulation_tm <- specaccum(mymat[mymat$METHODE=="Trimanuel", -c(1:3)], method = "random")
plot(site_accumulation_tm,  main = "SAC tri manuel")
site_accumulation_ba <- specaccum(mymat[mymat$METHODE=="Pitfall", -c(1:3)], method = "random")
plot(site_accumulation_ba,  main = "SAC Barber")

# Initialisation d'une liste pour stocker les résultats de courbes d'accumulation par méthode
accumulation_data <- list()

# Boucle pour calculer la courbe d'accumulation pour chaque méthode et chaque site
for (method in unique(mymat$METHODE)) {
  for (station in unique(mymat$STATION)) {
    # Sous-ensemble des données pour chaque méthode et chaque site
    site_data <- mymat %>%
      filter(METHODE == method, STATION == station) %>%
      select(-STATION, -REPETITION, -METHODE)
    
    # Calcul de la courbe d'accumulation pour le site et la méthode
    site_accumulation <- specaccum(site_data, method = "random")
    
    # Stockage des résultats dans un dataframe pour ggplot2
    accumulation_data[[paste(method, station, sep = "_")]] <- data.frame(
      Methode = method,
      Station = station,
      Sites = site_accumulation$sites,
      Richness = site_accumulation$richness
    )
  }
}

# Conversion de la liste en un seul dataframe
accumulation_df <- do.call(rbind, accumulation_data) %>%
  left_join(carab[, 1:9], join_by(Station==STATION))

# Graphique des courbes d'accumulation pour chaque méthode
ggplot(accumulation_df, aes(x = Sites, y = Richness, group = Station, colour = OBSERVATEUR)) +
  geom_line() +
  facet_wrap(REGION~., scales = "free_y") +
  labs(title = "Courbes d'Accumulation d'espèces par station",
       x = "Nombre d'échantillons",
       y = "Richesse spécifique",
       color = "Site") +
  theme_minimal() +
  theme(legend.position = "none")  # Supprime la légende si trop de sites



# Initialisation de la liste pour stocker les résultats du test de plateau pour chaque site et chaque méthode
plateau_status <- data.frame(STATION = unique(mymat$STATION), METHODE = unique(mymat$METHODE), PLATEAU = 0)

# 2. Boucle pour calculer la courbe d'accumulation et tester le plateau pour chaque méthode et chaque site
for (method in unique(mymat$METHODE)) {
  for (station in unique(mymat$STATION)) {
    # Sous-ensemble des données pour chaque méthode et chaque site
    site_data <- mymat %>% 
      filter(METHODE == method, STATION == station) %>%
      select(-STATION, -REPETITION, -METHODE)
    
    # Calcul de la courbe d'accumulation pour le site et la méthode
    site_accumulation <- specaccum(site_data, method = "random")
    n_samples <- site_accumulation$sites
    richness <- site_accumulation$richness
    
    # Ajustement d'un modèle asymptotique simple (Michaelis-Menten) pour estimer la richesse maximale
    model <- try(nls(richness ~ a * n_samples / (b + n_samples),
                     start = list(a = max(richness), b = max(n_samples) / 2)),
                 silent = TRUE)
    
    # Vérification si l'ajustement du modèle a réussi
    if (inherits(model, "nls")) {
      # Extraction de la richesse asymptotique estimée
      plateau <- coef(model)["a"]
      
      # Vérification si la courbe atteint 95 % du plateau estimé
      plateau_status$PLATEAU[plateau_status$STATION == station & plateau_status$METHODE == method] <- 
        ifelse(max(richness) >= 0.9 * plateau, 1, 0)
    } else {
      # Si l'ajustement du modèle a échoué, considérez comme non atteint
      plateau_status$PLATEAU[plateau_status$STATION == station & plateau_status$METHODE == method] <- 0
    }
  }
}

# 3. Création d'un boxplot pour visualiser les résultats du test de plateau par méthode
ggplot(plateau_status, aes(x = METHODE, y = PLATEAU)) +
  geom_jitter(width = 0.25, height =0) +
  labs(title = "Atteinte du Plateau d'Accumulation d'Espèces \npar Méthode et par Site (90% proba)",
       y = "Plateau atteint \n(1 = Oui, 0 = Non)",
       x = "Méthode") +
  theme_minimal()
