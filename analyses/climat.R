library(nasapower)
library(lubridate)
library(dplyr)
library(zoo)

# Importer les dates de prélèvement
sampl_date0 <- read.csv("data/raw-data/sampling_date.csv", h=T, sep = ";")
sampl_date <- sampl_date0 %>%
  mutate(date = as.Date(date, "%d/%m/%Y"))

# Initialiser une liste pour stocker les résultats
climate_data <- list()

# Boucle pour récupérer les données climatiques pour chaque station
for (i in 1:nrow(sampl_date)) { 

  # Récupérer les données pour une station
    data0 <- get_power(
      community = "ag",
      lonlat = c(sampl_date$Latitude[i], sampl_date$Longitude[i]), 
      pars = c("T2M_MAX", "T2M_MIN", "T2MDEW", "EVLAND"), #liste en bas de https://power.larc.nasa.gov/      
      dates = c("2024-01-01", "2024-08-31"), 
      temporal_api = "daily"
    )
    
    # Ajouter une colonne pour identifier la station
    clim_data <- data0 %>%
      mutate(STATION = sampl_date$station[i]) %>%
      as.data.frame()
    
    # Stocker dans la liste
    climate_data[[i]] <- clim_data
}

# Combiner toutes les données en une seule tibble
climate_data <- bind_rows(climate_data) 

# enregistrement
write.csv(x=final_data, file="data/derived-data/final_data.csv")

# calculer les degré-jours
    # Température de base
      T_base <- 10

      DD <- climate_data %>%
    # Calcul des températures moyennes
      mutate(T_avg = (climate_data$T2M_MAX + climate_data$T2M_MIN) / 2) %>%
    # Supprime les NA
      filter(!is.na(T_avg)) %>%
    # Calcul des degrés-jours
      group_by(STATION) %>%
      mutate(degree_days = pmax(0, T_avg - T_base)) %>%
    # Somme cumulative des degrés-jours
      mutate(DD_cum = cumsum(degree_days)) %>%
      ungroup() %>%
    # fusion avec les dates de prélèvements
      right_join(sampl_date[,1:2], by = join_by(STATION == station, YYYYMMDD == date)) %>%
      select(STATION,DD_cum)
      
# calculer le point de rosée moyen et l'évapotranspiration sur les 30 derniers jours 
      DEW_EV <- climate_data %>%
    # Supprime les NA
      filter(!is.na(T2MDEW)) %>%
    # Calcul du point de rosée
      group_by(STATION) %>%
      mutate(T2MDEW_moy = rollapply(T2MDEW, width = 30, 
                                    FUN = mean, align = "right", fill = NA),
             EVLAND_moy = rollapply(EVLAND, width = 30, 
                                    FUN = mean, align = "right", fill = NA),
             T2max_moy = rollapply(T2M_MAX, width = 30, 
                                    FUN = mean, align = "right", fill = NA)) %>%
      ungroup() %>%
    # fusion avec les dates de prélèvements
      right_join(sampl_date[,1:2], by = join_by(STATION == station, YYYYMMDD == date)) %>%
      select(STATION,T2MDEW_moy,EVLAND_moy, T2max_moy)

# grouper DD, EV et DEW
      clim_moy <- DD %>%
        left_join(DEW_EV)
       
# enregistrement
write.csv(x=clim_moy, file="data/derived-data/clim_moy.csv")




