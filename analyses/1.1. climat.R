library(nasapower)
library(lubridate)
library(dplyr)
library(zoo)

# Importer les dates de prélèvement
sampl_date0 <- read.csv("data/raw-data/sampling_date.csv", h=T, sep = ";")
sampl_date <- sampl_date0 %>%
  mutate(date = as.Date(date, "%d/%m/%Y"))%>%
  rename(sampl_date = date, LAT = Latitude, 
         LON = Longitude, STATION = station)


# Initialiser une liste pour stocker les résultats
climate_data <- list()

# Boucle pour récupérer les données climatiques pour chaque station
for (i in 1:nrow(sampl_date)) { 

  # Récupérer les données pour une station
  climate_data[[i]] <- get_power(
      community = "ag",
      lonlat = c(sampl_date$LON[i], sampl_date$LAT[i]), 
      pars = c("T2M_MAX", "T2M_MIN"), #liste en bas de https://power.larc.nasa.gov/      
      dates = c("2023-01-01", "2024-08-31"), 
      temporal_api = "daily") %>%
      as.data.frame() %>%
      mutate(STATION = sampl_date$STATION[i])
}

# Combiner toutes les données en une seule tibble
clim_data <- bind_rows(lst(climate_data)) %>%
  left_join(sampl_date) %>%
  mutate(T2M_MEAN = (T2M_MAX + T2M_MIN) / 2)

# enregistrement
write.csv(x=clim_data, file="data/derived-data/climate_data.csv")


#################
clim_data <- read.csv("data/derived-data/climate_data.csv", h=T, sep = ",")
# calculer les degré-jours
    # Température de base
      T_base <- 10

      DD <- clim_data %>%
    # Supprime les NA
      filter(!is.na(T2M_MEAN),
             YYYYMMDD > as.Date("2024-01-01"),
             YYYYMMDD < sampl_date) %>%
    # Calcul des degrés-jours
      group_by(STATION) %>%
      mutate(degree_days = pmax(0, T2M_MEAN - T_base)) %>%
    # Somme cumulative des degrés-jours
      mutate(DD_cum = cumsum(degree_days)) %>%
      filter(DD_cum == max(DD_cum)) 
      
# calculer les Températures sur les 30 derniers jours 
      Temp <- clim_data %>%
    # Calcul du point de rosée
      group_by(STATION) %>%
      mutate(T30_mean = rollapply(T2M_MEAN, width = 30, 
                                    FUN = mean, align = "right", fill = NA),
             T30_max = rollapply(T2M_MAX, width = 30, 
                                    FUN = mean, align = "right", fill = NA),
             T30_min = rollapply(T2M_MIN, width = 30, 
                                    FUN = mean, align = "right", fill = NA),
             T30_sd = rollapply(T2M_MEAN, width = 90, 
                                 FUN = sd, align = "right", fill = NA),
             T360_sd = rollapply(T2M_MEAN, width = 360, 
                                FUN = sd, align = "right", fill = NA),
             T360_mean = rollapply(T2M_MEAN, width = 360, 
                                  FUN = mean, align = "right", fill = NA),) %>%
      ungroup() %>%
        filter(YYYYMMDD == sampl_date) %>%
        select(STATION, T30_mean, T30_max, T30_min, T30_sd, T360_sd, T360_mean)

# grouper DD et Temp
      clim_stat <- DD %>%
        left_join(Temp)
       
# enregistrement
write.csv(x=clim_stat, file="data/derived-data/clim_stat.csv")




