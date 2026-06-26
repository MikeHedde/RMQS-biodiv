## Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, tidyverse, terra)

# Charger le raster du climat Köppen-Geiger
koppen_raster <- rast("data/raw-data/bioclimatic_variant.tif")  #https://doi.org/10.1016/j.dib.2020.105815

# Vérifier les propriétés du raster
print(koppen_raster)
plot(koppen_raster)  # Visualiser le raster

station_koppen <- read.csv("data/derived-data/all_env_variables.csv") %>%
   select(STATION, LAT, LON)

# Convertir en objet spatial (vecteur)
stations_sf <- vect(station_koppen, geom = c("LON", "LAT"), 
                    crs = crs(koppen_raster))

# Vérifier le système de projection
print(stations_sf)

# Extraire les valeurs du raster aux emplacements des stations
station_koppen$kopp <- extract(koppen_raster, stations_sf)[,2] 

station_koppen <- station_koppen%>% 
  unique()

# Vérifier les résultats
head(station_koppen)

write.csv(station_koppen, "data/derived-data/station_koppen.csv")
