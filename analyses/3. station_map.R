# Représentation des points prélevés
# Auteur : Mickael Hedde
# Date : 17-03-2025
# Objectif : Cartographier les points d'échantillonnage 

# Chargement des bibliothèques nécessaires
librarian::shelf(ggplot2, dplyr, sf, osmdata, ggspatial, maptiles, tidyterra, giscoR)

# Définition de la palette de couleurs
land_use_colors <- c("Annual crops" = "#121510FF",
                     "Pastures and meadows" = "#6D8325FF",
                     "Urban" = "#D6CFB7FF",
                     "Permanent crops" = "#E5AD4FFF", 
                     "Forests & tree plantations" = "#BD5630FF")


# Chargement des données (exemple)
stations <- read.csv("data/raw-data/1.faune/nematoda.csv", sep = ";", h = T, dec=",") %>%
  mutate(LAND_USE = factor(LAND_USE, levels = names(land_use_colors))) %>%
  select(STATION, LAND_USE, X, Y) %>%
  rename(LATITUDE = Y,
         LONGITUDE = X)%>%
  unique()

# Transformation en objet spatial si nécessaire
stations_sf <- st_as_sf(stations, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)


# Définition de la bounding box (bbox) sous forme d'un objet sf
bbox_france <- st_bbox(c(xmin = -5.2, ymin = 41.0, xmax = 9.7, ymax = 51.5), crs = st_crs(4326))
 
# Conversion de bbox en un objet SpatVector (utilisable par get_tiles)
bbox_extent <- vect(st_as_sfc(bbox_france))
 
# Charger un fond de carte OpenStreetMap avec get_tiles
tiles <- get_tiles(x = bbox_extent)

# Charger les limites administratives de la France depuis giscoR
france_admin <- gisco_get_countries(resolution = "3", country = "FR")


# Création de la carte
map_plot <- 
  ggplot() +
  geom_sf(data = france_admin, color = "black", size = 0.5, fill = "grey90") +
  geom_point(data = stations, aes(x = LONGITUDE, y = LATITUDE, color = LAND_USE), size = 5) +
  coord_sf(xlim = c(-5.2, 9.7), ylim = c(41.0, 51.5)) +  # Recadrage sur la France
  scale_color_manual(name = "Land Use", values = land_use_colors) +
  theme_minimal() +
  theme(legend.position = "left")

# Sauvegarde
ggsave("figures/article ASE/stations_map.png", map_plot, width = 12, height = 8, units = "in", dpi = 300)

