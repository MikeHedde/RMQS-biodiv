librarian::shelf(dplyr, ggplot2, MODISTools)

# Définir les points d'intérêt
ROIs <- read.csv("data/raw-data/sampling_date.csv", h = T, sep = ";") %>%
  rename(site_name = station, lat = Latitude, lon = Longitude) %>%
  select(site_name, lat, lon) %>%
  mutate(site_name = as.factor(site_name))

# Sélectionner le/les produits MODIS d'intérêt
list_products <- mt_products()[c(5, 6, 7, 8, 16),]$product
MODproduct = list_products[5]

# Sélectionner la plage de dates
NDVI_dates <- mt_dates(product = MODproduct,
                       lat = ROIs$lat[1],
                       lon = ROIs$lon[1])
startDate = "2024-01-01"
endDate = "2024-08-12"

#Selectionner la 'band' qui nous interesse (possibilité d'en sélectionner plusieurs)
bands = mt_bands(MODproduct)
bandsOfInterest = bands[5,]$band

# Téléchargement en batch
## NDVI data
NDVI = mt_batch_subset(df = ROIs,
                       product = MODproduct,
                       band = bandsOfInterest,
                       internal = T,
                       start = startDate,
                       end = endDate)
length(unique(NDVI$site))

## Graph de la dynamique temporelle du NDVI
hab <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  mutate(STATION = as.factor(STATION))

NDVI_hab <- NDVI %>%
  left_join(hab, join_by(site == STATION))
  
ggplot(NDVI_hab, aes(x = as.Date(calendar_date), y = value/10000, colour = HABITAT_CODE))+
  geom_line()+
  labs(x = "date", y = "NDVI")+
  facet_wrap(.~site)+
  theme_bw()

## Qualité du produit VI
VI_q = mt_batch_subset(df = ROIs,
                       product = MODproduct,
                       band = "250m_16_days_pixel_reliability",
                       internal = T,
                       start = startDate,
                       end = endDate)

## Binarization de la variable de qualité (0 = bonne qualité)
VI_q2 <- VI_q %>%
  rename(quality = value) %>%
  mutate(quality = case_when(quality %in% c(-1, 0) ~ 0,
                             quality %in% c(1, 2, 3) ~ 1)) %>%
  select(site, start, end, calendar_date, quality)

# Fusion des 2 data frames
NDVI_q <- NDVI_hab %>%
  rename(ndvi = value) %>%
  mutate(ndvi = ndvi/10000) %>%
  filter(ndvi>0) %>%
  select(site, start, end, calendar_date, HABITAT_CODE, HABITAT_TYPE, ndvi) %>%
  left_join(VI_q2)

## Graph de la dynamique temporelle du NDVI et de sa qualité
ggplot(subset(NDVI_q, quality == 0), aes(x = as.Date(calendar_date), y = ndvi, colour = HABITAT_CODE))+
  geom_line()+
  labs(x = "date", y = "NDVI (cleaned)")+
  facet_wrap(.~site)+
  theme_bw()

## Calcul de la moyenne, du maximum et de la déviation standard
NDVI_short <- NDVI_q %>%
  select(site, HABITAT_CODE, HABITAT_TYPE, ndvi) %>%
  group_by(site, HABITAT_CODE, HABITAT_TYPE) %>%
  summarise(ndvi_m = mean(ndvi), ndvi_max = max(ndvi), ndvi_sd = sd(ndvi))

ggplot(NDVI_short, aes(x = HABITAT_TYPE, y = ndvi_m))+
  geom_boxplot()+
  geom_jitter(aes(colour = HABITAT_CODE))+
  labs(x = "habitat", y = "NDVI (cleaned)")+
  theme_bw()

# Export données
write.csv(NDVI_short, "data/derived-data/ndvi.csv")
