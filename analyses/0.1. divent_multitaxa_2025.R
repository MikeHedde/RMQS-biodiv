librarian::shelf(ggplot2, dplyr, forcats, ade4)


data_car <- read.csv("data/raw-data/lab files/databases update/carabidae_250906.csv", h = T, sep = ";", fileEncoding="latin1")
data_ara <- read.csv("data/raw-data/lab files/databases update/araneae_250906.csv", h = T, sep = ";", fileEncoding="latin1")


list_hab <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper) 

source("analyses/0. fonctions maison/DIVERSITY.R")
divent_car_ba <- DIVERSITY(data = data_car[data_car$METHODE == "Pitfall",])
divent_car_gpd <- DIVERSITY(data = data_car[data_car$METHODE == "GPD",])
divent_ara_ba <- DIVERSITY(data = data_ara[data_ara$METHODE == "Pitfall",])
divent_ara_gpd <- DIVERSITY(data = data_ara[data$ara$METHODE == "GPD",])


# Species accumulation curves
acc  <- rbind(cbind(taxa = "carabidae", METHOD = "Pitfall", divent_car_ba$acc), 
              cbind(taxa = "carabidae", METHOD = "GPD", divent_car_gpd$acc),
              cbind(taxa = "araneae", METHOD = "Pitfall", divent_ara_ba$acc) 
              #cbind(taxa = "araneae", METHOD = "GPD", divent_ara_gpd$acc)
              ) %>%
  left_join(divent_car_ba$station) 

      # représentation
      ggplot(acc, aes(x = level, y = diversity, group = STATION)) +
        geom_path(aes(colour = STATION))+
        facet_grid(METHOD ~ taxa)+
        theme_bw()
      

  