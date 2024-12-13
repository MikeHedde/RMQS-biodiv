
data_car <- read.csv("data/raw-data/carabidae.csv", h = T, sep = ";", fileEncoding="latin1")
data_ara <- read.csv("data/raw-data/araneae.csv", h = T, sep = ";", fileEncoding="latin1")
data_iso <- read.csv("data/raw-data/isopoda.csv", h = T, sep = ";", fileEncoding="latin1")
data_diplo <- read.csv("data/raw-data/diplopoda.csv", h = T, sep = ";", fileEncoding="latin1")
data_coll <- read.csv("data/raw-data/collembola.csv", h = T, sep = ";", fileEncoding="latin1")

list_hab <- read.csv("data/derived-data/liste_habitat.csv", h = T, sep = ";") %>%
  rename_with(toupper) 

source("analyses/0. fonctions maison/DIVERSITY.R")
divent_car <- DIVERSITY(data = data_car)
divent_ara <- DIVERSITY(data = data_ara)
divent_iso <- DIVERSITY(data = data_iso)

autoplot(divent_ara$acc)+
  theme_bw()


