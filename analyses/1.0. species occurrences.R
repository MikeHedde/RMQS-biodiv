librarian::shelf(ggplot2, dplyr, tidyr, forcats, ade4, paletteer)


data_car <- read.csv("data/raw-data/1.faune/carabidae.csv", h = T, sep = ";", fileEncoding="latin1") %>%
  separate(ID_ECHANTILLON, "_", into = c("pj", "YEAR", "CITY", "STATION", "ECH")) %>%
  select(STATION, NOM_VALIDE, ABONDANCE_TOTALE, RANK) %>%
  mutate(STATION = as.integer(STATION), taxa = "Carabidae")
data_ara <- read.csv("data/raw-data/1.faune/araneae.csv", h = T, sep = ";", fileEncoding="latin1")  %>%
  select(STATION, NOM_VALIDE, ABONDANCE_TOTALE, RANK) %>%
  mutate(taxa = "Araneae")
data_iso <- read.csv("data/raw-data/1.faune/isopoda.csv", h = T, sep = ";", fileEncoding="latin1")  %>%
  select(STATION, NOM_VALIDE, ABONDANCE_TOTALE, RANK) %>%
  mutate(taxa = "Isopoda")
data_diplo <- read.csv("data/raw-data/1.faune/diplopoda.csv", h = T, sep = ";", fileEncoding="latin1")   %>%
  select(STATION, NOM_VALIDE, ABONDANCE_TOTALE, RANK) %>%
  mutate(taxa = "Diplopoda")
data_coll <- read.csv("data/raw-data/1.faune/collembola.csv", h = T, sep = ";", fileEncoding="latin1")  %>%
  select(STATION, NOM_VALIDE, ABONDANCE_TOTALE, RANK) %>%
  mutate(taxa = "Collembola")
data_form <- read.csv("data/raw-data/1.faune/formicidae.csv", h = T, sep = ";", fileEncoding="latin1")  %>%
  select(STATION, NOM_VALIDE, ABONDANCE_TOTALE, RANK) %>%
  mutate(taxa = "Formicidae")
data_vdt <- read.csv("data/raw-data/1.faune/vdt.csv", h = T, sep = ";", fileEncoding="latin1")  %>%
  select(STATION, NOM_VALIDE, ABONDANCE_TOTALE, RANK) %>%
  mutate(taxa = "Oligochaeta")

occ_sp <- bind_rows(list(data_car, data_ara, data_iso, data_diplo, data_coll, data_form, data_vdt)) %>%
  filter(RANK == "S") %>%
  group_by(taxa, NOM_VALIDE, STATION) %>%
  summarise(ab = sum(ABONDANCE_TOTALE)) %>%
  mutate(occ = ifelse(ab == 0, 0, 1))

nb_station = length(unique(occ_sp$STATION))

# calcul de l'occurence des espèces dans tout le dataset
occ_prop_sp <- occ_sp %>%
  group_by(taxa, NOM_VALIDE) %>%
  summarise(prop = sum(occ)/nb_station *100) %>%
  arrange(prop)

# Graphique
occ_prop_plot <- ggplot(occ_prop_sp, aes(x = reorder(NOM_VALIDE, -prop), y = prop)) + 
  geom_bar(stat = "identity", aes(fill = taxa))+
  facet_wrap(taxa~. , nrow = 2)+
  geom_hline(yintercept = 20, linetype = "dotted")+
  labs(x = "", y = "Precentage of occurrence")+
  scale_fill_paletteer_d("rcartocolor::Temps")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.length.x = unit(0,"cm"),
        legend.position="none")

ggsave("figures/article ASE/occurrence_proportion.png",
         width = 19, height = 10, units = "cm")
