
sol0 <- read.csv("data/raw-data/pc_sols.csv", h = T, sep = ";") %>%
  select(STATION, PARAM, valeur) %>%
  pivot_wider(id_cols = STATION, names_from = PARAM, values_from = valeur)

sol_camp1 <- read.csv("data/raw-data/pc_sols_campagne1.csv", h = T, sep = ";")%>%
  filter(no_couche == 1) %>%
  mutate(sand = rowSums(.[16:17])) %>%
  select("id_site", "argile", "sand",
         "cd_tot_hf",	"cd_ext_66_1", "cu_tot_hf",	"cu_ext_66_1", "pb_tot_hf",	"pb_ext_66_1")%>%
  rename(STATION = id_site,
         clay = argile, 
         Cd_tot = cd_tot_hf,	
         Cd_disp = cd_ext_66_1, 
         Cu_tot = cu_tot_hf,	
         Cu_disp = cu_ext_66_1, 
         Pb_tot = pb_tot_hf,	
         Pb_disp = pb_ext_66_1
         )

sol <- sol0 %>%
  left_join(sol_camp1)

write.csv(sol, "data/derived-data/pc_sols.csv")
         