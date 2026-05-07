# R/01_prep_all_clean.R
source("R/00_config.R")

prep_dat0 <- function(all_clean, cfg) {
  all_clean %>%
    filter(family == cfg$family_keep, method == cfg$method_keep) %>%
    mutate(
      site_id = paste(project, yr, city, station, sep = "_"),
      species = as.character(valid_name),
      yr      = as.integer(yr),
      rep     = as.integer(rep),
      det     = as.integer(tot_abund > 0)
    ) %>%
    filter(!is.na(site_id), !is.na(species), !is.na(yr), !is.na(rep))
}

# run
run_01 <- function(all_clean) {
  dat0 <- prep_dat0(all_clean, cfg)
  
  tag <- time_tag()
  f_out <- file.path(cfg$data_dir, "pitfall_segments", paste0("dat0_", tag, ".rds"))
  dir_create(path_dir(f_out))
  saveRDS(dat0, f_out)
  
  message("01_prep_all_clean OK: ", nrow(dat0), " rows -> ", f_out)
  invisible(list(dat0 = dat0, tag = tag, f_out = f_out))
}
