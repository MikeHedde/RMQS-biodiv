# R/02_build_pitfall_segments.R
source("R/00_config.R")

build_segments <- function(dat0, cfg) {
  
  # reps réellement présents = effort réalisé
  reps_realized <- dat0 %>%
    distinct(site_id, yr, rep)
  
  # det par rep
  det_rep <- dat0 %>%
    group_by(site_id, yr, rep, species) %>%
    summarise(det = as.integer(any(det == 1)), .groups = "drop")
  
  # compléter les 0 sur reps réalisés
  site_species <- dat0 %>% distinct(site_id, yr, species)
  
  det_rep_full <- reps_realized %>%
    inner_join(site_species, by = c("site_id","yr")) %>%
    left_join(det_rep, by = c("site_id","yr","rep","species")) %>%
    mutate(det = replace_na(det, 0L))
  
  # fonction segment: k premiers reps présents (triés par rep)
  seg_one <- function(k) {
    det_rep_full %>%
      group_by(site_id, yr, species) %>%
      group_modify(~{
        df <- .x %>% arrange(rep)
        n  <- nrow(df)  # nb reps réalisés (constant pour site×yr)
        tibble(det = if (n < k) NA_integer_ else as.integer(any(df$det[seq_len(k)] == 1)))
      }) %>%
      ungroup()
  }
  
  out6  <- seg_one(6L)  %>% mutate(seg = "Pitfall6",  effort = 6L)
  out10 <- seg_one(10L) %>% mutate(seg = "Pitfall10", effort = 10L)
  
  det_long <- bind_rows(out6, out10) %>%
    filter(!(seg == "Pitfall10" & yr == 2024L)) # optionnel : explicite
  
  det_long
}

run_02 <- function(dat0, tag) {
  det_long <- build_segments(dat0, cfg)
  
  f_csv <- file.path(cfg$data_dir, "pitfall_segments", paste0("det_segments_", tag, ".csv"))
  f_rds <- file.path(cfg$data_dir, "pitfall_segments", paste0("det_segments_", tag, ".rds"))
  dir_create(path_dir(f_csv))
  write_csv(det_long, f_csv)
  saveRDS(det_long, f_rds)
  
  message("02_build_pitfall_segments OK: ", nrow(det_long), " rows -> ", f_csv)
  invisible(list(det_long = det_long, f_csv = f_csv, f_rds = f_rds))
}
