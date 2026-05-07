# R/03_fit_occu_by_species.R
source("R/00_config.R")

fit_one_species <- function(df_sp, cfg) {
  sp <- unique(df_sp$species)
  
  # wide sites × seg
  Yw <- df_sp %>%
    select(site_id, seg, det) %>%
    pivot_wider(names_from = seg, values_from = det) %>%
    arrange(site_id)
  
  segs <- intersect(c("Pitfall6","Pitfall10"), names(Yw))
  if (length(segs) < 1) return(NULL)
  
  y <- as.matrix(Yw[, segs, drop = FALSE])
  
  # filtre minimal : assez de sites non-NA
  keep_sites <- rowSums(!is.na(y)) > 0
  y <- y[keep_sites, , drop = FALSE]
  if (nrow(y) < cfg$min_sites) return(NULL)
  
  # si jamais aucune détection
  if (sum(y == 1, na.rm = TRUE) == 0) return(NULL)
  
  # obsCov seg (facteur) : matrix nsites × nvisits
  seg_levels <- colnames(y)
  seg_mat <- matrix(rep(seg_levels, each = nrow(y)), nrow = nrow(y))
  seg_mat <- apply(seg_mat, 2, factor, levels = c("Pitfall6","Pitfall10"))
  
  umf <- unmarked::unmarkedFrameOccu(
    y = y,
    obsCovs = list(seg = seg_mat)
  )
  
  fit <- tryCatch(
    unmarked::occu(detformula = ~ seg, stateformula = ~ 1, data = umf),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  # predict p for each seg level
  newdata <- data.frame(seg = factor(c("Pitfall6","Pitfall10"), levels = c("Pitfall6","Pitfall10")))
  pr <- predict(fit, type = "det", newdata = newdata)
  
  ptab <- tibble(
    species = sp,
    seg = c("Pitfall6","Pitfall10"),
    p_hat = pr$Predicted,
    lcl   = pr$lower,
    ucl   = pr$upper
  )
  
  list(fit = fit, ptab = ptab)
}

run_03 <- function(det_long, tag) {
  species_list <- det_long %>% distinct(species) %>% pull(species)
  
  res <- tibble(species = species_list) %>%
    mutate(out = map(species, ~{
      df_sp <- det_long %>% filter(species == .x)
      fit_one_species(df_sp, cfg)
    }))
  
  tab_p <- res %>%
    filter(!map_lgl(out, is.null)) %>%
    transmute(species, ptab = map(out, "ptab")) %>%
    unnest(ptab)
  
  fits <- res %>%
    filter(!map_lgl(out, is.null)) %>%
    transmute(species, fit = map(out, "fit"))
  
  f_p <- file.path(cfg$data_dir, "p_hat", paste0("p_hat_pitfall6_vs_10_", tag, ".csv"))
  f_fit <- file.path(cfg$data_dir, "occ_fits", paste0("fits_pitfall6_vs_10_", tag, ".rds"))
  dir_create(path_dir(f_p)); dir_create(path_dir(f_fit))
  write_csv(tab_p, f_p)
  saveRDS(fits, f_fit)
  
  message("03_fit_occu_by_species OK: ", nrow(tab_p), " rows -> ", f_p)
  invisible(list(tab_p = tab_p, fits = fits, f_p = f_p, f_fit = f_fit))
}
