# ============================================================
# 004_resume_taxon_profiles_after_fix.R
# Resume after the metric-column repair without rerunning simulations
# ============================================================
# Run this in the SAME R session that produced `experiment`.
# It recomputes summaries and figures only; it does not simulate worlds again.

if (!exists("experiment", inherits = TRUE)) {
  stop(
    "Object `experiment` is not available. Start a fresh R session and rerun 003_run_taxon_profiles.R after replacing 13_taxon_specific_simulation.R with the corrected version.",
    call. = FALSE
  )
}
if (!exists("summarise_taxon_profile_experiment", inherits = TRUE)) {
  stop("Source 13_taxon_specific_simulation.R first.", call. = FALSE)
}

summary_out <- summarise_taxon_profile_experiment(experiment)

saveRDS(experiment, file.path(THEORY_OUT_DIR, "taxon_profile_experiment.rds"))
write_theory_csv(experiment$profiles, "taxon_profile_simulated_worlds.csv")
write_theory_csv(experiment$results, "taxon_profile_simulation_by_replicate.csv")
write_theory_csv(summary_out$blowes_summary, "taxon_profile_blowes_summary.csv")
write_theory_csv(summary_out$stability_summary, "taxon_profile_stability_summary.csv")

fig_gamma <- plot_cross_taxon_scenarios(summary_out$blowes_summary, response = "delta_gamma")
fig_occupancy <- plot_cross_taxon_scenarios(summary_out$blowes_summary, response = "delta_occupancy")
fig_blowes <- plot_cross_taxon_blowes(summary_out$blowes_summary)

ggsave(file.path(THEORY_OUT_DIR, "FigP2_cross_taxon_delta_gamma.pdf"), fig_gamma, width = 10, height = 6, units = "in")
ggsave(file.path(THEORY_OUT_DIR, "FigP3_cross_taxon_delta_occupancy.pdf"), fig_occupancy, width = 10, height = 6, units = "in")
ggsave(file.path(THEORY_OUT_DIR, "FigP4_cross_taxon_blowes.pdf"), fig_blowes, width = 9, height = 7, units = "in")

ggsave(file.path(THEORY_OUT_DIR, "FigP2_cross_taxon_delta_gamma.png"), fig_gamma, width = 10, height = 6, units = "in", dpi = 300)
ggsave(file.path(THEORY_OUT_DIR, "FigP3_cross_taxon_delta_occupancy.png"), fig_occupancy, width = 10, height = 6, units = "in", dpi = 300)
ggsave(file.path(THEORY_OUT_DIR, "FigP4_cross_taxon_blowes.png"), fig_blowes, width = 9, height = 7, units = "in", dpi = 300)

message("Taxon-profile summaries and figures written to: ", THEORY_OUT_DIR)

