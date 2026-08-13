# Master file
# Replication of Tables and Figures

directory <- "." # change working directory here

dir_R       <- file.path(directory, "R")
dir_data    <- file.path(directory, "data")
dir_network <- file.path(dir_data, "networks")
dir_out    <- file.path(directory, "output")
dir_fig     <- file.path(dir_out, "figures")
dir_table   <- file.path(dir_out, "tables")

for (d in c(dir_network, dir_fig, dir_table))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

results_NASC     <- file.path(dir_data, "data_ex_rho_2.rds")
results_SCM <- file.path(dir_data, "data_ex_rho_scm.rds")

dgp_path <- file.path(dir_R, "dgp_functions.R")
source(dgp_path)

for (p in c(0, 0.1, 0.4)) {
  f <- file.path(dir_network, sprintf("SAR_k2_p%.2f_N15.rds", p))
  if (!file.exists(f))
    saveRDS(generate_watts_strogatz_matrix(N = 15, k = 2, p = p, seed = 13), f)
}

PLOT_FONT <- "Times"
N_WORKERS <- max(1L, min(12L, parallel::detectCores() - 1L))

writeLines(capture.output(sessionInfo()), file.path(dir_out, "sessionInfo.txt"))

run_all <- function() {
  scripts <- c(
    "Figures/Plot_Ternary_ETDir.R",
    "Figures/Plot_MC_bias.R",
    "Figures/Plot_MC_SC_contam.R",
    "Figures/Plot_MC_bias_corr_factor.R",
    "Figures/Plot_example_MC_iter.R",
    "Tables/Table_MC_results.R",
    "Tables/Table_MC_SC_contam.R"
  )
  for (s in scripts) source(file.path(dir_R, s), local = new.env())
  invisible(NULL)
}
