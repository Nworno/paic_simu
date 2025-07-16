
source("env_variables.R")
weighting_plots_files <- list.files(file.path("results_simulations", DATE_EXPERIMENT),
                                    pattern = c("weighting"),
                                    recursive = TRUE,
                                    full.names = TRUE)
processed_results_files <- list.files(file.path("results_simulations", DATE_EXPERIMENT, "processed_results"),
                                      recursive = TRUE,
                                      full.names = TRUE)
average_outcome_files <- list.files(file.path("results_simulations", DATE_EXPERIMENT),
                                    pattern = "average_outcome_df",
                                    recursive = TRUE,
                                    full.names = TRUE)
parameters_files <- list.files(file.path("results_simulations", DATE_EXPERIMENT),
                               pattern = "parameters",
                               recursive = FALSE,
                               full.names = TRUE)
pngs <- list.files(file.path("results_simulations", DATE_EXPERIMENT),
                   pattern = "\\.png",
                   recursive = TRUE,
                   full.names = TRUE)
pdfs <- list.files(file.path("results_simulations", DATE_EXPERIMENT),
                   pattern = "\\.pdf",
                   recursive = TRUE,
                   full.names = TRUE)
jpegs <- list.files(file.path("results_simulations", DATE_EXPERIMENT),
                   pattern = "\\.jpeg",
                   recursive = TRUE,
                   full.names = TRUE)
logs <- "results_simulations/logs.txt"
file.remove("results_simulations/results.zip")
zip("results_simulations/results.zip", files = c(processed_results_files, weighting_plots_files, average_outcome_files, parameters_files, pngs, pdfs, logs, jpegs))
