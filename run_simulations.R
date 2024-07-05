options(mc.cores = 1)

source("estimators.R")
source("data_generation.R")

####################
# Parameters
####################
N_pop <- 10^6
N_BOOT_ITER <- 200
n_iter <- 500

###############
### SIMULATIONS
###############
time_start <- Sys.time()
time_start_string <- format(time_start, "%Y%m%d_%H%M%S")
print(time_start_string)
experiment_results_directory <- file.path("results_simulations", time_start_string)
if (!dir.exists(experiment_results_directory)) dir.create(experiment_results_directory)
saveRDS(df_population_parameters, file.path(experiment_results_directory, "df_population_parameters.RDS"))
saveRDS(df_estimators_parameters, file.path(experiment_results_directory, "df_estimators_parameters.RDS"))
for (row_population in 1:nrow(df_population_parameters)) {
  dir_sub_experiment <- file.path(experiment_results_directory, row_population)
  dir.create(dir_sub_experiment)
  list_simulation_parameters <- df_population_parameters[row_population, ] |> unlist()

  population <- creating_population(list_simulation_parameters) # pop initial
  pop_init <- population$pop_init # données simulées

  saveRDS(population$average_outcome_df, file.path(dir_sub_experiment, "average_outcome_df.RDS"))
  # Commented because huge file, so would take could much space if saved for every try
  # saveRDS(pop_init, file.path(dir_sub_experiment, "pop_init.RDS"))

  for (row_estimators in 1:nrow(df_estimators_parameters)) {
    list_estimators_parameters <- df_estimators_parameters[row_estimators, ] |> unlist(recursive = FALSE)
    results_simulations <- parallel::mclapply(1:n_iter, \(i) {
      time_start_iteration <- Sys.time()
      result_indirect_comparison <- indirect_comparisons(pop_init,
                                                         struct_results,
                                                         N_BOOT_ITER,
                                                         list_simulation_parameters[["N_RCT"]],
                                                         list_estimators_parameters[["outcome_regression_model"]],
                                                         list_estimators_parameters[["covariate_names"]],
                                                         list_estimators_parameters[["assignment_model"]],
                                                         list_simulation_parameters[["outcome_distribution"]])
      time_eluded <- Sys.time() - time_start_iteration
      cat("Experiment ", row_population, ".", row_estimators, ", Iteration ", i, ", length: ", time_eluded, " seconds\n", sep = "")
      return(result_indirect_comparison)
    })
    saveRDS(results_simulations, file = file.path(dir_sub_experiment,
                                                  paste0("experiment_", row_estimators, ".RDS")))
  }
}
cat("Simulation length: ")
print(Sys.time() - time_start)
print(time_start_string)
write(paste0("DATE_EXPERIMENT <- '", time_start_string, "'"), file = "env_variables.R", append = FALSE)
