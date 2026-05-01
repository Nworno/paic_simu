options(mc.cores = 8)

source("estimators.R")
source("data_generation.R")

####################
# Parameters
####################
N_pop <- 2*10^6 # Size of the superpopulation, from which trials are drawn for each iteration
N_BOOT_ITER <- 200 # Number of bootstrap iterations for variance estimation
n_iter <- 1000 # Number of Monte-Carlo iterations
retrieve_ps_weights <- TRUE # Whether to save each trial data, in order to be able to plot propensity score for each weighting estimator


###############
### SIMULATIONS
###############
time_start <- Sys.time()
time_start_string <- format(time_start, "%Y%m%d_%H%M%S")
print(time_start_string)
experiment_results_directory <- file.path("results_simulations", time_start_string)
if (!dir.exists(experiment_results_directory)) dir.create(experiment_results_directory, recursive = TRUE)
saveRDS(df_population_parameters, file.path(experiment_results_directory, "df_population_parameters.RDS"))
saveRDS(df_estimators_parameters, file.path(experiment_results_directory, "df_estimators_parameters.RDS"))
data.table::setDTthreads(1)  # disable data.table threading to avoid conflicts with mclapply forking
for (row_population in 1:nrow(df_population_parameters)) {
  list_simulation_parameters <- df_population_parameters[row_population, ] |> unlist()
  dir_sub_experiment <- file.path(experiment_results_directory, list_simulation_parameters[["population_parameters_num"]])
  dir.create(dir_sub_experiment)

  print("creating population")
  print(Sys.time())
  population <- creating_population(list_simulation_parameters, N_pop) # pop initial
  pop_init <- population$pop_init
  print("population created")
  print(Sys.time())

  saveRDS(population$average_outcome_df, file.path(dir_sub_experiment, "average_outcome_df.RDS"))

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
                                                         list_simulation_parameters[["outcome_distribution"]],
                                                         retrieve_ps_weights = retrieve_ps_weights)
      time_eluded <- Sys.time() - time_start_iteration
      cat("Experiment ", row_population, ".", row_estimators, ", Iteration ", i, ", length: \n", sep = "")
      print(time_eluded)

      return(result_indirect_comparison)
    })
    if (retrieve_ps_weights) {
      saveRDS(sapply(results_simulations, \(x) x$struct_ps_df, USE.NAMES = TRUE, simplify = FALSE),
              file = file.path(dir_sub_experiment,
                               paste0("experiment_dfs_", row_estimators, ".RDS")))
    }
    saveRDS(sapply(results_simulations, \(x) x$list_dfs_errors, USE.NAMES = TRUE, simplify = FALSE),
            file = file.path(dir_sub_experiment,
                             paste0("experiment_dfs_errors_", row_estimators, ".RDS")))
    saveRDS(sapply(results_simulations, \(x) x$list_warnings, USE.NAMES = TRUE, simplify = FALSE),
            file = file.path(dir_sub_experiment,
                             paste0("experiment_warnings_", row_estimators, ".RDS")))
    saveRDS(sapply(results_simulations, \(x) x$rectangle_results, USE.NAMES = TRUE, simplify = FALSE),
            file = file.path(dir_sub_experiment,
                             paste0("experiment_results_", row_estimators, ".RDS")))
  }
}
cat("Simulation length: ")
print(Sys.time() - time_start)
print(time_start_string)
write(paste0("DATE_EXPERIMENT <- '", time_start_string, "'"), file = "env_variables.R", append = FALSE)
print("processing results")
source("processing_results.R")
print("graph distribution covariates")
source("graph_distribution_covariates.R")
source("graphs_publication.R")
if (retrieve_ps_weights) source("graph_propensity_scores.R")
source("zipping_results.R") # Useful if executed on a remote server, to pull results easily
beepr::beep(3)
Sys.sleep(4)
