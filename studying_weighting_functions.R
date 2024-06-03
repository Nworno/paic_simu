######################################
##### Studying weighting distributions
######################################
library(data.table)
source("estimators.R")
options(mc.cores = 8)
N_pop <- 10^6
N_BOOT_ITER <- 100
# Parameters
####################

df_population_parameters <- list(
  prop_X1 = 0.5, # Variable binaire, prevalence dans la population
  bT_X1 = 0,   # Effet de la variable binaire sur la probabilité d'être dans l'essai AC
  bT_X2 = c(0, 1.5, 5),   # Effet de la variable continue X2...
  bT_X3 = 0,     # Idem, mais inutile pour le moment
  bT_X4 = 0,     # Idem, mais inutile pour le moment
  bY_X1 = 1.5,   # Effet de X1 sur l'outcome
  bY_X2 = 0.5,   # Effet de X2 sur l'outcome
  bY_X3 = 0,     # Effet de X3 sur l'outcome (inutile pour le moment)
  bY_X4 = 0.5,     # Effet de X4 sur l'outcome
  bY_A_X1 = c(0), # Interaction A et X1 dans le modèle outcome
  bY_A_X2 = c(0), # Interaction A et X2 dans le modèle outcome
  bY_A_X3 = c(0), # Interaction A et X3 dans le modèle outcome
  bY_A_X4 = c(0), # Interaction A et X4 dans le modèle outcome
  binary_marker = c(bquote(rbinom(N_pop, 1, 0.5))), # Utilisé pour la variable bimodale
  f_X1 = c(bquote(rbinom(N_pop, 1, 0.5))), # distribution de X1 (revoir car il faudrait utiliser prop_X1)
  f_X2 = c(bquote(rnorm(N_pop, 0.5, 1)),
           bquote(binary_marker * rnorm(N_pop, -2, 1) + (1 - binary_marker) * rnorm(N_pop, 2, 1))), # distribution de X2
  f_X3 = c(bquote(0)), # inutile pour le moment
  # f_X4 = c(bquote(binary_marker * rnorm(N_pop, -1.5, 1) + (1 - binary_marker) * rnorm(N_pop, 1.5, 1))),
  f_X4 = c(bquote(0)), # inutile pour le moment
  bY_A = 1.5,  # Effet de A par rapport à C
  bY_B = 1.5,  # Effet de B par rapport à C
  bY_C = 0,    # Pas d'effet de C sur l'outcome
  AC_trial_model = c(bquote(X1 * bT_X1 + X2 * bT_X2 + X3 * bT_X3 + X4 * bT_X4)), # Modèle d'attribution de l'essai AC
  BC_trial_model = c(bquote(0)),  # Modèle d'attribution de l'essai BC
  outcome_generation_formula =  c(bquote(
    bY_X1*X1 + bY_X2 * X2 + bY_X3 * X3 + bY_X4 * X4 + (bY_A + bY_A_X1*X1 + bY_A_X2*X2 + bY_A_X3*X3 + bY_A_X4*X4) * A +  bY_B*B + bY_C*C
  ))
) |> 
  expand.grid() |>
  as.data.table()
df_population_parameters[, population_parameters_num := 1:.N]
df_population_parameters <- df_population_parameters[population_parameters_num == 6,]

##############################################
########### Creating an overarching population
##############################################

## Génère une data.frame de 10^6 ou 7 lignes
creating_population <- function(list_simulation_parameters) {
  attach(list_simulation_parameters)
  binary_marker <- rbinom(N_pop, 1, 0.5)
  pop_init <- data.table(
    id = 1:N_pop,
    X1 = eval(f_X1),
    X2 = eval(f_X2),
    X3 = rlnorm(N_pop, 0.5, 0.5),
    X4 = binary_marker * rnorm(N_pop, -1.5, 1) + (1 - binary_marker) * rnorm(N_pop, 1.5, 1) # tentative d'une variable bimodale (mais pas utilisé finalement, coef à zéro)
  ) |>
    setkey("id")
  
  
  trial_assignement_prob <- function(trial_assignment_model, df) {
    predicted <- with(df, eval(trial_assignment_model))
    return(plogis(predicted))
  }
  
  predict_outcome <- function(outcome_model, df) {
    with(df, eval(outcome_model))
  }
  
  pop_init[, prob_w_trial_AC := trial_assignement_prob(AC_trial_model, df = pop_init)]
  pop_init[, prob_w_trial_BC := trial_assignement_prob(BC_trial_model, df = pop_init)]
  # AC_pop_init <- pop_init[sample(id, size = 10^6, prob = prob_w_trial_AC, replace = TRUE), ]
  # Drawing a very large population with BC characteristics to compute theoretical outcomes
  BC_pop_init <- pop_init[sample(id, size = 10^6, prob = prob_w_trial_BC, replace = TRUE), ] 
  
  covariate_names <- c("X1", "X2", "X3", "X4")
  df_outcomes <- sapply(list(A = BC_pop_init[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                             B = BC_pop_init[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                             C = BC_pop_init[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
                        predict_outcome,
                        outcome_model = outcome_generation_formula,
                        simplify = FALSE) |>
    c("id" = list(1:nrow(pop_init))) |>
    as.data.table() |>
    melt(id.vars = c("id"), variable.name = "ttt", value.name = "Y_theo")
  df_outcomes[, Y_obs := Y_theo + rnorm(n = length(Y_theo), mean = 0, sd = 1)]
  
  average_BC_pop_init <- BC_pop_init[, lapply(.SD, mean), .SDcols = covariate_names]
  
  average_conditional_outcome <- sapply(list("A" = average_BC_pop_init[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                                             "B" = average_BC_pop_init[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                                             "C" = average_BC_pop_init[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
                                        predict_outcome,
                                        outcome_model = outcome_generation_formula,
                                        simplify = FALSE) |>
    as.data.table() |>
    melt(measure.vars = c("A", "B", "C"), variable.name = "ttt", value.name = "outcome")
  marginal_outcome <- df_outcomes[, .(outcome = mean(Y_obs)), by = c("ttt")]
  average_outcome_df <- rbindlist(
    list("conditional" = average_conditional_outcome,
         "marginal" = marginal_outcome),
    use.names = TRUE,
    idcol = "outcome_type") |>
    dcast(outcome_type ~ ttt, value.var = "outcome")
  average_outcome_df[, AB := A - B]
  if (average_outcome_df[outcome_type == "conditional", AB] == 0) average_outcome_df[, AB := 0]
  
  return(list(
    "pop_init" = pop_init,
    "df_outcomes" = df_outcomes,
    "average_outcome_df" = average_outcome_df
  ))
}

# méthodes de comparaison indirecte
indirect_comparisons <- function(pop_init,
                                 df_outcomes,
                                 N_RCT,
                                 covariate_names) { # info score de propension
  
  #############################
  ############## Drawing trials
  #############################
  
  # selected_individuals_AC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_AC)][
  #   , ttt := rep_len(c("A", "C"), length.out = .N)]
  # selected_outcomes_AC <- df_outcomes[selected_individuals_AC, on = c("id", "ttt")][
  #   , c("id", "ttt", "Y_obs")]
  # trial_AC <- pop_init[selected_outcomes_AC, on = "id"][, ttt := factor(ttt, levels = c("C", "A"))]
  selected_individuals_AC_A <- pop_init[trial.f == "AC" & trt.f == "A"][sample(1:.N, N_RCT/2, replace = FALSE)]
  selected_individuals_AC_C <- pop_init[trial.f == "AC" & trt.f == "A"][sample(1:.N, N_RCT/2, replace = FALSE)]
  trial_AC <- rbindlist(list(selected_individuals_AC_A, selected_individuals_AC_C))[, .(id, X1, X2, ttt = trt.f, Y_obs = Y)]
  
  # selected_individuals_BC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_BC)][
  #   , ttt := rep_len(c("B", "C"), length.out = .N)]
  # selected_outcomes_BC <- df_outcomes[selected_individuals_BC, on = c("id", "ttt")][
  #   , c("id", "ttt", "Y_obs")]
  # trial_BC <- pop_init[selected_outcomes_BC, on = "id"][, ttt := factor(ttt, levels = c("C", "B"))]
  selected_individuals_BC_B <- pop_init[trial.f == "BC" & trt.f == "B"][sample(1:.N, N_RCT/2, replace = FALSE)]
  selected_individuals_BC_C <- pop_init[trial.f == "BC" & trt.f == "C"][sample(1:.N, N_RCT/2, replace = FALSE)]
  trial_BC <- rbindlist(list(selected_individuals_BC_B, selected_individuals_BC_C))[, .(id, X1, X2, ttt = trt.f, Y_obs = Y)]
  
  is_anchored <- FALSE
  
  ps_w <- propensity_score(trial_AC, 
                           trial_BC, 
                           covariate_names, 
                           anchored = is_anchored,
                           weight_estimation_method = "max_likelihood",
                           outcome_family = gaussian(link = "identity"),
                           studying_populations = TRUE)
  maic_w <- propensity_score(trial_AC, 
                           trial_BC, 
                           covariate_names, 
                           anchored = is_anchored,
                           weight_estimation_method = "moments",
                           outcome_family = gaussian(link = "identity"),
                           studying_populations = TRUE)
  if (!is_anchored) trial_AC = trial_AC[ttt == "A", ]
  if (!is_anchored) trial_AC = trial_AC[ttt == "A", ]
  trial_AC$maic_w <- maic_w
  trial_AC$ps_w <- ps_w
  # print(trial_AC)
  return(list("trial_AC" = trial_AC,
              "trial_BC" = trial_BC))
}

# draw_trials <- function(pop_init, N_RCT) {
#   selected_individuals_AC <- pop_init[sample(id, N_RCT, replace = TRUE, prob = prob_w_trial_AC)][
#     , ttt := rep_len(c("A", "C"), length.out = .N)]
# 
#   selected_individuals_BC <- pop_init[sample(id, N_RCT, replace = TRUE, prob = prob_w_trial_BC)][
#     , ttt := rep_len(c("B", "C"), length.out = .N)]
#   return(list(trial_AC = selected_individuals_AC, 
#               trial_BC = selected_individuals_BC))
# }

# Used to specify variables to use for "trial exposure" models, and unanchored STC
list_covariate_names <- c(
  combn(c("X1", "X2"), m = 1, simplify = FALSE),
  combn(c("X1", "X2"), m = 2, simplify = FALSE)
) |> as.vector()


df_estimators_parameters <- data.table(
  "covariate_names" = list_covariate_names,
  N_RCT = c(60, 200, 1000)
)
df_estimators_parameters[, estimator_num := 1:.N]
# Selecting only estimator with X2 only
# df_estimators_parameters = df_estimators_parameters[estimator_num == 2,]

n_iter = 400
time_start <- Sys.time()
time_start_string <- format(time_start, "%Y%m%d_%H%M%S")
print(time_start_string)
# experiment_results_directory <- file.path("studying_weighting_fn", time_start_string)
experiment_results_directory <- file.path("david_10052024/studying_weighting_fn")
if (!dir.exists(experiment_results_directory)) dir.create(experiment_results_directory)
# saveRDS(df_population_parameters, file.path(experiment_results_directory, "df_population_parameters.RDS"))
saveRDS(df_estimators_parameters, file.path(experiment_results_directory, "df_estimators_parameters.RDS"))
n_scenarios <- c("6", "10")
# for (row_population in 1:nrow(df_population_parameters)) {
for (n_scenario in n_scenarios) {
  # dir_sub_experiment <- file.path(experiment_results_directory, row_population)
  dir_sub_experiment <- file.path(experiment_results_directory, n_scenario)
  dir.create(dir_sub_experiment)
  # list_simulation_parameters <- df_population_parameters[row_population, ] |> unlist()
  # population <- creating_population(list_simulation_parameters) # pop initiale
  # df_outcomes <- population$df_outcomes
  # pop_init <- population$pop_init
  load(file = paste0("david_10052024/simulated_data_scenario", n_scenario, ".RDS"))
  pop_init <- data.table::rbindlist(list(dfAC, dfBC))
  
  pop_init[, .(id, X1, X2, )]
  # df_outcomes <- pop_init[, .(id, Y_obs = Y)]
  
  for (row_estimators in 1:nrow(df_estimators_parameters)) {
    list_estimators_parameters <- df_estimators_parameters[row_estimators, ] |> unlist(recursive = FALSE)
    results_simulations <- parallel::mclapply(1:n_iter, \(i) {
      time_start_iteration <- Sys.time()
      list_trials <- indirect_comparisons(pop_init,
                                          df_outcomes = df_outcomes,
                                          list_estimators_parameters[["N_RCT"]],
                                          list_estimators_parameters[["covariate_names"]])
      time_eluded <- Sys.time() - time_start_iteration
      cat("Experiment ", row_population, ".", row_estimators, ", Iteration ", i, ", length: ", time_eluded, " seconds\n", sep = "")
      return(list_trials)
    })
    saveRDS(results_simulations, file = file.path(dir_sub_experiment,
                                          paste0("experiment_", row_estimators, ".RDS")))
  }
}
cat("Simulation length: ")
print(Sys.time() - time_start)



