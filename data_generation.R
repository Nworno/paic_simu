library(data.table)

options(mc.cores = 1)
source("estimators.R")
N_pop <- 10^6
N_BOOT_ITER <- 2
# Parameters
####################

df_population_parameters <- list(
  prop_X1 = 0.5, # Variable binaire, prevalence dans la population
  bT_X1 = 0.5,   # Effet de la variable binaire sur la probabilité d'être dans l'essai AC
  bT_X2 = 0.2,   # Effet de la variable continue X2...
  bT_X3 = 0,     # Idem, mais inutile pour le moment
  bT_X4 = 0,     # Idem, mais inutile pour le moment
  bY_X1 = 1.5,   # Effet de X1 sur l'outcome
  bY_X2 = 0.5,   # Effet de X2 sur l'outcome
  bY_X3 = 0,     # Effet de X3 sur l'outcome (inutile pour le moment)
  bY_X4 = 0,     # Effet de X4 sur l'outcome (inutile pour le moment)
  bY_A_X1 = c(0, 1.2), # Interaction A et X1 dans le modèle outcome
  bY_A_X2 = c(0, 0.2), # Interaction A et X2 dans le modèle outcome
  f_X1 = c(bquote(rbinom(N_pop, 1, 0.5))), # distribution de X1 (revoir car il faudrait utiliser prop_X1)
  f_X2 = c(bquote(rnorm(N_pop, 0.5, 1)), bquote(rlnorm(N_pop, 0.5, 0.5))), # distribution de X2
  f_X3 = c(bquote(0)), # inutile pour le moment
  f_X4 = c(bquote(0)), # inutile pour le moment
  bY_A = 1.5,  # Effet de A par rapport à C
  bY_B = 1.5,  # Effet de B par rapport à C
  bY_C = 0,    # Pas d'effet de C sur l'outcome
  AC_trial_model = c(bquote(X1 * bT_X1 + X2 * bT_X2)), # Modèle d'attribution de l'essai AC
  BC_trial_model = c(bquote(0)),  # Modèle d'attribution de l'essai BC
  outcome_generation_formula = c(
  bquote(bY_X1*X1 + bY_X2 * X2 + bY_X3 * X3 + bY_X4 * X4 + (bY_A + bY_A_X1*X1 + bY_A_X2*X2) * A +  bY_B*B + bY_C*C))) |> 
  expand.grid() |>
  as.data.table()
df_population_parameters[, population_parameters_num := 1:.N]


######### Models
##################
## Note David: we could have
## one binary variable X1
## one continuous variable X2
## one parameter to switch the binary variable between prognostic only (bYA_X1 = 0) or effect modifier (bYA_X1 != 0)
## one parameter to switch the continuous variable between prognostic only (bYA_X2 = 0) or effect modifier (bYA_X2 != 0)
## one parameter to switch the continuous variable distribution: normal (symmetrical) or lognormal (asymmetrical)
## Overall, 8 scenarios here
##
## Then, run all estimators. For estimators taking into accounts covariates, run three estimations:
## - only X1
## - only X2
## - both X1 and X2
## For estimators taking into accounts moments, run with (to be discussed):
## - first moment only
## - first and second moments

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

  covariate_names <- c("X1", "X2", "X3", "X4")
  df_outcomes <- sapply(list(A = pop_init[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                             B = pop_init[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                             C = pop_init[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
                        predict_outcome,
                        outcome_model = outcome_generation_formula,
                        simplify = FALSE) |>
    c("id" = list(1:nrow(pop_init))) |>
    as.data.table() |>
    melt(id.vars = c("id"), variable.name = "ttt", value.name = "Y_theo")
  df_outcomes[, Y_obs := Y_theo + rnorm(n = length(Y_theo), mean = 0, sd = 1)]

  average_pop_init <- pop_init[, lapply(.SD, mean), .SDcols = covariate_names]

  average_conditional_outcome <- sapply(list("A" = average_pop_init[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                                             "B" = average_pop_init[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                                             "C" = average_pop_init[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
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

  return(list(
    "pop_init" = pop_init,
    "df_outcomes" = df_outcomes,
    "average_outcome_df" = average_outcome_df
  ))
}

# méthodes de comparaison indirecte
indirect_comparisons <- function(pop_init,
                                 df_outcomes,
                                 struct_results,
                                 N_BOOT_ITER,
                                 N_RCT,
                                 outcome_regression_model, # info modèle de régresion
                                 covariate_names) { # info score de propension

  #############################
  ############## Drawing trials
  #############################

  selected_individuals_AC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_AC)][
    , ttt := rep_len(c("A", "C"), length.out = .N)]
  selected_outcomes_AC <- df_outcomes[selected_individuals_AC, on = c("id", "ttt")][
    , c("id", "ttt", "Y_obs")]
  trial_AC <- pop_init[selected_outcomes_AC, on = "id"][, ttt := factor(ttt, levels = c("C", "A"))]

  selected_individuals_BC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_BC)][
      , ttt := rep_len(c("B", "C"), length.out = .N)]
  selected_outcomes_BC <- df_outcomes[selected_individuals_BC, on = c("id", "ttt")][
    , c("id", "ttt", "Y_obs")]
  trial_BC <- pop_init[selected_outcomes_BC, on = "id"][, ttt := factor(ttt, levels = c("C", "B"))]

  average_trial_BC_covariates <- trial_BC[, lapply(.SD, mean), .SDcols = covariate_names]
  centered_trial_AC <- trial_AC[, mapply(`-`, .SD, average_trial_BC_covariates), .SDcols = covariate_names] |>
    cbind(trial_AC[, .(ttt, Y_obs)])
  centered_trial_BC <- trial_BC[, mapply(`-`, .SD, average_trial_BC_covariates), .SDcols = covariate_names] |>
    cbind(trial_BC[, .(ttt, Y_obs)])

  ###############################
  ########## Unadjusted estimator
  ###############################
  struct_results$unadjusted$anchored$unadjusted <- run_unadjusted_estimator(trial_AC, trial_BC, anchored = TRUE)
  struct_results$unadjusted$unanchored$unadjusted <- run_unadjusted_estimator(trial_AC, trial_BC, anchored = FALSE)

  ##############################################################
  ########## REGRESSION BASED OUTCOME MODEL (both treatment IPD)
  ##############################################################

  struct_results$regression$anchored$glm <- run_anchored_conditional_estimation(centered_trial_AC,
                                                                                centered_trial_BC,
                                                                                outcome_regression_model,
                                                                                gaussian)
  struct_results$regression$unanchored$glm <- run_unanchored_conditional_estimation(centered_trial_AC,
                                                                                    centered_trial_BC,
                                                                                    outcome_regression_model,
                                                                                    gaussian)

  #################################################
  ########### PROPENSITY SCORE (both treatment IPD)
  #################################################

  struct_results$iptw$anchored$ml <- run_propensity_score(trial_AC, trial_BC, anchored = TRUE, covariate_names)
  struct_results$iptw$unanchored$ml <- run_propensity_score(trial_AC, trial_BC, anchored = FALSE, covariate_names)

  #########
  ### MAIC
  #########
  struct_results$iptw$anchored$maic <- run_maic(trial_AC, trial_BC, covariate_names, anchored = TRUE, outcome_family = gaussian)
  struct_results$iptw$unanchored$maic <- run_maic(trial_AC, trial_BC, covariate_names, anchored = FALSE, outcome_family = gaussian)

  ##########
  ###### STC
  ##########
  struct_results$regression$anchored$stc <- run_stc(centered_trial_AC,
                                                    trial_BC,
                                                    anchored = TRUE,
                                                    outcome_regression_model,
                                                    covariate_names = NULL,
                                                    outcome_family = gaussian)
  struct_results$regression$unanchored$stc <- run_stc(centered_trial_AC,
                                                      trial_BC,
                                                      anchored = FALSE,
                                                      outcome_regression_model = NULL,
                                                      covariate_names,
                                                      outcome_family = gaussian)


  ##############################
  ############ COMPILING RESULTS
  ##############################

  rectangle_results <- struct_results |> tibble::enframe() |>
    tidyr::unnest_longer(value, indices_to = "anchored") |>
    tidyr::unnest_longer(value, indices_to = "model") |>
    tidyr::unnest_longer(value, indices_to = "indicator") |>
    dplyr::rename(adjustment = name) |>
    tidyr::pivot_wider(names_from = indicator, values_from = value)

  return(rectangle_results)
}




struct_results <- list(
  unadjusted = list(
    anchored = list(agd = NULL),
    unanchored = list(agd = NULL)
  ),
  regression = list(
    anchored = list(glm = NULL, stc = NULL),
    unanchored = list(glm = NULL, stc = NULL)
  ),
  iptw = list(
    anchored = list(ml = NULL, maic = NULL),
    unanchored = list(ml = NULL, maic = NULL)
  )
)

###############
#### ESTIMATORS
###############

# Used to specify variables to use for "trial exposure" models, and unanchored STC
list_covariate_names <- c(
  combn(c("X1", "X2"), m = 1, simplify = FALSE),
  combn(c("X1", "X2"), m = 2, simplify = FALSE)
) |> as.vector()

list_outcome_regression_models <- paste0("Y_obs ~ ", c(
  "X1*ttt",
  "X2*ttt",
  "X1*ttt + X2*ttt"
))

df_estimators_parameters <- data.table(
  "covariate_names" = list_covariate_names,
  "outcome_regression_model" = list_outcome_regression_models,
  N_RCT = 200
)
df_estimators_parameters[, estimator_num := 1:.N]


# Intermediate data frame used for creating all combinations of
# combination_experiments <- expand.grid(
#   "estimator_num" = df_estimators_parameters$estimator_num,
#   "population_parameters_num" = df_population_parameters$population_parameters_num
# ) |> as.data.table()
#
# df_experiments <- combination_experiments[df_population_parameters, , on = "population_parameters_num"][
#   df_estimators_parameters, , on = "estimator_num"
# ]



###############
### SIMULATIONS
###############
n_iter = 3
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
  population <- creating_population(list_simulation_parameters) # pop initiale
  pop_init <- population$pop_init # données simulées
  df_outcomes <- population$df_outcomes # outcome théorique par patient
  average_outcome_df <- population$average_outcome_df # moyenne de ces outcomes (estimation empirique de l'effet marginal et conditionnel)

  saveRDS(average_outcome_df, file.path(dir_sub_experiment, "average_outcome_df.RDS"))
  for (row_estimators in 1:nrow(df_estimators_parameters)) {
    list_estimators_parameters <- df_estimators_parameters[row_estimators, ] |> unlist(recursive = FALSE)
    results_simulations <- parallel::mclapply(1:n_iter, \(i) {
      time_start_iteration <- Sys.time()
      result_indirect_comparison <- indirect_comparisons(pop_init,
                                                         df_outcomes,
                                                         struct_results,
                                                         N_BOOT_ITER,
                                                         list_estimators_parameters[["N_RCT"]],
                                                         list_estimators_parameters[["outcome_regression_model"]],
                                                         list_estimators_parameters[["covariate_names"]])
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

