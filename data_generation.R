library(data.table)

##############################################
########### Data Generation Parameters
##############################################

df_default_parameters <- list(
  N_RCT  = c(500),
  bT = 0.5,
  bT_X1 = 1,   # Effet de la variable binaire sur la probabilité d'être dans l'essai BC
  bT2_X1 = -0.5,   # Effet de la variable binaire sur la probabilité d'être dans l'essai BC
  bT_X2 = 1,   # Effet de la variable continue X2...
  bT2_X2 = -0.5,   # Effet de la variable continue X2...
  fbT = 1,
  bT_X3 = 0,
  bT_X4 = 0,
  bT2_X3 = 0,
  bT2_X4 = 0,
  bY_X1 = 1,   # Effet de X1 sur l'outcome
  bY_X2 = 2,   # Effet de X2 sur l'outcome
  bY_X3 = 0,     # Effet de X3 sur l'outcome
  bY_X4 = 0,     # Effet de X4 sur l'outcome
  bY_A_X1 = c(2), # Interaction A et X1 dans le modèle outcome
  bY_A_X2 = c(0), # Interaction A et X2 dans le modèle outcome
  bY_A_X3 = c(0), # Interaction A et X3 dans le modèle outcome
  bY_A_X4 = c(0), # Interaction A et X4 dans le modèle outcome
  bY_B_X1 = c(0), # Interaction A et X1 dans le modèle outcome
  bY_B_X2 = c(0), # Interaction A et X2 dans le modèle outcome
  bY_B_X3 = c(0), # Interaction A et X3 dans le modèle outcome
  bY_B_X4 = c(0), # Interaction A et X4 dans le modèle outcome
  # binary_marker = c(bquote(rbinom(N_pop, 1, 0.5))), # Utilisé pour la variable bimodale : introduit corrélation entre les variables, car le marker est le même pour tous les individus
  f_X1 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1.5))))),
  f_X2 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1.5))))),
  f_X3 = c(bquote(0)),
  f_X4 = c(bquote(0)),
  bY_A = 1,  # Effet de A par rapport à C
  bY_B = 1,  # Effet de B par rapport à C
  bY_C = 0,    # Pas d'effet de C sur l'outcome
  BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2 + bT_X3 * X3 + bT2_X3 * X3^2 + bT_X4 * X4 + bT2_X4 * X4^2)), # Modèle d'attribution de l'essai BC
  outcome_distribution = c("normal"),
  outcome_generation_formula =  c(bquote(
    bY_X1*X1 + bY_X2 * X2 + bY_X3 * X3 + bY_X4 * X4 + (bY_A + bY_A_X1*X1 + bY_A_X2*X2 + bY_A_X3*X3 + bY_A_X4*X4) * A +  (bY_B + bY_B_X1*X1 + bY_B_X2*X2 + bY_B_X3*X3 + bY_B_X4*X4) *B + bY_C*C
  ))
)


list_changing_parameters <- list(
  "1" = list(
    f_X1 = c(bquote(rnorm(N_pop, 0, 1))),
    f_X2 = c(bquote(rnorm(N_pop, 0, 1))),
    bT = 0.5,
    bT_X1 = 1,
    bT2_X1 = -0.5,
    bT_X2 = 1,
    bT2_X2 = -0.5,
    fbT = 1,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "2" = list(
    f_X1 = c(bquote(rnorm(N_pop, 0, 1))),
    f_X2 = c(bquote(rnorm(N_pop, 0, 1))),
    bT = 0.5,
    bT_X1 = 1,
    bT2_X1 = -0.5,
    bT_X2 = 1,
    bT2_X2 = -0.5,
    fbT = -1,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
    ),
  "3" = list(
    f_X1 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    f_X2 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    bT = 1,
    bT_X1 = 1,
    bT2_X1 = -1,
    bT_X2 = 1,
    bT2_X2 = -1,
    fbT = 0.5,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "4" = list(
    f_X1 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    f_X2 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    bT = 2,
    bT_X1 = 1,
    bT2_X1 = -1,
    bT_X2 = 1,
    bT2_X2 = -1,
    fbT = -0.5,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "5" = list(
    f_X1 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    f_X2 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    bT = 1,
    bT_X1 = 2,
    bT2_X1 = 0,
    bT_X2 = 2,
    bT2_X2 = 0,
    fbT = 1,
    bY_X1 = 0,
    bY_X2 = 0,
    bY_A_X1 = 2,
    bY_A_X2 = 0,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "6" = list(
    f_X1 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    f_X2 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    bT = 1,
    bT_X1 = 2,
    bT2_X1 = 0,
    bT_X2 = 2,
    bT2_X2 = 0,
    fbT = -1,
    bY_X1 = 0,
    bY_X2 = 0,
    bY_A_X1 = 2,
    bY_A_X2 = 0,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "7" = list(
    f_X1 = c(bquote(sample(c(rnorm(N_pop/2, 0, 0.5), rnorm(N_pop/2, 3, 0.5))))),
    f_X2 = c(bquote(sample(c(rnorm(N_pop/2, 0, 0.5), rnorm(N_pop/2, 3, 0.5))))),
    bT = 0.5,
    bT_X1 = 2,
    bT2_X1 = -2,
    bT_X2 = 2,
    bT2_X2 = -2,
    fbT = 1,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^3 + bT_X2 * X2 + bT2_X2 * X2^3))
  ),
  "8" = list(
    f_X1 = c(bquote(sample(c(rnorm(N_pop/2, 0, 0.5), rnorm(N_pop/2, 3, 0.5))))),
    f_X2 = c(bquote(sample(c(rnorm(N_pop/2, 0, 0.5), rnorm(N_pop/2, 3, 0.5))))),
    bT = 0.5,
    bT_X1 = 2,
    bT2_X1 = -2,
    bT_X2 = 2,
    bT2_X2 = -2,
    fbT = -1,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^3 + bT_X2 * X2 + bT2_X2 * X2^3))
  )
)

list_parameters <- lapply(list_changing_parameters, \(x) {
  list_parameter <- c(x, df_default_parameters[!names(df_default_parameters) %in% names(x)], recursive = TRUE)
  list_parameter <- list_parameter[names(df_default_parameters)]
})

df_population_parameters <- list_parameters |> tibble::as_tibble() |> t()
colnames(df_population_parameters) <- names(df_default_parameters)
df_population_parameters <- data.table::as.data.table(df_population_parameters)
df_population_parameters[, population_parameters_num := 1:.N]


##############################################
########### Creating an overarching population
##############################################

## Génère une data.frame de 10^6 ou 7 lignes
creating_population <- function(list_simulation_parameters) {
  attach(list_simulation_parameters)

  # Modifying the coefficient values using fbT
  bT <- bT*fbT
  bT_X1 <- bT_X1*fbT
  bT2_X1 <- bT2_X1*fbT
  bT_X2 <- bT_X2*fbT
  bT2_X2 <- bT2_X2*fbT
  # bT_X3 <- bT_X3*fbT
  # bT2_X3 <- bT2_X3*fbT
  # bT_X4 <- bT_X4*fbT
  # bT2_X4 <- bT2_X4*fbT

  pop_init <- data.table(
    id = 1:N_pop,
    X1 = eval(f_X1),
    X2 = eval(f_X2),
    X3 = eval(f_X3),
    X4 = eval(f_X4)
  ) |>
    setkey("id")


  trial_assignement_prob <- function(trial_assignment_model, df) {
    predicted <- with(df, eval(trial_assignment_model))
    return(plogis(predicted))
  }

  predict_outcome <- function(outcome_model, df) {
    with(df, eval(outcome_model))
  }

  covariate_names <- c("X1", "X2", "X3", "X4")

  df_outcomes_pop_init <- sapply(list(A = pop_init[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                                      B = pop_init[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                                      C = pop_init[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
                                 predict_outcome,
                                 outcome_model = outcome_generation_formula,
                                 simplify = FALSE) |>
    c("id" =  list(pop_init$id),
      "trial" = list(pop_init$trial)) |>
    as.data.table() |>
    melt(id.vars = c("id"), variable.name = "ttt", value.name = "Y_theo")

  if (outcome_distribution == "normal") {
    df_outcomes_pop_init[, Y_obs := Y_theo + rnorm(n = .N, mean = 0, sd = 1)]
  } else if (outcome_distribution == "binomial") {
    df_outcomes_pop_init[, Y_obs := rbinom(n = .N, size = 1, prob = plogis(Y_theo))] # equivalent to ifelse(rnorm(.N, 0, pi/sqrt(3)) > 0, 1, 0) bc of the variance of the logistic function, thus normal and binomial outcome distributions are not interchangeable using plogis and 1/(1 + exp(-x))
    moy_outcomes <- tapply(df_outcomes_pop_init$Y_obs, df_outcomes_pop_init$ttt, mean, simplify = FALSE)
    if (any(moy_outcomes < 0.02 | moy_outcomes > 0.98)) { # arbitrary thresholds, to avoid downstreams problem with model fitting
      stop("Too extreme outcomes")
    }
  } else {
    stop("Unknown outcome distribution")
  }


  pop_init <- df_outcomes_pop_init[, Y_theo:= NULL][pop_init, on = "id"] |> data.table::dcast(formula = ... ~ ttt, value.var = "Y_obs")

  pop_init[, prob_BC := trial_assignement_prob(BC_trial_model, df = pop_init)]
  pop_init[, trial := rbinom(.N, 1, prob_BC) |>
             factor(levels = c(0, 1), labels = c("AC", "BC"))]

  print(table(pop_init$trial)/N_RCT)
  if (any(table(pop_init$trial)/N_RCT < 5)) stop("One of the trial's superpopulation size is less than 5 times the sample size per trial") # warning if propensity scores too extreme
  if (any(table(pop_init$trial)/N_RCT < 10)) warning("One of the trial's superpopulation size is less than 10 times the sample size per trial") # warning if propensity scores too extreme

  pop_BC <- pop_init[trial == "BC"][, ttt := rep_len(c("C", "B"), length.out = .N)] # 06/05/25 : no duplication of patients anymore: with a large enough general population and a not too imbalanced population, should be fine
  pop_AC <- pop_init[trial == "AC"][, ttt := rep_len(c("C", "A"), length.out = .N)]

  all_individuals <- data.table::rbindlist(list(pop_BC, pop_AC), use.names = TRUE)
  average_all_individuals <- all_individuals[, lapply(.SD, mean), .SDcols = covariate_names, by = trial]


  average_conditional_outcome_all_individuals <- sapply(list("A" = average_all_individuals[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                                                             "B" = average_all_individuals[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                                                             "C" = average_all_individuals[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
                                                        predict_outcome,
                                                        outcome_model = outcome_generation_formula,
                                                        simplify = FALSE) |>
    c("trial" = list(average_all_individuals$trial)) |>
    as.data.table() |>
    melt(measure.vars = c("A", "B", "C"), variable.name = "ttt", value.name = "outcome")
  # if (outcome_distribution == "binomial") average_conditional_outcome_all_individuals[, outcome :=  plogis(outcome)]

  marginal_outcome_all_individuals <- all_individuals[, lapply(.SD, mean), .SDcols = c("A", "B", "C"), by = c("trial")] |>
    data.table::melt(id.vars = "trial", measure.vars = c("A", "B", "C"), value.name = "outcome", variable.name = "ttt")
  if (outcome_distribution == "binomial") marginal_outcome_all_individuals <- marginal_outcome_all_individuals[ ,.(trial, ttt, outcome  = log(outcome/(1 - outcome)))] # setting back to linear scale to be able to estimate AB as A - B

  average_outcome_df <- rbindlist(
    list("conditional" = average_conditional_outcome_all_individuals,
         "marginal" = marginal_outcome_all_individuals),
    use.names = TRUE,
    idcol = "outcome_type") |>
    dcast(trial + outcome_type ~ ttt, value.var = "outcome")
  average_outcome_df[, AB := A - B] # linear scale

  # population_variance <- df_outcomes[, .(var_Y_obs = var(Y_obs)), by = c("ttt")] |>
  #   dcast(. ~ ttt, value.var = "var_Y_obs") |>
  #   dplyr::rename(var_population = `.`) |>
  #   dplyr::mutate(var_AB = A + B)
  #
  # df_outcomes |> ggplot() + geom_violin(aes(x = Y_obs, y = ttt))
  # diff_AB <- df_outcomes[ttt == "A", Y_obs] - df_outcomes[ttt == "B", Y_obs]
  # mean(diff_AB)
  # mean((diff_AB - mean(diff_AB))^2)
  #
  # average_outcome_df <- merge(average_outcome_df, population_variance, by = "ttt")

  # Correcting theoretical marginal effect, so that it is set to 0 when there is actually
  # no difference between theoretical conditional and marginal, as it should be
  # Useful to quantify estimators alpha and beta nominal risk level
  # if (average_outcome_df[outcome_type == "conditional", AB] == 0) average_outcome_df[, AB := 0]

  detach(list_simulation_parameters)
  return(list(
    "pop_init" = pop_init,
    "average_outcome_df" = average_outcome_df
  ))
}

# méthodes de comparaison indirecte
indirect_comparisons <- function(pop_init,
                                 struct_results,
                                 N_BOOT_ITER,
                                 N_RCT,
                                 outcome_regression_model,
                                 covariate_names,
                                 assignment_model,
                                 outcome_distribution,
                                 retrieve_ps_weights = FALSE) {

  #############################
  ############## Drawing trials
  #############################

  trial_AC <- pop_init[trial == "AC"][sample(1:.N, N_RCT, replace = TRUE), ][
    , ttt := rep_len(c("C", "A"), length.out = .N) |> factor(levels = c("C", "A"))]
  trial_BC <- pop_init[trial == "BC"][sample(1:.N, N_RCT, replace = TRUE), ][
    , ttt := rep_len(c("C", "B"), length.out = .N) |> factor(levels = c("C", "B"))]

  trial_AC$Y_obs <- with(trial_AC, dplyr::case_match(as.character(ttt), "A" ~ A, "B" ~ B, "C" ~ C))
  trial_BC$Y_obs <- with(trial_BC, dplyr::case_match(as.character(ttt), "A" ~ A, "B" ~ B, "C" ~ C))
  stopifnot(all(levels(trial_AC$ttt)[[1]] == "C",
                levels(trial_BC$ttt)[[1]] == "C"))

  is.binary <- function(x) length(unique(x) |> Filter(f = \(y) !is.na(y))) <= 2
  are_binary_covariates <- sapply(trial_AC[, ..covariate_names], is.binary)

  glm_family <- switch(outcome_distribution,
                       normal = gaussian(link = "identity"),
                       binomial = binomial(link = "logit"))

  stopifnot(names(trial_AC) == names(trial_BC))
  ###############################
  ########## Unadjusted estimator
  ###############################
  struct_results$unadjusted$anchored$unadjusted <- run_unadjusted_estimator(trial_AC, trial_BC, anchored = TRUE, glm_family)
  struct_results$unadjusted$unanchored$unadjusted <- run_unadjusted_estimator(trial_AC, trial_BC, anchored = FALSE, glm_family)

  ##############################################################
  ########## REGRESSION BASED OUTCOME MODEL (both treatment IPD)
  ##############################################################

  struct_results$regression$anchored$glm <- run_regression_model(trial_AC,
                                                                 trial_BC,
                                                                 outcome_regression_model,
                                                                 covariate_names,
                                                                 full_ipd = TRUE,
                                                                 anchored = TRUE,
                                                                 outcome_family = glm_family)
  struct_results$regression$unanchored$glm <-  run_regression_model(trial_AC,
                                                                    trial_BC,
                                                                    outcome_regression_model,
                                                                    covariate_names,
                                                                    full_ipd = TRUE,
                                                                    anchored = FALSE,
                                                                    outcome_family = glm_family)

  #################################################
  ########### PROPENSITY SCORE (both treatment IPD)
  #################################################

  struct_results$iptw$anchored$ml <- run_propensity_score(trial_AC,
                                                          trial_BC,
                                                          covariate_names,
                                                          assignment_model,
                                                          anchored = TRUE,
                                                          weight_estimation_method = "max_likelihood",
                                                          retrieve_ps_weights = retrieve_ps_weights,
                                                          outcome_family = glm_family,
                                                          are_binary_covariates = are_binary_covariates)
  struct_results$iptw$unanchored$ml <- run_propensity_score(trial_AC,
                                                            trial_BC,
                                                            covariate_names,
                                                            assignment_model,
                                                            anchored = FALSE,
                                                            weight_estimation_method = "max_likelihood",
                                                            retrieve_ps_weights = retrieve_ps_weights,
                                                            outcome_family = glm_family,
                                                            are_binary_covariates = are_binary_covariates)
  ##################
  ### MAIC Moments 1
  ##################
  struct_results$iptw$anchored$maic_1 <- run_propensity_score(trial_AC,
                                                            trial_BC,
                                                            covariate_names,
                                                            assignment_model,
                                                            anchored = TRUE,
                                                            weight_estimation_method = "moments_1",
                                                            retrieve_ps_weights = retrieve_ps_weights,
                                                            outcome_family = glm_family,
                                                            are_binary_covariates = are_binary_covariates)

  struct_results$iptw$unanchored$maic_1 <- run_propensity_score(trial_AC,
                                                              trial_BC,
                                                              covariate_names,
                                                              anchored = FALSE,
                                                              assignment_model,
                                                              weight_estimation_method = "moments_1",
                                                              retrieve_ps_weights = retrieve_ps_weights,
                                                              outcome_family = glm_family,
                                                              are_binary_covariates = are_binary_covariates)
  ##################
  ### MAIC Moments 1
  ##################
  struct_results$iptw$anchored$maic_2 <- run_propensity_score(trial_AC,
                                                            trial_BC,
                                                            covariate_names,
                                                            assignment_model,
                                                            anchored = TRUE,
                                                            weight_estimation_method = "moments_2",
                                                            retrieve_ps_weights = retrieve_ps_weights,
                                                            outcome_family = glm_family,
                                                            are_binary_covariates = are_binary_covariates)
  struct_results$iptw$unanchored$maic_2 <- run_propensity_score(trial_AC,
                                                              trial_BC,
                                                              covariate_names,
                                                              anchored = FALSE,
                                                              assignment_model,
                                                              weight_estimation_method = "moments_2",
                                                              retrieve_ps_weights = retrieve_ps_weights,
                                                              outcome_family = glm_family,
                                                              are_binary_covariates = are_binary_covariates)

  if (retrieve_ps_weights) {
    print("retrieving ps weights")
     list_dfs <- sapply(struct_results[["iptw"]],
                           \(sublist) sapply(sublist, \(subsublist) {
                             subsublist[["df"]]
                           }, simplify = FALSE, USE.NAMES = TRUE),
                           simplify = FALSE,
                           USE.NAMES = TRUE)
     dfs_wo_weights <- sapply(list_dfs, \(l) l$ml[, names(l$ml) != "trial_weights", with = FALSE], simplify = FALSE, USE.NAMES = TRUE)
     weights <- sapply(list_dfs, \(l) sapply(l, \(x) x$trial_weights, simplify = TRUE, USE.NAMES = TRUE), simplify = FALSE, USE.NAMES = TRUE)

     struct_ps_df <- mapply(FUN = cbind,
                            dfs_wo_weights,
                            weights,
                            SIMPLIFY = FALSE,
                            USE.NAMES = FALSE)

  }
  # Retrieving errors
  list_dfs_boot_errors <- sapply(struct_results[["iptw"]],
                     \(sublist) sapply(sublist, \(subsublist) {
                       subsublist[["boot_errors"]]
                     }, simplify = FALSE, USE.NAMES = TRUE),
                     simplify = FALSE,
                     USE.NAMES = TRUE)
  list_boot_warnings <- sapply(struct_results[["iptw"]],
                          \(sublist) sapply(sublist, \(subsublist) {
                            subsublist[["boot_warnings"]]
                          }, simplify = FALSE, USE.NAMES = TRUE),
                          simplify = FALSE,
                          USE.NAMES = TRUE)




  struct_results[["iptw"]] <- sapply(struct_results[["iptw"]],
                                     \(sublist) sapply(sublist,
                                                       \(subsublist) {
                                                         subsublist[["df"]] <- NULL
                                                         subsublist[["boot_errors"]] <- NULL
                                                         subsublist[["boot_warnings"]] <- NULL
                                                         subsublist[["estimate"]] <- ifelse(is.numeric(subsublist[["estimate"]]),
                                                                                            subsublist[["estimate"]],
                                                                                            NA)
                                                         subsublist[["variance"]] <- ifelse(is.numeric(subsublist[["variance"]]),
                                                                                            subsublist[["variance"]],
                                                                                            NA)
                                                         return(subsublist)
                                                       }, simplify = FALSE, USE.NAMES = TRUE),
                                     simplify = FALSE,
                                     USE.NAMES = TRUE)


  ##########
  ###### STC
  ##########
  struct_results$regression$anchored$stc <- run_regression_model(trial_AC,
                                                                 trial_BC,
                                                                 outcome_regression_model,
                                                                 covariate_names,
                                                                 full_ipd = FALSE,
                                                                 anchored = TRUE,
                                                                 outcome_family = glm_family)
  struct_results$regression$unanchored$stc <- run_regression_model(trial_AC,
                                                                   trial_BC,
                                                                   outcome_regression_model,
                                                                   covariate_names,
                                                                   full_ipd = FALSE,
                                                                   anchored = FALSE,
                                                                   outcome_family = glm_family)


  ##############################
  ############ COMPILING RESULTS
  ##############################

  rectangle_results <- struct_results |> tibble::enframe() |>
    tidyr::unnest_longer(value, indices_to = "anchored") |>
    tidyr::unnest_longer(value, indices_to = "model") |>
    tidyr::unnest_longer(value, indices_to = "indicator") |>
    dplyr::rename(adjustment = name) |>
    tidyr::pivot_wider(names_from = indicator, values_from = value)

  list_results <- list(rectangle_results = rectangle_results)
  if (retrieve_ps_weights) list_results$struct_ps_df <- struct_ps_df
  list_results$list_dfs_boot_errors <- list_dfs_boot_errors
  list_results$list_boot_warnings <- list_boot_warnings
  return(list_results)
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
    anchored = list(ml = NULL, maic_1 = NULL, maic_2 = NULL),
    unanchored = list(ml = NULL, maic_1 = NULL, maic_2 = NULL)
  )
)

###############
#### ESTIMATORS
###############

# Used to specify variables to use for "trial exposure" models, and unanchored STC
list_covariate_names <- c(
  combn(c("X1", "X2"), m = 1, simplify = FALSE),
  combn(c("X1", "X2"), m = 2, simplify = FALSE)
  # combn(c("X1", "X2", "X3", "X4"), m = 3, simplify = FALSE),
  # combn(c("X1", "X2", "X3", "X4"), m = 4, simplify = FALSE)
) |> as.vector()

quadratic = TRUE
if (quadratic == TRUE) {
  list_assignment_model <- sapply(
    list_covariate_names,
    \(li) sapply(li, \(x) paste0("poly(", x, ", 2, raw = TRUE)")) |> paste0(collapse = " + ")
  )
} else {
  list_assignment_model <- sapply(list_covariate_names, paste0, collapse = " + ")
}
list_outcome_regression_models <- sapply(
  list_covariate_names,
  \(li) sapply(li, \(x) paste0(x, "*ttt")) |> paste0(collapse = " + ")
)

df_estimators_parameters <- data.table(
  "covariate_names" = list_covariate_names,
  "assignment_model" = list_assignment_model,
  "outcome_regression_model" = list_outcome_regression_models
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


