library(data.table)

##############################################
########### Data Generation Parameters
##############################################

df_default_parameters <- list(
  N_RCT  = c(500),
  bT_X1 = 1,   # Effet de la variable binaire sur la probabilité d'être dans l'essai BC
  bT2_X1 = -0.5,   # Effet de la variable binaire sur la probabilité d'être dans l'essai BC
  bT_X2 = 1,   # Effet de la variable continue X2...
  bT2_X2 = -0.5,   # Effet de la variable continue X2...
  fbT = 1,
  f_X1 = bquote(stop("f_X1 must be specified for each scenario")),
  f_X2 = bquote(stop("f_X2 must be specified for each scenario")),
  bY_X1 = 1,   # Effet de X1 sur l'outcome
  bY_X2 = 2,   # Effet de X2 sur l'outcome
  bY_A_X1 = c(2), # Interaction A et X1 dans le modèle outcome
  bY_A_X1_2 = c(0), # Interaction quadratique A et X1 dans le modèle outcome
  bY_A_X2 = c(0), # Interaction A et X2 dans le modèle outcome
  bY_B_X1 = c(0), # Interaction A et X1 dans le modèle outcome
  bY_B_X2 = c(0), # Interaction A et X2 dans le modèle outcome
  bY_A = 1,  # Effet de A par rapport à C
  bY_B = 1,  # Effet de B par rapport à C
  bY_C = 0,    # Pas d'effet de C sur l'outcome
  BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2)), # Modèle d'attribution de l'essai BC
  outcome_distribution = c("normal"),
  outcome_generation_formula =  c(bquote(
    bY_X1*X1 + bY_X2 * X2 + (bY_A + bY_A_X1*X1 + bY_A_X1_2*X1^2 + bY_A_X2*X2) * A +  (bY_B + bY_B_X1*X1 + bY_B_X2*X2) *B + bY_C*C
  ))
)


list_changing_parameters <- list(
  "1" = list(
    f_X1 = c(bquote(rnorm(N_pop, 0, 1))),
    f_X2 = c(bquote(rnorm(N_pop, 0, 1))),
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
    bT_X1 = 1,
    bT2_X1 = -0.5,
    bT_X2 = 1,
    bT2_X2 = -0.5,
    fbT = -1,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
    ),
  "3" = list(
    f_X1 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    f_X2 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    bT_X1 = 2,
    bT2_X1 = 0,
    bT_X2 = 2,
    bT2_X2 = 0,
    fbT = -1,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "4" = list(
    f_X1 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    f_X2 = c(bquote(pmin(rlnorm(N_pop, 0, 0.5), 5))),
    bT_X1 = 2,
    bT2_X1 = 0,
    bT_X2 = 2,
    bT2_X2 = 0,
    fbT = 1,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "5" = list(
    f_X1 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    f_X2 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    bT_X1 = 1,
    bT2_X1 = -1,
    bT_X2 = 1,
    bT2_X2 = -1,
    fbT = 0.5,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  "6" = list(
    f_X1 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    f_X2 = c(bquote(sample(c(rnorm(N_pop/2, 0, 1), rnorm(N_pop/2, 3, 1))))),
    bT_X1 = 1,
    bT2_X1 = -1,
    bT_X2 = 1,
    bT2_X2 = -1,
    fbT = -0.5,
    BC_trial_model = c(bquote(bT + bT_X1 * X1 + bT2_X1 * X1^2 + bT_X2 * X2 + bT2_X2 * X2^2))
  ),
  # DGM-7 : interaction quadratique dans l'essai AC - bon overlap
  "7" = list(
    f_X1 = c(bquote(rnorm(N_pop, 0, 1))),
    f_X2 = c(bquote(rnorm(N_pop, 0, 1))),
    bY_A_X1_2 = 1
  ),

  # DGM-8 : interaction quadratique dans l'essai AC - mauvais overlap
  "8" = list(
    f_X1 = c(bquote(rnorm(N_pop, 0, 1))),
    f_X2 = c(bquote(rnorm(N_pop, 0, 1))),
    bY_A_X1_2 = 1,
    fbT = -1
  ),

  # DGM-9 : outcome binaire - bon overlap
  "9" = list(
    f_X1 = c(bquote(rnorm(N_pop, 0, 1))),
    f_X2 = c(bquote(rnorm(N_pop, 0, 1))),
    outcome_distribution = "binomial",
    fbT = 1
  ),

  # DGM-10 : outcome binaire - mauvais overlap
  "10" = list(
    f_X1 = c(bquote(rnorm(N_pop, 0, 1))),
    f_X2 = c(bquote(rnorm(N_pop, 0, 1))),
    outcome_distribution = "binomial",
    fbT = -1
  )
)

list_parameters <- lapply(list_changing_parameters, \(x) {
  list_parameter <- c(x, df_default_parameters[!names(df_default_parameters) %in% names(x)], recursive = TRUE)
  list_parameter <- list_parameter[names(df_default_parameters)]
})

df_population_parameters <- list_parameters |> tibble::as_tibble() |> t()
colnames(df_population_parameters) <- names(df_default_parameters)
df_population_parameters <- data.table::as.data.table(df_population_parameters)
df_population_parameters[, population_parameters_num := as.integer(names(list_changing_parameters))]


##############################################
########### Creating an overarching population
##############################################

creating_population <- function(list_simulation_parameters, N_pop) {
  attach(list_simulation_parameters)
  # Modifying the coefficient values using fbT
  bT_X1 <- bT_X1*fbT
  bT2_X1 <- bT2_X1*fbT
  bT_X2 <- bT_X2*fbT
  bT2_X2 <- bT2_X2*fbT

  pop_init <- data.table(
    id = 1:N_pop,
    X1 = eval(f_X1),
    X2 = eval(f_X2)
  ) |>
    setkey("id")

  # Finding out bT value which provides balanced probabilities
  bT = 0
  qPtrial <- with(pop_init, eval(BC_trial_model))
  fct <- function(x, qP, prev) {
    mean(plogis(x+qP)-prev)
  }
  tryCatch({
    bT <- uniroot(fct, interval = c(mean(qPtrial) - 100, mean(qPtrial) + 100), qP = qPtrial, prev = 0.5)$root
  }, error = function(e) {
    print("Error in uniroot")
  })

  trial_assignement_prob <- function(trial_assignment_model, df) {
    predicted <- with(df, eval(trial_assignment_model))
    return(plogis(predicted))
  }

  predict_outcome <- function(outcome_model, df) {
    with(df, eval(outcome_model))
  }

  covariate_names <- c("X1", "X2")

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
    df_outcomes_pop_init[, Y_obs := rbinom(n = .N, size = 1, prob = plogis(Y_theo))]
    moy_outcomes <- tapply(df_outcomes_pop_init$Y_obs, df_outcomes_pop_init$ttt, mean, simplify = FALSE)
    if (any(moy_outcomes < 0.02 | moy_outcomes > 0.98)) { # arbitrary thresholds, to avoid downstream problems
      stop("Too extreme outcomes")
    }
  } else {
    stop("Wrong outcome_distribution parameter value")
  }


  pop_init <- df_outcomes_pop_init[, Y_theo:= NULL][pop_init, on = "id"] |> data.table::dcast(formula = ... ~ ttt, value.var = "Y_obs")

  pop_init[, prob_BC := trial_assignement_prob(BC_trial_model, df = pop_init)]
  pop_init[, trial := rbinom(.N, 1, prob_BC) |>
             factor(levels = c(0, 1), labels = c("AC", "BC"))]

  print(table(pop_init$trial)/N_RCT)
  if (any(table(pop_init$trial)/N_RCT < 5)) stop("One of the trial's superpopulation size is less than 5 times the sample size per trial") # warning if propensity scores too extreme
  if (any(table(pop_init$trial)/N_RCT < 10)) warning("One of the trial's superpopulation size is less than 10 times the sample size per trial") # warning if propensity scores too extreme

  pop_BC <- pop_init[trial == "BC"][, ttt := rep_len(c("C", "B"), length.out = .N)]
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

  marginal_outcome_all_individuals <- all_individuals[, lapply(.SD, mean), .SDcols = c("A", "B", "C"), by = c("trial")] |>
    data.table::melt(id.vars = "trial", measure.vars = c("A", "B", "C"), value.name = "outcome", variable.name = "ttt")
  # setting back to linear scale to be able to estimate AB as A - B
  if (outcome_distribution == "binomial") marginal_outcome_all_individuals <- marginal_outcome_all_individuals[ ,.(trial, ttt, outcome  = log(outcome/(1 - outcome)))]

  average_outcome_df <- rbindlist(
    list("conditional" = average_conditional_outcome_all_individuals,
         "marginal" = marginal_outcome_all_individuals),
    use.names = TRUE,
    idcol = "outcome_type") |>
    dcast(trial + outcome_type ~ ttt, value.var = "outcome")
  average_outcome_df[, AB := A - B] # linear scale

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

  # struct_results$regression$anchored$glm <- run_regression_model(trial_AC,
  #                                                                trial_BC,
  #                                                                outcome_regression_model,
  #                                                                covariate_names,
  #                                                                full_ipd = TRUE,
  #                                                                anchored = TRUE,
  #                                                                outcome_family = glm_family)
  # struct_results$regression$unanchored$glm <-  run_regression_model(trial_AC,
  #                                                                   trial_BC,
  #                                                                   outcome_regression_model,
  #                                                                   covariate_names,
  #                                                                   full_ipd = TRUE,
  #                                                                   anchored = FALSE,
  #                                                                   outcome_family = glm_family)

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
  ### MAIC Moments 2
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
                                                         # Convert condition objects to their message string so that pivot_wider
                                                         # does not create list-columns (which would break as.numeric on ess/estimate/variance)
                                                         err <- subsublist[["error_estimate"]]
                                                         subsublist[["error_estimate"]] <- if (inherits(err, "condition")) conditionMessage(err) else NA_character_
                                                         return(subsublist)
                                                       }, simplify = FALSE, USE.NAMES = TRUE),
                                     simplify = FALSE,
                                     USE.NAMES = TRUE)


  ##########
  ###### STC
  ##########
  # struct_results$regression$anchored$stc <- run_regression_model(trial_AC,
  #                                                                trial_BC,
  #                                                                outcome_regression_model,
  #                                                                covariate_names,
  #                                                                full_ipd = FALSE,
  #                                                                anchored = TRUE,
  #                                                                outcome_family = glm_family)
  # struct_results$regression$unanchored$stc <- run_regression_model(trial_AC,
  #                                                                  trial_BC,
  #                                                                  outcome_regression_model,
  #                                                                  covariate_names,
  #                                                                  full_ipd = FALSE,
  #                                                                  anchored = FALSE,
  #                                                                  outcome_family = glm_family)


  ##############################
  ############ COMPILING RESULTS
  ##############################

  rectangle_results <- struct_results |> tibble::enframe() |>
    tidyr::unnest_longer(value, indices_to = "anchored") |>
    tidyr::unnest_longer(value, indices_to = "model") |>
    dplyr::rename(adjustment = name) |>
    dplyr::mutate(
      estimate       = sapply(value, \(x) { e <- x[["estimate"]];       if (is.null(e) || !is.numeric(e)) NA_real_      else as.numeric(e) }),
      variance       = sapply(value, \(x) { e <- x[["variance"]];       if (is.null(e) || !is.numeric(e)) NA_real_      else as.numeric(e) }),
      ess            = sapply(value, \(x) { e <- x[["ess"]];            if (is.null(e))                    NA_real_      else as.numeric(e) }),
      error_estimate = sapply(value, \(x) { e <- x[["error_estimate"]]; if (is.null(e))                    NA_character_ else as.character(e) })
    ) |>
    dplyr::select(-value)

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
  # regression = list(
  #   anchored = list(glm = NULL, stc = NULL),
  #   unanchored = list(glm = NULL, stc = NULL)
  # ),
  iptw = list(
    anchored = list(ml = NULL, maic_1 = NULL, maic_2 = NULL),
    unanchored = list(ml = NULL, maic_1 = NULL, maic_2 = NULL)
  )
)

# Used to specify variables to use for "trial exposure" models, and unanchored STC
list_covariate_names <- c(
  combn(c("X1", "X2"), m = 1, simplify = FALSE),
  combn(c("X1", "X2"), m = 2, simplify = FALSE)
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
