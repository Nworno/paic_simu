library(data.table)

options(mc.cores = 1)
# Parameters
####################
N_pop <- 10^6

prop_X1 <- 0.5
bT_X1 <- 0.5
bT_X2 <- 0.2
bT_X3 <- -0.5
bT_X4 <- 0.3

bY_X1 <- 1.5
bY_X2 <- -0.5
bY_X3 <- 0.5
bY_X4 <- 2

bY_A_X1 <- 1.2
bY_A_X4 <- 0.2

bY_A <- 1.5
bY_B <- 1.5
bY_C <- 0

N_RCT <- 200
N_BOOT_ITER <- 200

list_simulation_parameters <- list(
  N_pop = N_pop,
  prop_X1 = prop_X1,
  bT_X1 = bT_X1,
  bY_X1 = bY_X1,
  bYA_X1 = bY_A_X1,
  bY_A = bY_A,
  bY_B = bY_B,
  bY_C = bY_C,
  N_RCT = N_RCT,
  N_BOOT_ITER = N_BOOT_ITER
)
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
trial_assignment_model_BC <- bquote(0) ## Note David: The value for BC trial doesn't matter here, as long as it is the same for all individuals, ensuring that the prevalence of X1 in the BC trial will be the same as the source population
list_outcome_generation_formula <- list(
  "1" = bquote(bY_A_X1*X1*A + bY_X1*X1 + bY_A*A + bY_B*B + bY_C*C),
  "2" = bquote(bY_X1*X1 + bY_X2 * X2 + (bY_A + bY_A_X1*X1) * A + bY_B*B + bY_C*C),
  "3" = bquote(bY_X1*X1 + bY_X2 * X2 + bY_X3 * X3 + bY_X4 * X4 + (bY_A + bY_A_X1*X1) * A + bY_B*B + bY_C*C),
  "4" = bquote(bY_X1*X1 + bY_X2 * X2 + bY_X3 * X3 + bY_X4 * X4 + (bY_A + bY_A_X4*X4) * A + bY_B*B + bY_C*C)
)

list_AC_trial_models <- list(
  "1" = bquote(X1 * bT_X1),
  "2" = bquote(X1 * bT_X1 + X2 * bT_X2),
  "3" = bquote(X1 * bT_X1 + X2 * bT_X2 + X3 * bT_X3 + X4*bT_X4),
  "4" = bquote(X1 * bT_X1 + X2 * bT_X2 + X3 * bT_X3 + X4*bT_X4)
)
list_BC_trial_models <- list(
  "1" = bquote(0)
)

# Used to specify variables to use for "trial exposure" models
list_covariate_names <- list(
  "1" = c("X1"),
  "2" = c("X1", "X2"),
  "3" = c("X1", "X2", "X3", "X4"),
  "4" = c("X1", "X2", "X3", "X4")
)
list_outcome_regression_models <- list(
  "1" = formula("Y_obs ~ X1*ttt"),
  "2" = formula("Y_obs ~ X1*ttt + X2"),
  "3" = formula("Y_obs ~ X1*ttt + X2 + X3 + X4"),
  "4" = formula("Y_obs ~ X1 + X2 + X3 + X4*ttt")
)

##############################################
########### Creating an overarching population
##############################################
pop_init <- tibble::tibble(
  id = 1:N_pop,
  X1 = rbinom(N_pop, 1, prop_X1),
  X2 = rnorm(N_pop, 0, 1),
  X3 = rlnorm(N_pop, 0, 0.5),
  binary_marker = rbinom(N_pop, 1, prop_X1),
  X4 = binary_marker * rnorm(N_pop, -1.5, 1) + (1 - binary_marker) * rnorm(N_pop, 1.5, 1)
) |>
  as.data.table() |>
  setkey("id")


trial_assignement_prob <- function(trial_assignment_model, df) {
  predicted <- with(df, eval(trial_assignment_model))
  return(plogis(predicted))
}

predict_outcome <- function(outcome_model, df) {
  with(df, eval(outcome_model))
}

## Note David: theoretical prevalence of X1 in the AC trial: prop_X1*plogis(bT_X1)/(prop_X1*plogis(bT_X1) + (1-prop_X1)*plogis(0))
df_proba_AC <- sapply(list_AC_trial_models,
                      trial_assignement_prob,
                      df = pop_init,
                      simplify = FALSE) |>
  c("id" = list(1:nrow(pop_init))) |>
  as.data.table() |>
  melt(id.vars = c("id"), variable.name = "trial_AC_model", value.name = "prob_w_trial_AC") |>
  setkey("id")
## Note David: theoretical prevalence of X1 in the BC trial: prop_X1*plogis(0)/(prop_X1*plogis(0) + (1-prop_X1)*plogis(0))
df_proba_BC <- sapply(list_BC_trial_models, trial_assignement_prob, df = pop_init, simplify = FALSE) |>
  c("id" = list(1:nrow(pop_init))) |>
  as.data.table() |>
  melt(id.vars = c("id"), variable.name = "trial_BC_model", value.name = "prob_w_trial_BC") |>
  setkey("id")


theoretical_outcomes_A <- lapply(list_outcome_generation_formula,
                                 predict_outcome,
                                 df = pop_init[, c("A", "B", "C") := .(1L, 0L, 0L)]) |>
  setNames(1:length(list_outcome_generation_formula)) |>
  c("ttt" = "A", "id" = list(1:nrow(pop_init))) |>
  as.data.table() |>
  melt(id.vars = c("id", "ttt"), variable.name = "outcome_model", value.name = "Y_theo")
theoretical_outcomes_B <- lapply(list_outcome_generation_formula,
                                 predict_outcome,
                                 df = pop_init[, c("A", "B", "C") := .(0L, 1L, 0L)]) |>
  setNames(1:length(list_outcome_generation_formula)) |>
  c("ttt" = "B", "id" = list(1:nrow(pop_init))) |>
  as.data.table() |>
  melt(id.vars = c("id", "ttt"), variable.name = "outcome_model", value.name = "Y_theo")
theoretical_outcomes_C <- lapply(list_outcome_generation_formula,
                                 predict_outcome,
                                 df = pop_init[, c("A", "B", "C") := .(0L, 0L, 1L)]) |>
  setNames(1:length(list_outcome_generation_formula)) |>
  c("ttt" = "C", "id" = list(1:nrow(pop_init))) |>
  as.data.table() |>
  melt(id.vars = c("id", "ttt"), variable.name = "outcome_model", value.name = "Y_theo")

pop_init[, c("A", "B", "C") := NULL]
df_outcomes <- rbindlist(list(theoretical_outcomes_A, theoretical_outcomes_B, theoretical_outcomes_C)) |>
  setkey("id")
df_outcomes[, Y_obs := Y_theo + rnorm(n = length(Y_theo), mean = 0, sd = 1)]

average_pop_init <- pop_init[, lapply(.SD, mean), .SDcols = c("X1", "X2", "X3", "X4")]


average_conditional_outcome_A <- sapply(list_outcome_generation_formula,
                                        predict_outcome,
                                        average_pop_init[,c("A", "B", "C") := .(1L, 0L, 0L)])


average_conditional_outcome_B <- sapply(list_outcome_generation_formula,
                                        predict_outcome,
                                        average_pop_init[,c("A", "B", "C") := .(0L, 1L, 0L)])


average_conditional_outcome_C <- sapply(list_outcome_generation_formula,
       predict_outcome,
       average_pop_init[,c("A", "B", "C") := .(0L, 0L, 1L)])

average_conditional_outcome <- data.frame("A" = average_conditional_outcome_A,
                                          "B" = average_conditional_outcome_B,
                                          "C" = average_conditional_outcome_C,
                                          "outcome_model" = names(list_outcome_generation_formula)) |>
  as.data.table() |>
  melt(measure.vars = c("A", "B", "C"), variable.name = "ttt", value.name = "outcome")
marginal_outcome <- df_outcomes[, .(outcome = mean(Y_obs)), by = c("ttt", "outcome_model")]

average_outcome_df <- rbindlist(
  list("conditional" = average_conditional_outcome,
       "marginal" = marginal_outcome),
  use.names = TRUE,
  idcol = "outcome_type") |>
  dcast(outcome_type + outcome_model ~ ttt, value.var = "outcome")
average_outcome_df[, AB := A - B]
################################
########## TRUE TREATMENT EFFECT
################################

# Marginal true effect in pop_init
# Works because the distribution of the initial population is similar to the
# distribution of the BC trial population, i.e. pop_init$prop_BC is constant
# true_marginal_AB <- pop_init[, mean(Y_theo_A)] - pop_init[, mean(Y_theo_B)]
#
# # (Average) Conditional true effect in pop_init (at the average levels of the confounders/treatment effect modifiers)
# true_conditional_AB <- eval(outcome_model, envir = list(X1 = prop_X1, A = 1, B = 0, C = 0)) -
#   eval(outcome_model, envir = list(X1 = prop_X1, A = 0, B = 1, C = 0))

################################
##### TREATMENT EFFECTS variance
################################
# Under a given treatment, the terms of the outcome model are independent,
# and the variance of the outcome is therefore the sum of the variances
# of the individual terms (taking into account that treatment is constant, therefore coefficients for a given treatment does not need to be taken into account)
# If there are dependencies between terms, then the outcome variance has to take into account covariances
# var_theo_A <- (bY_X1 + bYA_X1)^2 * (prop_X1) * (1 - prop_X1) # Var(aX) = a^2 Var(X)
# var_theo_B <- (bY_X1)^2 * (prop_X1) * (1 - prop_X1)
# var_theo_C <- (bY_X1)^2 * (prop_X1) * (1 - prop_X1)
#
# var_theo_AB <- var_theo_A/N_RCT + var_theo_B/N_RCT

###############################
########## Unadjusted estimator
##############################

# (Anchored) naive observed effect in pop_init
unadjusted_estimator <- function(trial_AC, trial_BC, anchored) {
  if (anchored) { # ie two steps
    naive_conditional_model_AC <- glm(Y_obs ~ ttt, family = "gaussian", data = trial_AC)
    naive_conditional_model_BC <- glm(Y_obs ~ ttt, family = "gaussian", data = trial_BC)
    naive_AB <- naive_conditional_model_AC$coefficients[["tttA"]] - naive_conditional_model_BC$coefficients[["tttB"]]
  } else {
    both_trials <- rbind(trial_AC[ttt == "A", .(Y_obs, ttt)], trial_BC[ttt == "B", .(Y_obs, ttt)])
    both_trials[ ,ttt := relevel(ttt, ref = "B")]
    naive_conditional_model_AB <- glm(Y_obs ~ ttt, family = "gaussian", data = both_trials)
    naive_AB <- naive_conditional_model_AB$coefficients[["tttA"]]
    # naive_variance_AB <- vcov(naive_conditional_model_AC)["tttA", "tttA"] + vcov(naive_conditional_model_BC)["tttB", "tttB"]
    ## Equivalent to
    # naive_AB <- trial_AC[ttt == "A", mean(Y_obs)] - trial_BC[ttt == "B", mean(Y_obs)]
    # var_naive_unanchored_AB <- trial_AC[ttt == "A", var(Y_obs) / .N] + trial_BC[ttt == "B", var(Y_obs) / .N]
  }
  return(naive_AB)
}

# (Unanchored) unadjusted observed effect in pop_init
run_unadjusted_estimator <- function(trial_AC, trial_BC, anchored) {
  estimate <- unadjusted_estimator(trial_AC, trial_BC, anchored)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) unadjusted_estimator(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                                    trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                                    anchored))
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var(na.rm = TRUE)
  return(list("estimate" = estimate, "variance" = variance))
}


##############################################################
########## REGRESSION BASED OUTCOME MODEL (both treatment IPD)
##############################################################

##### Anchored
# Anchored observed conditional effect with IPD (in population similar to BC trial)
# Of note, the model is provided with the true outcome model, easier
# different implementations methodologies of "anchored" conditional model

# centring covariates used in a model predicts the average outcome effect at the level of the value used to center the covariates

anchored_conditional_estimation <- function(centered_trial_AC, centered_trial_BC, regression_model, glm_family) {
  # two-steps individual patient data network meta analysis --> would be interesting to compare performance differences of this method
  # as compared to random effect NMA
  # potentially less biased, but systematically less precise as compared to "one-step unanchored" approach,
  # because takes into account unobserved confounding between A and C, and B and C
  # Equivalent of random effect meta analysis
  model_AC <- glm(regression_model, glm_family, data = centered_trial_AC) # adjusted conditional effect
  model_BC <- glm(regression_model, glm_family, data = centered_trial_BC) # adjusted conditional effect
  estimate_AC <- model_AC$coefficients[["tttA"]]
  estimate_BC <- model_BC$coefficients[["tttB"]]
  estimate_AB <- estimate_AC - estimate_BC

  ## Equivalent to
  # model_AC <- glm(formula(Y_obs ~ X1*ttt), glm_family, trial_AC)
  # estimate_AC <- model_AC$coefficients[["tttA"]] + model_AC$coefficients[["X1:tttA"]] * trial_BC[, mean(X1)]
  # model_BC <- glm(formula(Y_obs ~ X1*ttt), glm_family, data = trial_BC)
  # estimate_BC <- model_BC$coefficients[["tttB"]] + model_BC$coefficients[["X1:tttB"]] * trial_BC[, mean(X1)]
  # estimate_AB <- estimate_AC - estimate_BC
  return(estimate_AB)
}

run_anchored_conditional_estimation <- function(centered_trial_AC, trial_BC, outcome_regression_model, glm_family) {
  estimate <- anchored_conditional_estimation(centered_trial_AC, trial_BC, outcome_regression_model, glm_family)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) {
    anchored_conditional_estimation(
      centered_trial_AC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
      trial_BC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
      outcome_regression_model, gaussian)
  })
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}

### Unanchored
#### Note: not really following either anchored or unanchored scheme, basically
#### unanchored analysis when a common comparator is available: make the assumption
#### that there are no unobserved imbalance in terms of prognostic factors or
#### TEM between trials. Will be more precise in such a case, but biased in case
#### of residual confounding. Not used currently
botharms_unanchored_conditional_estimation <- function(centered_trial_AC, centered_trial_BC, regression_model, glm_family) {
  both_trials <- rbind(centered_trial_AC, centered_trial_BC)
  both_trials[, `:=`(ttt = factor(ttt, levels = c("B", "A", "C")))]
  model_AB <- glm(regression_model, data = both_trials, family = glm_family) # conditional effect
  estimate_AB <- model_AB$coefficients[["tttA"]]
}

unanchored_conditional_estimation <- function(centered_trial_AC, centered_trial_BC, outcome_regression_model, glm_family) {
  both_trials <- rbind(centered_trial_AC, centered_trial_BC)
  both_trials <- both_trials[ttt %in% c("A", "B")][, `:=`(ttt = factor(ttt, levels = c("B", "A")))]

  model_AB <- glm(outcome_regression_model, data = both_trials, family = glm_family) # conditional effect
  estimate_AB <- model_AB$coefficients[["tttA"]]

  # equivalent to
  # both_trials <- rbind(trial_AC, trial_BC)
  # both_trials[, `:=`(ttt = factor(ttt, levels = c("B", "A", "C")))]
  # model_AB <- glm(regression_model, data = both_trials) # conditional effect
  # estimate_AB <- model_AB$coefficients[["tttA"]] + model_AB$coefficients[["X1:tttA"]]*mean(trial_BC[["X1"]])

}


run_unanchored_conditional_estimation <- function(centered_trial_AC, centered_trial_BC, outcome_regression_model, glm_family) {
  estimate <- unanchored_conditional_estimation(centered_trial_AC, centered_trial_BC, outcome_regression_model, glm_family)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) unanchored_conditional_estimation(
    centered_trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
    centered_trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
    outcome_regression_model,
    glm_family
  ))
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}



########## PROPENSITY SCORE
propensity_score <- function(trial_AC, trial_BC, anchored, covariate_names) {
  df_anchored <- rbindlist(list("AC" = trial_AC, "BC" = trial_BC), idcol = "trial", fill = TRUE)
  df_anchored[, trial := as.factor(trial)]
  df_unanchored <- df_anchored[ttt %in% c("A", "B")]
  stopifnot(levels(df_anchored$trial)[[1]] == "AC")
  stopifnot(levels(df_unanchored$trial)[[1]] == "AC")

  trial_BC_assigment_model <- paste0("trial ~ ", paste0(covariate_names, collapse = " + "))

  # /!\ a major difference between anchored and unanchored here is that anchored has 2x more patients!
  if (anchored) {
    # predicting belonging to the AC trial
    df_anchored$PS_BC_trial <- glm(trial_BC_assigment_model, df_anchored, family = binomial(link = "logit"))$fitted.values
    df_anchored$ATC_w <- df_anchored[, (trial == "BC") + (trial == "AC") * PS_BC_trial / (1 - PS_BC_trial)]
  } else {
    df_unanchored$PS_B_trial <- glm(trial_BC_assigment_model, data = df_unanchored, family = binomial(link = "logit"))$fitted.values
    df_unanchored$ATC_w <- df_unanchored[, (trial == "BC") + (trial == "AC") * PS_B_trial / (1 - PS_B_trial)]
  }
  if (anchored) {
    # obs_marginal_AC <- df_anchored[trial == "AC" & ttt == "A", weighted.mean(Y_obs, ATC_w)] -
    #   df_anchored[trial == "AC" & ttt == "C", weighted.mean(Y_obs, ATC_w)]
    # obs_marginal_BC <- df_anchored[trial == "BC" & ttt == "B", mean(Y_obs)] -
    #   df_anchored[trial == "BC" & ttt == "C", mean(Y_obs)]
    # estimate_AB <- obs_marginal_AC - obs_marginal_BC
    # var_anchored_marginal_AB <- trials_combined[trial == "AC" & ttt == "A", var(Y_obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "AC" & ttt == "C", var(Y_obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "BC" & ttt == "B", var(Y_obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "BC" & ttt == "C", var(Y_obs * ATT_BC_w) / (.N * mean(ATT_BC_w))]
    # ### Alternative implementation
    model_marginal_AC <- glm(Y_obs ~ ttt, family = gaussian, data = df_anchored[trial == "AC"], weights = ATC_w)
    model_marginal_BC <- glm(Y_obs ~ ttt, family = gaussian, data = df_anchored[trial == "BC"])
    estimate_AB <- model_marginal_AC$coefficients[["tttA"]] - model_marginal_BC$coefficients[["tttB"]]
    # var_anchored_marginal_AB <- sandwich::vcovHC(model_marginal_AC)["tttA", "tttA"] +
    #   sandwich::vcovHC(model_marginal_BC)["tttB", "tttB"]

    # Unanchored implementation
    # model_marginal_AB <- glm(Y_obs ~ ttt, family = gaussian, data = df_anchored[, ttt := relevel(ttt, ref = "B")], weights = ATC_w)
    # estimate_AB <- model_marginal_AB$coefficients[["tttA"]]

  } else {
    estimate_AB <- df_unanchored[ttt == "A", weighted.mean(Y_obs, ATC_w)] - df_unanchored[ttt == "B", mean(Y_obs)]
    # var_unanchored_marginal_AB <- trials_combined[trial == "AC" & ttt == "A", var(Y_obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "BC" & ttt == "B", var(Y_obs * ATT_BC_w) / (.N * mean(ATT_BC_w))]
  }
  return(estimate_AB)
}

run_propensity_score <- function(trial_AC, trial_BC, anchored, trial_assignment_model) {
  estimate <- propensity_score(trial_AC, trial_BC, anchored, trial_assignment_model)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) propensity_score(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                          trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                          anchored,
                                                          trial_assignment_model))
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}

#########
### MAIC
#########

# MAIC
mm_obj_fun <- function(params, X) {
  sum(exp(X %*% params))
}

mm_grad_fun <- function(params, X) {
  colSums(sweep(X, 1, exp(X %*% params), FUN = "*"))
}

maic <- function(trial_AC, trial_BC, covariate_names, anchored, outcome_family) {
  if (anchored) {
    mean_trial_BC <- trial_BC[, lapply(.SD, mean), .SDcols = covariate_names] |> data.matrix()
    centered_covariates <- sweep(trial_AC[, ..covariate_names] |> data.matrix(),
                                 2,
                                 mean_trial_BC,
                                 FUN = "-")
  } else {
    mean_trial_B <- trial_BC[ttt == "B", lapply(.SD, mean), .SDcols = covariate_names] |> data.matrix()
    centered_covariates <- sweep(trial_AC[ttt == "A", ..covariate_names] |> data.matrix(),
                                 2,
                                 mean_trial_B,
                                 FUN = "-")
  }
  random_init <- rep(0, ncol(centered_covariates))
  params <- optim(random_init, fn = mm_obj_fun, gr = mm_grad_fun, method = "BFGS", X = centered_covariates)$par
  weights <- exp(centered_covariates %*% params)
  # ess <- sum(weights)^2 / sum(weights^2)
  if (anchored) {
    trial_AC$weights <- weights
    fitted_AC_w <- glm(Y_obs ~ ttt, data = trial_AC, family = outcome_family, weights = weights)
    anchored_MAIC_AC <- fitted_AC_w$coefficients[["tttA"]]
    # Equivalent to
    # trial_AC[ttt == "A", weighted.mean(Y_obs, weights)] - trial_AC[ttt == "C", weighted.mean(Y_obs, weights)]
    obs_marginal_BC <- trial_BC[ttt == "B", mean(Y_obs)] - trial_BC[ttt == "C", mean(Y_obs)]
    estimate_AB <- anchored_MAIC_AC - obs_marginal_BC
  } else {
    estimate_AB <- trial_AC[ttt == "A", weighted.mean(Y_obs, weights)] - trial_BC[ttt == "B", mean(Y_obs)]
  }
  return(estimate_AB)
}

run_maic <- function(trial_AC, trial_BC, covariate_names, anchored, outcome_family) {
  stopifnot(levels(trial_AC$ttt)[[1]] == "C")
  estimate <- maic(trial_AC, trial_BC, covariate_names, anchored, outcome_family)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) maic(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                              trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                              covariate_names,
                                              anchored,
                                              outcome_family))
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}

##########
###### STC
##########


stc <- function(centered_trial_AC, trial_BC, anchored, outcome_regression_model, covariate_names, outcome_family = gaussian) {
  if (anchored) {
    average_observed_BC <- trial_BC[ttt == "B", mean(Y_obs)] - trial_BC[ttt == "C", mean(Y_obs)]
    stc_model <- glm(outcome_regression_model,
                     family = outcome_family,
                     data = centered_trial_AC)
    simulated_AC <- stc_model$coefficients[["tttA"]]
    # Equivalent to
    # stc_model <- glm(Y_obs ~ ttt*X1, family = gaussian, data = trial_AC)
    # estimate_AC <- stc_model$coefficients[["tttA"]] + stc_model$coefficients[["tttA:X1"]] * Ag_trial_BC[, sum(mean_X1*n)/sum(n)]
    estimate_AB <- simulated_AC - average_observed_BC
  } else {
    outcome_regression_model <- paste0("Y_obs ~ ", paste0(covariate_names, collapse = " + "))
    average_observed_B <- trial_BC[ttt == "B", mean(Y_obs)]
    stc_model <- glm(outcome_regression_model,
                     family = outcome_family,
                     data = centered_trial_AC[ttt == "A",])
    simulated_A <- stc_model$coefficients[["(Intercept)"]]
    # Equivalent to
    # ttt_A <- trial_AC[ttt == "A"]
    # stc_model <- glm(Y_obs ~ X1, family = gaussian, data = ttt_A)
    # simulated_A <- stc_model$coefficients[["(Intercept)"]] + stc_model$coefficients[["X1"]] * trial_BC[ttt == "B", mean(X1)]
    estimate_AB <- simulated_A - average_observed_B
  }
  return(estimate_AB)
}

run_stc <- function(centered_trial_AC, trial_BC, anchored, outcome_regression_model, covariate_names, outcome_family) {
  estimate <- stc(centered_trial_AC, trial_BC, anchored, outcome_regression_model, covariate_names, outcome_family)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) stc(centered_trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                   trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                   anchored,
                                                   outcome_regression_model,
                                                   covariate_names,
                                                   outcome_family))
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}


#############################
############## MC SIMULATIONS
#############################

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

comparison <- function(pop_init,
                       df_proba_AC,
                       df_proba_BC,
                       df_outcomes,
                       struct_results,
                       N_RCT,
                       N_BOOT_ITER,
                       outcome_generation_formula_number,
                       AC_trial_model_number,
                       BC_trial_model_number,
                       outcome_regression_model,
                       covariate_names) {

  #############################
  ############## Drawing trials
  #############################

  # Could do all the draws at once, and transform to long df with bind_rows, and then do one massive join
  selected_individuals_AC <- df_proba_AC[trial_AC_model == AC_trial_model_number,][
    sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_AC)][
      , ttt := rep_len(c("A", "C"), length.out = .N)]
  selected_outcomes_AC <- df_outcomes[selected_individuals_AC, on = c("id", "ttt")][
    outcome_model == outcome_generation_formula_number, c("id", "ttt", "outcome_model", "Y_obs")]
  trial_AC <- pop_init[selected_outcomes_AC, on = "id"][, ttt := factor(ttt, levels = c("C", "A"))]

  selected_individuals_BC <- df_proba_BC[trial_BC_model == BC_trial_model_number,][
    sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_BC)][
      , ttt := rep_len(c("B", "C"), length.out = .N)]
  selected_outcomes_BC <- df_outcomes[selected_individuals_BC, on = c("id", "ttt")][
    outcome_model == outcome_generation_formula_number, c("id", "ttt", "outcome_model", "Y_obs")]
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
  # |>
  #   dplyr::bind_rows(data.frame(adjustment = "true", model = c("conditional", "marginal"), anchored = NA, estimate = c(true_conditional_AB, true_marginal_AB), variance = var_theo_AB))
  return(rectangle_results)
}

n_iter = 300
list_experiments <- tibble::tibble(
  outcome_generation_formula = list_outcome_generation_formula,
  outcome_generation_formula_number = names(list_outcome_generation_formula),
  outcome_regression_model = list_outcome_regression_models,
  AC_trial_model = list_AC_trial_models,
  AC_trial_model_number = names(list_AC_trial_models),
  BC_trial_model = list_BC_trial_models,
  BC_trial_model_number = names(list_BC_trial_models),
  covariate_names = list_covariate_names
)

time_start <- Sys.time()
time_start_string <- format(time_start, "%Y%m%d_%H%M%S")
experiment_results_directory <- file.path("results_simulations", time_start_string)
if (!dir.exists(experiment_results_directory)) dir.create(experiment_results_directory)

for (row in 1:nrow(list_experiments)) {
  outcome_regression_model <- list_experiments[row, "outcome_regression_model", drop = TRUE][[1]]
  outcome_generation_formula_number <- list_experiments[row, "outcome_generation_formula_number", drop = TRUE][[1]]
  AC_trial_model_number <- list_experiments[row, "AC_trial_model_number", drop = TRUE][[1]]
  BC_trial_model_number <- list_experiments[row, "BC_trial_model_number", drop = TRUE][[1]]
  covariate_names <- list_experiments[row, "covariate_names", drop = TRUE][[1]]
  results_simulations <- parallel::mclapply(1:n_iter, \(i) {
    time_start_iteration <- Sys.time()
    result_comparison <- comparison(pop_init,
                                    df_proba_AC,
                                    df_proba_BC,
                                    df_outcomes,
                                    struct_results,
                                    N_RCT,
                                    N_BOOT_ITER,
                                    outcome_generation_formula_number,
                                    AC_trial_model_number,
                                    BC_trial_model_number,
                                    outcome_regression_model,
                                    covariate_names)
    time_eluded <- Sys.time() - time_start_iteration
    cat("Experiment ", row, ", Iteration ", i, ", length: ", time_eluded, "seconds\n")
    return(result_comparison)
  })

  saveRDS(results_simulations, file = file.path(experiment_results_directory,
                                                paste0("experiment_", row, ".RDS")))

}
  cat("Simulation length ", Sys.time() - time_start, "seconds \n")
