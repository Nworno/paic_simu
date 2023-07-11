library(ggplot2)
library(data.table)
library(dplyr)
library(tictoc)
library(fastglm)

ggthemr::ggthemr("flat")
logit <- function(x) log(x/(1 - x))

######## Parameters
####################
N_pop <- 10^6
N_RCT <- 200

prop_X1 <- 0.5
bT_X1 <- 2

bY_X1 <- 1.5
bYA_X1 <- 1.2

bY_A <- 1.5
bY_B <- 1.5
bY_C <- 0

######### Models
##################
trial_assignment_model_AC <- expr(X1 * bT_X1)
trial_assignment_model_BC <- expr(0) # required that prob of trial assignment independent of baseline characteristics of overarching population for the ATT to represent the marginal effect in the initial (overall) population
outcome_model <- expr(bYA_X1*X1*A + bY_X1*X1 + bY_A*A + bY_B*B + bY_C*C)

########### Creating an overarching population
##############################################
pop_init <- data.table(
  id = 1:N_pop,
  X1 = rbinom(N_pop, 1, prop_X1)
) |>
  setkey("id")

##### TREATMENT EFFECTS variance
#############################
# Under a given treatment, treatment assignment is independent on baseline characteristics
# by definition. Therefore, the terms of the outcome model are independent,
# and the variance of the total is the sum of the variances
# of the individual terms (that do vary, that is treatment effect do not influence the variance of the outcome under a given treatment)
# If there are dependencies between terms, then the outcome variance has to take into account covariances
# these are the variances of the outcome under the different treatments
# Var theoretical outcomes
var_theo_A <- (bY_X1 + bYA_X1)^2 * (prop_X1) * (1 - prop_X1)
var_theo_B <- (bY_X1)^2 * (prop_X1) * (1 - prop_X1)
var_theo_C <- (bY_X1)^2 * (prop_X1) * (1 - prop_X1) # variance of the outcome under treatment C without deviates


# var_theo_AB <- var_theo_A + var_theo_B
trial_assignement_prob <- function(df, trial_assignment_model, deviates = FALSE) {
  predicted <- with(df, eval(trial_assignment_model))
  # variance of a logistic function is set ((s^2\pi^2)/3), with s=1 when using normal law -->
  # making an assumption about bias and residual heterogeneity
  return(plogis(predicted))
}

predict_outcome <- function(df, outcome_model, deviates = FALSE) {
  predicted <- with(df, eval(outcome_model))
  if (deviates) return(predicted + rnorm(n = length(predicted), mean = 0, sd = 1)) else return(predicted)
}

pop_init$prob_AC <- trial_assignement_prob(pop_init, trial_assignment_model_AC)
pop_init$prob_BC <- trial_assignement_prob(pop_init, trial_assignment_model_BC)


pop_init$Y_theo_A <- copy(pop_init)[, c("A", "B", "C") := .(1L, 0L, 0L)] |>
  predict_outcome(outcome_model)
pop_init$Y_theo_B <- copy(pop_init)[, c("A", "B", "C") := .(0L, 1L, 0L)] |>
  predict_outcome(outcome_model)
pop_init$Y_theo_C <- copy(pop_init)[, c("A", "B", "C") := .(0L, 0L, 1L)] |>
  predict_outcome(outcome_model)

ttt_names <- c("A", "B", "C")
cols_obs <- paste0("Y_obs_", ttt_names)
cols_theo <- paste0("Y_theo_", ttt_names)
# Assigning "observed" values for all the individual in the population
pop_init[, eval(expr(cols_obs)) := lapply(.SD, \(x) x + rnorm(n = length(x))), .SD = cols_theo]
# equivalent
# pop_init[, sub("theo", "obs", cols_theo, fixed = TRUE) := lapply(.SD, \(x) x + rnorm(n = length(x))), .SD = cols_theo]


pop_init_long <- pop_init |> tidyr::pivot_longer(
  cols = starts_with("Y_"),
  names_to = c("measure_type", "ttt"),
  names_prefix = "Y_",
  names_sep = "_"
) |>
  tidyr::pivot_wider(names_from = "measure_type", values_from = "value") |>
  # C ref levels for subsequent models
  dplyr::mutate(ttt = factor(ttt, levels = c("C", "A", "B"))) |>
  data.table::setDT()
# may prefer alternative straightforward implementation using data.table
# pop_init_long <- melt(pop_init,
#                       id.vars = "id",
#                       measure.vars = patterns("theo", "obs"),
#                       value.name = c("theo", "obs"),
#                       variable.name = "Y",
#                       variable.factor = FALSE
# )

################################
########## TRUE TREATMENT EFFECT
################################

# Marginal true effect in pop_init
# Works because the distribution of the initial population is similar to the
# distribution of the BC trial population, i.e. pop_init[, mean(prop_BC)] = 0.5
true_marginal_AB <- pop_init[, mean(Y_theo_A)] - pop_init[, mean(Y_theo_B)]

# (Average) Conditional true effect in pop_init (at the average levels of the confounders/treatment effect modifiers)
true_conditional_AB <- eval(outcome_model, envir = list(X1 = prop_X1, A = 1, B = 0, C = 0)) -
  eval(outcome_model, envir = list(X1 = prop_X1, A = 0, B = 1, C = 0))

############## Drawing trials
#############################
# Could do all the draws at once, and transform to long df with bind_rows, and then do one massive join

trial_AC <- pop_init[
  sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_AC), # use integer based indexing
  .(id, ttt = factor(rep_len(c("A", "C"), length.out = N_RCT), levels = c("C", "A", "B")))][
    pop_init_long[,c("id", "X1", "ttt", "prob_AC", "theo", "obs")], on = .(id, ttt), nomatch = NULL
  ]
trial_BC <- pop_init[
  sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_BC), # use integer based indexing
  .(id, ttt = factor(rep_len(c("B", "C"), length.out = N_RCT), levels = c("C","A", "B")))][
    pop_init_long[,c("id", "X1", "ttt", "prob_BC", "theo", "obs")], on = .(id, ttt), nomatch = NULL]
trials_combined <- data.table::rbindlist(list("AC" = trial_AC, "BC" = trial_BC), fill = TRUE, idcol = "trial")[, trial := factor(trial, levels = c("AC", "BC"))] |> setkey("id")
# trials_combined[, X1_centred :=]


##############################################################
########## REGRESSION BASED OUTCOME MODEL (both treatment IPD)
##############################################################

##### Anchored
# Anchored observed conditional effect with IPD (in population similar to BC trial)
# Of note, the model is provided with the true outcome model, easier
  # different implementations methodologies of "anchored" conditional model
  # 1. two separate models, comparing the main treatment effects: adjusted network meta analysis using "aggregated" data
  # 2. two separate models, comparing the average outcome at the mean baseline values of the covariates from one of the trials: adjusted network meta analysis using individual data
  # 3. one model using all the data, comparing the main effects only (fixed effect regression): adjusted network meta analysis using individual data
  # 4. one model using all the data, comparing the predictions ... WIP (random effect model): one-step adjusted network meta analysis using individual data, assumes heterogeneity in treatment effects, but heterogeneity which will not be quantified, so doesn't make sense if the goal is to estimate this heterogeneity
# there is a need for average conditional treatment effect
# Methodological choice: predicting using the average conditional treatment effect, and not just centring the covariates
# calculating the variance of the estimate by bootstrapping, easier this way

conditional_model <- function(trial, outcome_model, glm_family, ttt_name, ttt_interaction_name, average_interaction_level) {
  model <- glm(formula = outcome_model, family = glm_family, data = trial)
  # model <- fastglm::fastglm(formula = outcome_model, family = glm_family, data = trial)
  average_ttt_estimate <- coefficients(model)[ttt_name] + coefficients(model)[ttt_interaction_name] * average_interaction_level ### /!\ using only the main effect, so parameters need to be centred, equivalent to predict when parameters value are all 0
  # var_ttt_estimate <- vcov(model)[ttt_name, ttt_name]
  # + vcov(model)[ttt_interaction_name, ttt_interaction_name] * average_interaction_level ### /!\ do we need the variance of the interaction too? not used by Phillippo
  return(list("ttt_estimate" = average_ttt_estimate))
}

twosteps_anchored_conditional_effect <- function(trial_AC, trial_BC, outcome_model, glm_family) {
  # two-steps individual patient data meta analysis
  model_AC <- conditional_model(trial_AC, outcome_model, glm_family, "tttA", "X1:tttA", trial_BC[, mean(X1)])
  model_BC <- conditional_model(trial_BC, outcome_model, glm_family, "tttB", "X1:tttB", trial_BC[, mean(X1)])
  estimate_AB <- model_AC$ttt_estimate - model_BC$ttt_estimate
  # var_estimate_AB <- model_AC$var_ttt_estimate + model_BC$var_ttt_estimate
  return(list("estimate_AB" = estimate_AB))
}


estimate <- twosteps_anchored_conditional_effect(trial_AC, trial_BC, formula(obs ~ X1*ttt), gaussian)$estimate
variance <- sapply(1:300, \(x) twosteps_anchored_conditional_effect(trial_AC[sample(1:.N, size = .N, replace = TRUE),],
                                                                    trial_BC[sample(1:.N, size = .N, replace = TRUE),],
                                                                    formula(obs ~ X1*ttt), gaussian)$estimate_AB) |>
  var()



onestep_anchored_conditional_effect <- function(outcome_model, glm_family, data) {
  # one-step meta-analysis with or without accounting for trial clusterization (which can be with fixed or random effect)
  # /!\ requires that reference level for ttt variable is B!
  # In the absence of residual confounding in the DGM (that is all confounders are correctly included in the adjustment model)
  # trial-level clusterization shouldn't change the results
  #
  # outcome_model_w_cluster <- update.formula(test, ~ . + trial)
  model <- conditional_model(data, outcome_model, glm_family, "tttA", "X1:tttA", data[trial == "BC", mean(X1)])
  estimate_AB <- model$ttt_estimate
  return(list("estimate_AB" = estimate_AB))
}

trials_combined[, ttt := relevel(ttt, ref = "B")]
stopifnot(levels(trials_combined)[[1]] == "B")
estimate <- onestep_anchored_conditional_effect(formula(obs ~ X1*ttt), gaussian, data = trials_combined)$estimate_AB
variance <- sapply(1:300, \(x)
                   onestep_anchored_conditional_effect(
                     formula(obs ~ X1*ttt),
                     gaussian,
                     data = trials_combined[sample(1:.N, size = .N, replace = TRUE), ]
                     )$estimate_AB) |>
  var()



# Same but with mixed model (ie random effect meta analysis on individual patients)
# mixed_model <- lme4::lmer(formula = obs ~ X1*ttt + (1|trial), data = trials_combined)
# fixed_coeff <- lme4::fixef(mixed_model)
# obs_conditional_AC <- fixed_coeff["tttA"] + fixed_coeff["X1:tttA"] * trials_combined[trial == "BC", mean(X1)]
# obs_conditional_BC <- fixed_coeff["tttB"] + fixed_coeff["X1:tttB"] * trials_combined[trial == "BC", mean(X1)]
# obs_conditional_AB <- obs_conditional_AC - obs_conditional_BC

################## TESTS
# (model_unanchored <- glm(formula = obs ~ X1*ttt, family = "gaussian", data = trials_combined))
# (model_unanchored <- glm(formula = obs ~ X1 + ttt, family = "gaussian", data = trials_combined))
# (model_unanchored <- glm(formula = obs ~ I(X1 - 3) + ttt, family = "gaussian", data = trials_combined))
# (model_unanchored <- glm(formula = obs ~ I(X1 - trials_combined[trial == "BC", mean(X1)])*ttt, family = "gaussian", data = trials_combined))
# (model_unanchored <- glm(formula = obs ~ I(X1 - 1)*ttt, family = "gaussian", data = trials_combined))
# (model_unanchored <- glm(formula = obs ~ I(X1 - 1)*ttt, family = "gaussian", data = trials_combined))
### Graphical exploration
# Why offseting covariates change only some treatment estimates (ie treatment effectin that case)
# ggplot(trials_combined, mapping = aes(I(X1 -1), obs)) +
#   geom_point(aes(color = ttt)) +
#   # facet_wrap(facets = c("ttt")) +
#   geom_smooth(aes(color = ttt, fill = ttt), method = "lm")
###################################

##### Unanchored

unanchored_conditional_effect <- function(outcome_model, glm_family, trial_AC, trial_BC) {
  observed_B <- trial_BC[ttt == "B", mean(obs)]
  model_A <- glm(outcome_model, family = glm_family, data = trial_AC[ttt == "A", ])
  # predicted A at the average level of predictors in treatment arm B
  predicted_A <- coefficients(model_A)[["(Intercept)"]] + coefficients(model_A)[["X1"]] * trial_BC[ttt == "B", mean(X1)]
  return(list("estimate_AB" = predicted_A - observed_B))
}

estimate <- unanchored_conditional_effect(formula(obs ~ X1), gaussian, trial_AC, trial_BC)$estimate_AB
variance <- sapply(1:300, \(x) unanchored_conditional_effect(
  formula(obs ~ X1),
  gaussian,
  trial_AC[sample(1:.N, .N, replace = TRUE), ],
  trial_BC[sample(1:.N, .N, replace = TRUE), ])$estimate_AB
  ) |>
  var()


#################################################
########### PROPENSITY SCORE (both treatment IPD)
#################################################

# (Anchored) observed marginal effect (distribution of X1 similar to the observed BC trial)

propensity_score <- function(df = trials_combined, anchored) {
  df$PS_trial <- glm(trial ~ X1, df, family = binomial(link = "logit"))$fitted.values
  df[, AT_BC_w := (trial == "BC") + (trial == "AC") * PS_trial / (1 - PS_trial)]
  if (anchored) {
    # obs_marginal_AC <- df[trial == "AC" & ttt == "A", sum(obs * AT_BC_w) / sum(AT_BC_w)] - df[trial == "AC" & ttt == "C", sum(obs * AT_BC_w) / sum(AT_BC_w)]
    # obs_marginal_BC <- df[trial == "BC" & ttt == "B", mean(obs)] - df[trial == "BC" & ttt == "C", mean(obs)]
    # estimate_AB <- obs_marginal_AC - obs_marginal_BC
    # var_anchored_marginal_AB <- trials_combined[trial == "AC" & ttt == "A", var(obs * AT_BC_w) / (.N * mean(AT_BC_w))] + trials_combined[trial == "AC" & ttt == "C", var(obs * AT_BC_w) / (.N * mean(AT_BC_w))] + trials_combined[trial == "BC" & ttt == "B", var(obs * AT_BC_w) / (.N * mean(AT_BC_w))] + trials_combined[trial == "BC" & ttt == "C", var(obs * AT_BC_w) / (.N * mean(AT_BC_w))]
    # ### Alternative implementation
    model_marginal_AC <- glm(obs ~ ttt, family = gaussian, data = df[trial == "AC"], weights = AT_BC_w)
    model_marginal_BC <- glm(obs ~ ttt, family = gaussian, data = df[trial == "BC"])
    estimate_AB <- model_marginal_AC$coefficients[["tttA"]] - model_marginal_BC$coefficients[["tttB"]]
    # var_anchored_marginal_AB <- sandwich::vcovHC(model_marginal_AC)["tttA", "tttA"] +
    #   sandwich::vcovHC(model_marginal_BC)["tttB", "tttB"]
  } else {
    estimate_AB <- df[trial == "AC" & ttt == "A", sum(obs * AT_BC_w) / sum(AT_BC_w)] -
      df[trial == "BC" & ttt == "B", mean(obs)]
    # var_unanchored_marginal_AB <- trials_combined[trial == "AC" & ttt == "A", var(obs * AT_BC_w) / (.N * mean(AT_BC_w))] + trials_combined[trial == "BC" & ttt == "B", var(obs * AT_BC_w) / (.N * mean(AT_BC_w))]
  }
  return(list("estimate_AB" = estimate_AB))
}

trials_combined[, ttt := relevel(ttt, ref = "C")]
stopifnot(levels(trials_combined)[[1]] == "C")
estimate <- propensity_score(trials_combined, anchored = TRUE)$estimate_AB
variance <- sapply(1:300, \(x) propensity_score(trials_combined[sample(1:.N, .N, replace = TRUE)], anchored = TRUE)$estimate_AB) |> var()

estimate <- propensity_score(trials_combined, anchored = FALSE)$estimate_AB
variance <- sapply(1:300, \(x) propensity_score(trials_combined[sample(1:.N, .N, replace = TRUE)], anchored = FALSE)$estimate_AB) |> var()


##############################################
########## Naive estimator
##############################################

# (Anchored) naive observed (or marginal, identical in that case) effect in pop_init

naive_estimator <- function(trial_AC, trial_BC, anchored = TRUE) {
  if (anchored) {
    naive_conditional_model_AC <- glm(obs ~ ttt, family = "gaussian", data = trial_AC)
    naive_conditional_model_BC <- glm(obs ~ ttt, family = "gaussian", data = trial_BC)
    naive_AB <- naive_conditional_model_AC$coefficients[["tttA"]] - naive_conditional_model_BC$coefficients[["tttB"]]
    # naive_variance_AB <- vcov(naive_conditional_model_AC)["tttA", "tttA"] + vcov(naive_conditional_model_BC)["tttB", "tttB"]
  } else {
    naive_AB <- trial_AC[ttt == "A", mean(obs)] - trial_BC[ttt == "B", mean(obs)]
    # var_naive_unanchored_AB <- trial_AC[ttt == "A", var(obs) / .N] + trial_BC[ttt == "B", var(obs) / .N]
  }
  return(list("estimate_AB" = naive_AB))
}

# (Unanchored) unadjusted observed effect in pop_init

estimate <- naive_estimator(trial_AC, trial_BC, anchored = TRUE)
variance <- sapply(1:300, \(x) naive_estimator(trial_AC[sample(1:.N, .N, replace = TRUE), ],
                                               trial_BC[sample(1:.N, .N, replace = TRUE), ],
                                               anchored = TRUE)$estimate_AB) |> var()
estimate <- naive_estimator(trial_AC, trial_BC, anchored = FALSE)
variance <- sapply(1:300, \(x) naive_estimator(trial_AC[sample(1:.N, .N, replace = TRUE), ],
                                               trial_BC[sample(1:.N, .N, replace = TRUE), ],
                                               anchored = FALSE)$estimate_AB) |> var()

###########################
### MAIC
###########################

Ag_trial_BC <- trial_BC[, .("mean_X1" = mean(X1), "sd_X1" = sd(X1), "n" = length(id),
                             "mean_theo" = mean(theo), "var_theo" = var(theo),
                             "mean_obs" = mean(obs), "var_obs" = var(obs)), by = ttt]

# MAIC
mm_obj_fun <- function(params, X) {
  sum(exp(X %*% params))
}

mm_grad_fun <- function(params, X) {
  colSums(sweep(X, 1, exp(X %*% params), FUN = "*"))
}

maic <- function(trial_AC, Ag_trial_BC, anchored = TRUE) {
  names_covariates <- c("X1")
  if (anchored) {
    mean_trial_BC <- Ag_trial_BC[, sum(mean_X1*n) / sum(n)]
    centered_covariates <- sweep(trial_AC[, ..names_covariates], 2, mean_trial_BC, FUN = "-")
  } else {
    mean_ttt_B <- Ag_trial_BC[ttt == "B", mean_X1]
    centered_covariates <- sweep(trial_AC[ttt == "A", ..names_covariates], 2, mean_ttt_B, FUN = "-")
  }
  random_init <- rep(0, ncol(centered_covariates))
  params <- optim(random_init, fn = mm_obj_fun, gr = mm_grad_fun, method = "BFGS", X = data.matrix(centered_covariates))$par
  weights <- exp(data.matrix(centered_covariates) %*% matrix(params, nrow = 1))
  # ess <- sum(weights)^2 / sum(weights^2)
  if (anchored) {
    trial_AC$weights <- weights
    fitted_AC_w <- glm(obs ~ ttt, data = trial_AC, family = gaussian, weights = weights)
    anchored_MAIC_AC <- fitted_AC_w$coefficients["tttA"]
    obs_marginal_BC <- Ag_trial_BC[ttt == "B", mean_obs] - Ag_trial_BC[ttt == "C", mean_obs]
    estimate_AB <- anchored_MAIC_AC - obs_marginal_BC
    # design <- survey::svydesign(ids = ~ 1, weights = weights, data = trial_AC)
    # svymodel <- survey::svyglm(formula = obs ~ ttt, design = design)
    #
    # (var_adjusted_marginal_MAIC_AC <- sandwich::vcovHC(fitted_AC_w, type = "HC3")["tttA", "tttA"])
    # (var(trial_AC[ttt == "A", "obs"])/nrow(trial_AC[ttt == "A", ]) + var(trial_BC[ttt == "B", "obs"])/nrow(trial_BC[ttt == "B", ]))
    # (vcov(svymodel)["tttA", "tttA"])
  } else {
    ttt_A <- trial_AC[ttt == "A", .(weights = weights, obs)]
    estimate_AB <- ttt_A[, sum(obs*weights / sum(weights))] - Ag_trial_BC[ttt == "B", mean_obs]
# var_unanchored_marginal_MAIC_AB <- trial_AC[ttt == "A", var(obs*weights) / (.N * mean(weights))] + trial_BC[ttt == "B", var(obs)/.N]
  }
  return(list("estimate_AB" = estimate_AB))
}

estimate <- maic(trial_AC, Ag_trial_BC, anchored = TRUE)$estimate_AB
variance <- sapply(1:300, \(x) maic(trial_AC[sample(1:.N, .N, replace = TRUE), ],
                                    Ag_trial_BC,
                                    anchored = TRUE)$estimate_AB) |>
  var() + Ag_trial_BC[ttt == "B", var_obs/n] + Ag_trial_BC[ttt == "C", var_obs/n]

estimate <- maic(trial_AC, Ag_trial_BC, anchored = FALSE)$estimate_AB
variance <- sapply(1:300, \(x) maic(trial_AC[sample(1:.N, .N, replace = TRUE), ],
                                    Ag_trial_BC,
                                    anchored = FALSE)$estimate_AB) + Ag_trial_BC[ttt == "B", var_obs/n] |>
  var() + Ag_trial_BC[ttt == "B", var_obs/n]


##########
###### STC
##########

average_conditional_BC <- conditional_model(trial_BC, formula(obs ~ X1*ttt), gaussian, "tttB", "X1:tttB", trial_BC[, mean(X1)])$ttt_estimate

stc <- function(trial_AC, Ag_trial_BC, anchored, average_conditional_BC = NULL) {
  if (anchored) {
    trial_AC$X1_centred <- trial_AC$X1 - Ag_trial_BC[, sum(mean_X1 * n) / sum(n)]
    stc_model <- glm(obs ~ ttt*X1_centred, family = gaussian, data = trial_AC)
    estimate_AC <- stc_model$coefficients[["tttA"]]
    estimate_AB <- estimate_AC - average_conditional_BC
    # Equivalent to
    # stc_model <- glm(obs ~ ttt*X1, family = gaussian, data = trial_AC)
    # estimate_AC <- stc_model$coefficients[["tttA"]] + stc_model$coefficients[["tttA:X1"]] * Ag_trial_BC[, sum(mean_X1*n)/sum(n)]
  } else {
    ttt_A <- trial_AC[, X1_centred := X1 - Ag_trial_BC[ttt == "B", mean_X1]][ttt == "A"]
    stc_model <- glm(obs ~ X1_centred, family = gaussian, data = ttt_A)
    simulated_A <- stc_model$coefficients[["(Intercept)"]]
    estimate_AB <- simulated_A - Ag_trial_BC[ttt == "B", mean_obs]
    # Equivalent to
    # ttt_A <- trial_AC[ttt == "A"]
    # stc_model <- glm(obs ~ X1, family = gaussian, data = ttt_A)
    # simulated_A <- stc_model$coefficients[["(Intercept)"]] + stc_model$coefficients[["X1"]] * trial_BC[ttt == "B", mean(X1)]
  }
  return(list("estimate_AB" = estimate_AB))
}

estimate <- stc(trial_AC, Ag_trial_BC, anchored = TRUE, average_conditional_BC = average_conditional_BC)
variance <- sapply(1:300, \(x) stc(trial_AC[sample(1:.N, .N, replace = TRUE), ],
                                   Ag_trial_BC,
                                   anchored = TRUE,
                                   average_conditional_BC = average_conditional_BC)$estimate_AB) |>
  var() + Ag_trial_BC[ttt == "B", var_obs/n] + Ag_trial_BC[ttt == "C", var_obs/n]


estimate <- stc(trial_AC, Ag_trial_BC, anchored = FALSE)$estimate_AB
variance <- sapply(1:300, \(x) stc(trial_AC[sample(1:.N, .N, replace = TRUE), ], Ag_trial_BC, anchored = FALSE)$estimate_AB) |>
  var() + Ag_trial_BC[ttt == "B", var_obs/n]


#####################
### Gathering results
#####################
#TODO: start from here
