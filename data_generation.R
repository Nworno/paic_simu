library(data.table)

# Parameters
####################
N_pop <- 10^6

prop_X1 <- 0.5
bT_X1 <- 2

bY_X1 <- 1.5
bYA_X1 <- 1.2

bY_A <- 1.5
bY_B <- 1.5
bY_C <- 0

N_RCT <- 200
N_BOOT_ITER <- 100
######### Models
##################
names_covariates <- c("X1")

trial_assignment_model_AC <- bquote(X1 * bT_X1)
trial_assignment_model_BC <- bquote(0) # required that prob of trial assignment independent of baseline characteristics of overarching population for the ATT to represent the marginal effect in the initial (overall) population
outcome_model <- bquote(bYA_X1*X1*A + bY_X1*X1 + bY_A*A + bY_B*B + bY_C*C)

##############################################
########### Creating an overarching population
##############################################
pop_init <- data.table(
  id = 1:N_pop,
  X1 = rbinom(N_pop, 1, prop_X1)
) |>
  setkey("id")


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


pop_init$Y_theo_A <- pop_init[, c("A", "B", "C") := .(1L, 0L, 0L)] |>
  predict_outcome(outcome_model)
pop_init$Y_theo_B <- pop_init[, c("A", "B", "C") := .(0L, 1L, 0L)] |>
  predict_outcome(outcome_model)
pop_init$Y_theo_C <- pop_init[, c("A", "B", "C") := .(0L, 0L, 1L)] |>
  predict_outcome(outcome_model)
pop_init[, c("A", "B", "C") := NULL]

ttt_names <- c("A", "B", "C")
cols_obs <- paste0("Y_obs_", ttt_names)
cols_theo <- paste0("Y_theo_", ttt_names)
# Assigning "observed" values for all the individual in the population
pop_init[, eval(bquote(cols_obs)) := lapply(.SD, \(x) x + rnorm(n = length(x))), .SD = cols_theo]
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

################################
##### TREATMENT EFFECTS variance
################################
# Under a given treatment, the terms of the outcome model are independent,
# and the variance of the outcome is therefore the sum of the variances
# of the individual terms (taking into account that treatment is constant, therefore coefficients for a given treatment does not need to be taken into account)
# If there are dependencies between terms, then the outcome variance has to take into account covariances
var_theo_A <- (bY_X1 + bYA_X1)^2 * (prop_X1) * (1 - prop_X1) # Var(aX) = a^2 Var(X)
var_theo_B <- (bY_X1)^2 * (prop_X1) * (1 - prop_X1)
var_theo_C <- (bY_X1)^2 * (prop_X1) * (1 - prop_X1)

var_theo_AB <- var_theo_A/N_RCT + var_theo_B/N_RCT

###############################
########## Unadjusted estimator
##############################

# (Anchored) naive observed effect in pop_init
unadjusted_estimator <- function(trial_AC, trial_BC, names_covariates, anchored) {
  if (anchored) { # ie two steps
    naive_conditional_model_AC <- glm(obs ~ ttt, family = "gaussian", data = trial_AC)
    naive_conditional_model_BC <- glm(obs ~ ttt, family = "gaussian", data = trial_BC)
    naive_AB <- naive_conditional_model_AC$coefficients[["tttA"]] - naive_conditional_model_BC$coefficients[["tttB"]]
  } else {
    df <- rbind(trial_AC[, .SD, .SDcols = c("obs", "ttt", names_covariates)],
                trial_BC[, .SD, .SDcols = c("obs", "ttt", names_covariates)])
    df[ ,ttt := relevel(ttt, ref = "B")]
    naive_conditional_model_AB <- glm(obs ~ ttt, family = "gaussian", data = df)
    naive_AB <- naive_conditional_model_AB$coefficients[["tttA"]]
    # naive_variance_AB <- vcov(naive_conditional_model_AC)["tttA", "tttA"] + vcov(naive_conditional_model_BC)["tttB", "tttB"]
    ## Equivalent to
    # naive_AB <- trial_AC[ttt == "A", mean(obs)] - trial_BC[ttt == "B", mean(obs)]
    # var_naive_unanchored_AB <- trial_AC[ttt == "A", var(obs) / .N] + trial_BC[ttt == "B", var(obs) / .N]
  }
  return(naive_AB)
}

# (Unanchored) unadjusted observed effect in pop_init
run_unadjusted_estimator <- function(trial_AC, trial_BC, names_covariates, anchored) {
  estimate <- unadjusted_estimator(trial_AC, trial_BC, names_covariates, anchored)
  boot_estimates <- parallel::mclapply(1:N_BOOT_ITER, \(x) unadjusted_estimator(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                                               trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                                               names_covariates,
                                                                               anchored), 
  mc.cores = 1L)
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
# 2. two separate models, comparing the average outcome at the mean baseline values of the covariates from one of the trials: adjusted network meta analysis using individual data
# 4. one model using all the data, comparing the predictions ... WIP (random effect model): one-step adjusted network meta analysis using individual data, assumes heterogeneity in treatment effects, but heterogeneity which will not be quantified, so doesn't make sense if the goal is to estimate this heterogeneity

# centring covariates used in a model predicts the average outcome effect at the level of the value used to center the covariates

anchored_conditional_estimation <- function(trial_AC, trial_BC, outcome_model, glm_family) {
  # two-steps individual patient data meta analysis
  model_AC <- glm(outcome_model, glm_family, data = trial_AC[, X1_centred := X1 - trial_BC[, mean(X1)]])
  estimate_AC <- model_AC$coefficients[["tttA"]]
  model_BC <- glm(outcome_model, glm_family, data = trial_BC[, X1_centred := X1 - trial_BC[, mean(X1)]])
  estimate_BC <- model_BC$coefficients[["tttB"]]
  estimate_AB <- estimate_AC - estimate_BC

  ## Equivalent to
  # model_AC <- glm(formula(obs ~ X1*ttt), glm_family, trial_AC)
  # estimate_AC <- model_AC$coefficients[["tttA"]] + model_AC$coefficients[["X1:tttA"]] * trial_BC[, mean(X1)]
  # model_BC <- glm(formula(obs ~ X1*ttt), glm_family, data = trial_BC)
  # estimate_BC <- model_BC$coefficients[["tttB"]] + model_BC$coefficients[["X1:tttB"]] * trial_BC[, mean(X1)]
  # estimate_AB <- estimate_AC - estimate_BC
  return(estimate_AB)
}

run_anchored_conditional_estimation <- function(trial_AC, trial_BC, outcome_model, glm_family) {
  estimate <- anchored_conditional_estimation(trial_AC, trial_BC, outcome_model, glm_family)
  boot_estimates <- parallel::mclapply(1:N_BOOT_ITER, \(x) {
    anchored_conditional_estimation(
      trial_AC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
      trial_BC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
      outcome_model, gaussian)
  }, mc.cores = 1L)
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}


### Unanchored

unanchored_conditional_estimation <- function(trial_AC, trial_BC, outcome_model, glm_family) {
  observed_B <- trial_BC[ttt == "B", mean(obs)] # unadjusted B effect
  # equivalent to (with collapsible outcome at least)
  # model_B <- glm(outcome_model, family = glm_family, data = trial_BC[ttt == "B", ])
  # observed_B <- model_B$coefficients[["(Intercept)"]] + model_B$coefficients[["X1"]] * trial_BC[ttt == "B", mean(X1)]

  model_A <- glm(outcome_model,
                 family = gaussian,
                 data = trial_AC[ttt == "A"][, X1_centred := X1 - trial_BC[ttt == "B", mean(X1)]])
  predicted_A <- coefficients(model_A)[["(Intercept)"]]
  # equivalent to
  # model_A <- glm(obs ~ X1, family = glm_family, data = trial_AC[ttt == "A", ])
  # predicted_A <- coefficients(model_A)[["(Intercept)"]] + coefficients(model_A)[["X1"]] * trial_BC[ttt == "B", mean(X1)] # make the assumption that the effect modification is identical for all treatments
  estimate_AB <- predicted_A - observed_B

  #### Equivalent to
  # one-step meta-analysis without accounting for trial clusterization (which can be with fixed or random effect)
  # In the absence of residual confounding in the DGM (that is all confounders are correctly included in the adjustment model)
  # trial-level clusterization shouldn't change the results
  # outcome_model_w_cluster <- update.formula(test, ~ . + trial)
  #
  # may allow taking into account the fact that all treatment modifications may not be the same?
  # stopifnot(levels(trials_combined$ttt)[[1]] == "B")
  # model <- glm(obs ~ X1_centred*ttt, glm_family, trials_combined[ttt %in% c("A", "B")][, X1_centred := X1 - trial_BC[ttt == "B", mean(X1)]])
  # estimate_AB <- model$coefficients[["tttA"]]
  # Equivalent to
  # model <- glm(obs ~ X1 * ttt, glm_family, data)
  # estimate_AB <- model$coefficients[["tttA"]] + model$coefficients[["X1:tttA"]] * data[trial == "BC", mean(X1)]

  return(estimate_AB)
}

run_unanchored_conditional_estimation <- function(trial_AC, trial_BC, outcome_model, glm_family) {
  estimate <- unanchored_conditional_estimation(trial_AC, trial_BC, outcome_model, glm_family)
  boot_estimates <- parallel::mclapply(1:N_BOOT_ITER, \(x) unanchored_conditional_estimation(
    trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
    trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
    outcome_model,
    gaussian
  ), mc.cores = 1L)
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}



########## PROPENSITY SCORE
propensity_score <- function(trial_AC, trial_BC, anchored) {
  df_anchored <- rbindlist(list("AC" = trial_AC, "BC" = trial_BC), idcol = "trial", fill = TRUE)
  df_anchored[, trial := as.factor(trial)]
  df_unanchored <- df_anchored[ttt %in% c("A", "B")]
  stopifnot(levels(df_anchored$trial)[[1]] == "AC")
  stopifnot(levels(df_unanchored$trial)[[1]] == "AC")
  # /!\ a major difference between anchored and unanchored here is that anchored has 2x more patients!
  if (anchored) {
    # predicting belonging to the AC trial
    df_anchored$PS_BC_trial <- glm(trial ~ X1, df_anchored, family = binomial(link = "logit"))$fitted.values
    df_anchored$ATC_w <- df_anchored[, (trial == "BC") + (trial == "AC") * PS_BC_trial / (1 - PS_BC_trial)]
  } else {
    df_unanchored$PS_B_trial <- glm(trial ~ X1, data = df_unanchored, family = binomial(link = "logit"))$fitted.values
    df_unanchored$ATC_w <- df_unanchored[, (trial == "BC") + (trial == "AC") * PS_B_trial / (1 - PS_B_trial)]
  }
  if (anchored) {
    # obs_marginal_AC <- df_anchored[trial == "AC" & ttt == "A", weighted.mean(obs, ATC_w)] -
    #   df_anchored[trial == "AC" & ttt == "C", weighted.mean(obs, ATC_w)]
    # obs_marginal_BC <- df_anchored[trial == "BC" & ttt == "B", mean(obs)] -
    #   df_anchored[trial == "BC" & ttt == "C", mean(obs)]
    # estimate_AB <- obs_marginal_AC - obs_marginal_BC
    # var_anchored_marginal_AB <- trials_combined[trial == "AC" & ttt == "A", var(obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "AC" & ttt == "C", var(obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "BC" & ttt == "B", var(obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "BC" & ttt == "C", var(obs * ATT_BC_w) / (.N * mean(ATT_BC_w))]
    # ### Alternative implementation
    model_marginal_AC <- glm(obs ~ ttt, family = gaussian, data = df_anchored[trial == "AC"], weights = ATC_w)
    model_marginal_BC <- glm(obs ~ ttt, family = gaussian, data = df_anchored[trial == "BC"])
    estimate_AB <- model_marginal_AC$coefficients[["tttA"]] - model_marginal_BC$coefficients[["tttB"]]
    # var_anchored_marginal_AB <- sandwich::vcovHC(model_marginal_AC)["tttA", "tttA"] +
    #   sandwich::vcovHC(model_marginal_BC)["tttB", "tttB"]
  } else {
    estimate_AB <- df_unanchored[ttt == "A", weighted.mean(obs, ATC_w)] - df_unanchored[ttt == "B", mean(obs)]
    # var_unanchored_marginal_AB <- trials_combined[trial == "AC" & ttt == "A", var(obs * ATT_BC_w) / (.N * mean(ATT_BC_w))] + trials_combined[trial == "BC" & ttt == "B", var(obs * ATT_BC_w) / (.N * mean(ATT_BC_w))]
  }
  return(estimate_AB)
}

run_propensity_score <- function(trial_AC, trial_BC, anchored) {
  estimate <- propensity_score(trial_AC, trial_BC, anchored)
  boot_estimates <- parallel::mclapply(1:N_BOOT_ITER, \(x) propensity_score(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                          trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                          anchored), mc.cores = 1L)
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

maic <- function(trial_AC, trial_BC, names_covariates, anchored) {
  if (anchored) {
    mean_X1_trial_BC <- trial_BC[, mean(X1)]
    centered_covariates <- sweep(trial_AC[, ..names_covariates], 2, mean_X1_trial_BC, FUN = "-")
  } else {
    mean_X1_ttt_B <- trial_BC[ttt == "B", mean(X1)]
    centered_covariates <- sweep(trial_AC[ttt == "A", ..names_covariates], 2, mean_X1_ttt_B, FUN = "-")
  }
  random_init <- rep(0, ncol(centered_covariates))
  params <- optim(random_init, fn = mm_obj_fun, gr = mm_grad_fun, method = "BFGS", X = data.matrix(centered_covariates))$par
  weights <- exp(data.matrix(centered_covariates) %*% matrix(params, nrow = 1))
  # ess <- sum(weights)^2 / sum(weights^2)
  if (anchored) {
    trial_AC$weights <- weights
    fitted_AC_w <- glm(obs ~ ttt, data = trial_AC, family = gaussian, weights = weights)
    anchored_MAIC_AC <- fitted_AC_w$coefficients[["tttA"]]
    # Equivalent to
    # trial_AC[ttt == "A", weighted.mean(obs, weights)] - trial_AC[ttt == "C", weighted.mean(obs, weights)]
    obs_marginal_BC <- trial_BC[ttt == "B", mean(obs)] - trial_BC[ttt == "C", mean(obs)]
    estimate_AB <- anchored_MAIC_AC - obs_marginal_BC
  } else {
    estimate_AB <- trial_AC[ttt == "A", weighted.mean(obs, weights)] - trial_BC[ttt == "B", mean(obs)]
  }
  return(estimate_AB)
}

run_maic <- function(trial_AC, trial_BC, names_covariates, anchored) {
  stopifnot(levels(trial_AC$ttt)[[1]] == "C")
  estimate <- maic(trial_AC, trial_BC, names_covariates, anchored = anchored)
  boot_estimates <- parallel::mclapply(1:N_BOOT_ITER, \(x) maic(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                              trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                              names_covariates,
                                              anchored = anchored), mc.cores = 1L)
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}

##########
###### STC
##########


stc <- function(trial_AC, trial_BC, anchored) {
  if (anchored) {
    average_observed_BC <- trial_BC[ttt == "B", mean(obs)] - trial_BC[ttt == "C", mean(obs)]
    stc_model <- glm(obs ~ ttt*X1_centred,
                     family = gaussian,
                     data = trial_AC[, X1_centred := X1 - trial_BC[, mean(X1)]])
    estimate_AC <- stc_model$coefficients[["tttA"]]
    # Equivalent to
    # stc_model <- glm(obs ~ ttt*X1, family = gaussian, data = trial_AC)
    # estimate_AC <- stc_model$coefficients[["tttA"]] + stc_model$coefficients[["tttA:X1"]] * Ag_trial_BC[, sum(mean_X1*n)/sum(n)]
    estimate_AB <- estimate_AC - average_observed_BC
  } else {
    average_observed_B <- trial_BC[ttt == "B", mean(obs)]
    stc_model <- glm(obs ~ X1_centred,
                     family = gaussian,
                     data = trial_AC[ttt == "A"][, X1_centred := X1 - trial_BC[ttt == "B", mean(X1)]])
    simulated_A <- stc_model$coefficients[["(Intercept)"]]
    # Equivalent to
    # ttt_A <- trial_AC[ttt == "A"]
    # stc_model <- glm(obs ~ X1, family = gaussian, data = ttt_A)
    # simulated_A <- stc_model$coefficients[["(Intercept)"]] + stc_model$coefficients[["X1"]] * trial_BC[ttt == "B", mean(X1)]
    estimate_AB <- simulated_A - average_observed_B
  }
  return(estimate_AB)
}

run_stc <- function(trial_AC, trial_BC, anchored) {
  estimate <- stc(trial_AC, trial_BC, anchored)
  boot_estimates <- parallel::mclapply(1:N_BOOT_ITER, \(x) stc(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                             trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                             anchored), mc.cores = 1L)
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

comparison <- function(pop_init, struct_results, N_RCT, N_BOOT_ITER) {

  #############################
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

  # Not used currently
  Ag_trial_BC <- rbind(
    trial_BC[, .("mean_X1" = mean(X1), "var_X1" = var(X1), "n" = length(id),
                 "mean_theo" = mean(theo), "var_theo" = var(theo),
                 "mean_obs" = mean(obs), "var_obs" = var(obs)), by = ttt],
    trial_BC[, .("ttt" = "overall",
                 "mean_X1" = mean(X1), "var_X1" = var(X1), "n" = length(id),
                 "mean_theo" = mean(theo), "var_theo" = var(theo),
                 "mean_obs" = mean(obs), "var_obs" = var(obs))]
  )

  ###############################
  ########## Unadjusted estimator
  ###############################
  struct_results$unadjusted$anchored$unadjusted <- run_unadjusted_estimator(trial_AC, trial_BC, names_covariates, anchored = TRUE)
  struct_results$unadjusted$unanchored$unadjusted <- run_unadjusted_estimator(trial_AC, trial_BC, names_covariates, anchored = FALSE)

  ##############################################################
  ########## REGRESSION BASED OUTCOME MODEL (both treatment IPD)
  ##############################################################

  struct_results$regression$anchored$glm <- run_anchored_conditional_estimation(trial_AC, trial_BC, formula(obs ~ X1_centred*ttt), gaussian)
  struct_results$regression$unanchored$glm <- run_unanchored_conditional_estimation(trial_AC, trial_BC, formula(obs ~ X1_centred), gaussian)

  #################################################
  ########### PROPENSITY SCORE (both treatment IPD)
  #################################################

  struct_results$iptw$anchored$ml <- run_propensity_score(trial_AC, trial_BC, anchored = TRUE)
  struct_results$iptw$unanchored$ml <- run_propensity_score(trial_AC, trial_BC, anchored = FALSE)

  #########
  ### MAIC
  #########
  struct_results$iptw$anchored$maic <- run_maic(trial_AC, trial_BC, names_covariates, anchored = TRUE)
  struct_results$iptw$unanchored$maic <- run_maic(trial_AC, trial_BC, names_covariates, anchored = FALSE)

  ##########
  ###### STC
  ##########

  struct_results$regression$anchored$stc <- run_stc(trial_AC, trial_BC, anchored = TRUE)
  struct_results$regression$unanchored$stc <- run_stc(trial_AC, trial_BC, anchored = FALSE)


  ##############################
  ############ COMPILING RESULTS
  ##############################
  rectangle_results <- struct_results |> tibble::enframe() |>
    tidyr::unnest_longer(value, indices_to = "anchored") |>
    tidyr::unnest_longer(value, indices_to = "model") |>
    tidyr::unnest_longer(value, indices_to = "indicator") |>
    tidyr::pivot_wider(names_from = indicator, values_from = value) |>
    dplyr::bind_rows(data.frame(name = "true", model = c("conditional", "marginal"), anchored = NA, estimate = c(true_conditional_AB, true_marginal_AB), variance = var_theo_AB))
  return(rectangle_results)
}

n_iter = 10

# pb <- progress::progress_bar$new(total = n_iter)
time_start <- Sys.time()
results_simulations <- lapply(1:n_iter, \(i) {
  # pb$tick()
  time_start_iteration <- Sys.time()
  cat("Iteration ", i, "\n")
  comparison(pop_init, struct_results, N_RCT, N_BOOT_ITER)
  time_eluded <- Sys.time() - time_start_iteration
  cat("Iteration lenth: ", time_eluded, "seconds\n")
})
cat("Simulation length ", time_start - Sys.time(), "seconds \n")

saveRDS(results_simulations, file = file.path("results_simulations", paste0("results_simulations", ".RDS")))
