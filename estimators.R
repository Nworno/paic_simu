#######################
## ESTIMATORS DEFINTION
#######################

###############################
########## Unadjusted estimator
##############################

# (Anchored) naive observed effect in pop_init
unadjusted_estimator <- function(trial_AC,
                                 trial_BC,
                                 anchored, 
                                 glm_family) {
  if (anchored) { # ie two steps
    naive_conditional_model_AC <- glm(Y_obs ~ ttt, family = glm_family, data = trial_AC)
    naive_conditional_model_BC <- glm(Y_obs ~ ttt, family = glm_family, data = trial_BC)
    naive_AB <- naive_conditional_model_AC$coefficients[["tttA"]] - naive_conditional_model_BC$coefficients[["tttB"]]
  } else {
    both_trials <- rbind(trial_AC[ttt == "A", .(Y_obs, ttt)], trial_BC[ttt == "B", .(Y_obs, ttt)])
    both_trials[ ,ttt := relevel(ttt, ref = "B")]
    naive_conditional_model_AB <- glm(Y_obs ~ ttt, family = glm_family, data = both_trials)
    naive_AB <- naive_conditional_model_AB$coefficients[["tttA"]]
  }
  return(naive_AB)
}

# (Unanchored) unadjusted observed effect in pop_init
run_unadjusted_estimator <- function(trial_AC, trial_BC, anchored, glm_family) {
  estimate <- unadjusted_estimator(trial_AC, trial_BC, anchored, glm_family)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) unadjusted_estimator(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                                    trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
                                                                    anchored, 
                                                                    glm_family))
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var(na.rm = TRUE)
  return(list("estimate" = estimate, "variance" = variance))
}


###########################
########## PROPENSITY SCORE
###########################

#### Estimating wegihts with maximum likelihood logistic regression (propensity score) 
max_likelihood <- function(df, model, dependent_variable) {
  PS_BC_trial <- glm(model, df, family = binomial(link = "logit"))$fitted.values
  ATC_w <- (df[[dependent_variable]] %in% "BC") + (df[[dependent_variable]] == "AC") * PS_BC_trial / (1 - PS_BC_trial)
  return(ATC_w)
}

### Estimating weights with method of moments 
# Defining method of moments functions
mm_obj_fun <- function(params, X) {
  sum(exp(X %*% params))
}
mm_grad_fun <- function(params, X) {
  colSums(sweep(X, 1, exp(X %*% params), FUN = "*"))
}
mm <- function(df, mean_covariates) {
  centered_IPD <- sweep(df, 2, mean_covariates, "-")
  random_init <- rep(0, ncol(df))
  
  params <- optim(random_init, fn = mm_obj_fun, gr = mm_grad_fun, method = "BFGS", X = df)$par
  weights <- exp(df %*% params) 
  return(weights)
}

propensity_score <- function(trial_AC,
                             trial_BC, 
                             covariate_names, 
                             anchored, 
                             weight_estimation_method,
                             studying_populations = FALSE) {
  if (anchored) {
    df_anchored <- data.table::rbindlist(list("AC" = trial_AC, "BC" = trial_BC), idcol = "trial", fill = TRUE)
    df_anchored[, ttt := relevel(as.factor(ttt), ref = "B")]
    stopifnot(levels(df_anchored$ttt)[[1]] == "B")
    if (weight_estimation_method == "max_likelihood") {
      df_anchored[, trial := relevel(as.factor(trial), ref = "AC")]
      stopifnot(levels(df_anchored$trial)[[1]] == "AC")
      trial_assigment_model <- paste0("trial ~ ", paste0(covariate_names, collapse = " + "))
      PS_BC_trial <- glm(trial_assigment_model, df_anchored, family = binomial(link = "logit"))$fitted.values
      AC_weights <- (df_anchored[["trial"]] == "BC") + (df_anchored[["trial"]] == "AC") * PS_BC_trial / (1 - PS_BC_trial)
    } else if (weight_estimation_method == "moments") {
      AC_weights <- mm(trial_AC[, covariate_names, with = FALSE] |> data.matrix(),
                    colMeans(trial_BC[, covariate_names, with = FALSE])) |>
        c(rep(1, nrow(trial_BC)))
    } else {
      stop("No weight estimation method provided")
    }
    model_marginal_AC <- glm(Y_obs ~ ttt, family = gaussian, data = df_anchored, weights = AC_weights) 
    estimate_AB <- model_marginal_AC$coefficients[["tttA"]]
  } else {
    trial_A <- trial_AC[ttt == "A"]
    trial_B <- trial_BC[ttt == "B"]
    if (weight_estimation_method == "max_likelihood") {
      df_unanchored <- data.table::rbindlist(list(trial_A, trial_B))
      df_unanchored[, ttt := relevel(as.factor(ttt), ref = "A")]
      stopifnot(levels(df_unanchored$ttt)[[1]] == "A")
      ttt_assigment_model <- paste0("ttt ~ ", paste0(covariate_names, collapse = " + "))
      PS_B_ttt <- glm(ttt_assigment_model, df_unanchored, family = binomial(link = "logit"))$fitted.values
      A_weights <- (PS_B_ttt / (1 - PS_B_ttt))[df_unanchored[["ttt"]] == "A"]
    } else if (weight_estimation_method == "moments") {
      A_weights <- mm(trial_A[, covariate_names, with = FALSE] |> data.matrix(), 
                      colMeans(trial_B[, covariate_names, with = FALSE]))
    } else {
      stop("No weighting method provided")
    }
    estimate_AB <- trial_A[, weighted.mean(Y_obs, A_weights)] - trial_B[, mean(Y_obs)]
  }
  return(estimate_AB)
}

run_propensity_score <- function(trial_AC, trial_BC, covariate_names, anchored, weight_estimation_method, studying_populations = FALSE) {
  estimate <- propensity_score(trial_AC,
                               trial_BC, 
                               covariate_names, 
                               anchored, 
                               weight_estimation_method, 
                               studying_populations)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) propensity_score(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                                trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                                covariate_names, 
                                                                anchored, 
                                                                weight_estimation_method, 
                                                                studying_populations))
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}

############################
########## REGRESSION MODELS
############################
regression_model <- function(trial_AC,
                             trial_BC, 
                             outcome_regression_model, 
                             covariate_names,
                             anchored, 
                             full_ipd,
                             glm_family) {
    if (full_ipd) {
      df_full_ipd <- data.table::rbindlist(list("AC" = trial_AC, "BC" = trial_BC),
                                           idcol = "trial",
                                           fill = TRUE)
      if (!anchored) df_full_ipd <- df_full_ipd[ttt %in% c("A", "B"), ]
      df_full_ipd[, ttt := relevel(factor(ttt), ref = "B")]
      stopifnot(levels(df_full_ipd$ttt)[[1]] == "B")
      fitted_model <- glm(outcome_regression_model, data = df_full_ipd, family = glm_family)
      AB_relative_effect <- fitted_model$coefficients[["tttA"]] 
    } else {
      centered_trial_AC <- copy(trial_AC)
      centered_trial_AC[, ttt := relevel(factor(ttt), ref = "C")]
      stopifnot(levels(centered_trial_AC$ttt)[[1]] == "C")
      if (!anchored) {
        centered_trial_AC <- copy(centered_trial_AC[ttt == "A", ])
        trial_BC <- trial_BC[ttt == "B", ]
        outcome_regression_model <- gsub("*ttt", "", outcome_regression_model, fixed = TRUE)
      }
      centered_trial_AC[, (covariate_names) := Map(function(x, y) x - y, .SD, colMeans(trial_BC[, ..covariate_names])), .SDcols = covariate_names]
      fitted_model <- glm(outcome_regression_model, data = centered_trial_AC, family = glm_family)
      AB_relative_effect <- ifelse(anchored, 
                                   fitted_model$coefficients[["tttA"]], 
                                   fitted_model$coefficients[["(Intercept)"]] - mean(trial_BC$Y_obs))
    }
  return(AB_relative_effect)
}

run_regression_model <- function(trial_AC, trial_BC, outcome_regression_model, covariate_names, full_ipd, anchored, glm_family = gaussian) {
  estimate <- regression_model(trial_AC, trial_BC, outcome_regression_model, covariate_names, anchored, full_ipd, glm_family)
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) {
    regression_model(trial_AC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt], 
                     trial_BC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
                     outcome_regression_model, 
                     covariate_names,
                     anchored, 
                     full_ipd,
                     glm_family)
  })
  variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = estimate, "variance" = variance))
}


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



#########
### MAIC
#########
# 
# maic <- function(trial_AC, trial_BC, covariate_names, anchored, outcome_family, studying_populations = FALSE) {
#   if (anchored) {
#     mean_trial_BC <- trial_BC[, lapply(.SD, mean), .SDcols = covariate_names] |> data.matrix()
#     centered_covariates <- sweep(trial_AC[, ..covariate_names] |> data.matrix(),
#                                  2,
#                                  mean_trial_BC,
#                                  FUN = "-")
#   } else {
#     mean_trial_B <- trial_BC[ttt == "B", lapply(.SD, mean), .SDcols = covariate_names] |> data.matrix()
#     centered_covariates <- sweep(trial_AC[ttt == "A", ..covariate_names] |> data.matrix(),
#                                  2,
#                                  mean_trial_B,
#                                  FUN = "-")
#   }
#   weights <- mm(centered_covariates)
#   # ess <- sum(weights)^2 / sum(weights^2)
#   if (studying_populations) {
#     return(weights)
#   }
#   if (anchored) {
#     trial_AC$weights <- weights
#     fitted_AC_w <- glm(Y_obs ~ ttt, data = trial_AC, family = outcome_family, weights = weights)
#     anchored_MAIC_AC <- fitted_AC_w$coefficients[["tttA"]]
#     # Equivalent to
#     # trial_AC[ttt == "A", weighted.mean(Y_obs, weights)] - trial_AC[ttt == "C", weighted.mean(Y_obs, weights)]
#     obs_marginal_BC <- trial_BC[ttt == "B", mean(Y_obs)] - trial_BC[ttt == "C", mean(Y_obs)]
#     estimate_AB <- anchored_MAIC_AC - obs_marginal_BC
#   } else {
#     estimate_AB <- trial_AC[ttt == "A", weighted.mean(Y_obs, weights)] - trial_BC[ttt == "B", mean(Y_obs)]
#   }
#   return(estimate_AB)
# }
# 
# run_maic <- function(trial_AC, trial_BC, covariate_names, anchored, outcome_family) {
#   stopifnot(levels(trial_AC$ttt)[[1]] == "C")
#   estimate <- maic(trial_AC, trial_BC, covariate_names, anchored, outcome_family)
#   boot_estimates <- lapply(1:N_BOOT_ITER, \(x) maic(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
#                                                     trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
#                                                     covariate_names,
#                                                     anchored,
#                                                     outcome_family))
#   variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
#   return(list("estimate" = estimate, "variance" = variance))
# }

# 
# anchored_conditional_estimation <- function(centered_trial_AC, centered_trial_BC, regression_model, glm_family) {
#   # two-steps individual patient data network meta analysis --> would be interesting to compare performance differences of this method
#   # as compared to random effect NMA
#   # potentially less biased, but systematically less precise as compared to "one-step unanchored" approach,
#   # because takes into account unobserved confounding between A and C, and B and C
#   # Equivalent to random effect meta analysis
#   browser()
#   model_AC <- glm(regression_model, glm_family, data = centered_trial_AC) # adjusted conditional effect
#   # BC shouldn't be adjusted in most clinical trials, so this is a more favorable situation than what is usually done
#   # in practice, because usually comparing a conditional effect to a marginal one 
#   model_BC <- glm(regression_model, glm_family, data = centered_trial_BC) # adjusted conditional effect
#   estimate_AC <- model_AC$coefficients[["tttA"]]
#   estimate_BC <- model_BC$coefficients[["tttB"]]
#   estimate_AB <- estimate_AC - estimate_BC
#   
#   ## Equivalent to
#   # model_AC <- glm(formula(Y_obs ~ X1*ttt), glm_family, trial_AC)
#   # estimate_AC <- model_AC$coefficients[["tttA"]] + model_AC$coefficients[["X1:tttA"]] * trial_BC[, mean(X1)]
#   # model_BC <- glm(formula(Y_obs ~ X1*ttt), glm_family, data = trial_BC)
#   # estimate_BC <- model_BC$coefficients[["tttB"]] + model_BC$coefficients[["X1:tttB"]] * trial_BC[, mean(X1)]
#   # estimate_AB <- estimate_AC - estimate_BC
#   return(estimate_AB)
# }

# run_anchored_conditional_estimation <- function(centered_trial_AC, trial_BC, outcome_regression_model, glm_family) {
#   estimate <- anchored_conditional_estimation(centered_trial_AC, trial_BC, outcome_regression_model, glm_family)
#   boot_estimates <- lapply(1:N_BOOT_ITER, \(x) {
#     anchored_conditional_estimation(
#       centered_trial_AC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
#       trial_BC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
#       outcome_regression_model, gaussian)
#   })
#   variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
#   return(list("estimate" = estimate, "variance" = variance))
# }

### Unanchored
#### Note: not really following either anchored or unanchored scheme, basically
#### unanchored analysis when a common comparator is available: make the assumption
#### that there are no unobserved imbalance in terms of prognostic factors or
#### TEM between trials. Will be more precise in such a case, but biased in case
#### of residual confounding. Not used currently
# botharms_unanchored_conditional_estimation <- function(centered_trial_AC, centered_trial_BC, regression_model, glm_family) {
#   both_trials <- rbind(centered_trial_AC, centered_trial_BC)
#   both_trials[, `:=`(ttt = factor(ttt, levels = c("B", "A", "C")))]
#   model_AB <- glm(regression_model, data = both_trials, family = glm_family) # conditional effect
#   estimate_AB <- model_AB$coefficients[["tttA"]]
# }
# 


# run_unanchored_conditional_estimation <- function(centered_trial_AC, centered_trial_BC, outcome_regression_model, glm_family) {
#   estimate <- unanchored_conditional_estimation(centered_trial_AC, centered_trial_BC, outcome_regression_model, glm_family)
#   boot_estimates <- lapply(1:N_BOOT_ITER, \(x) unanchored_conditional_estimation(
#     centered_trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
#     centered_trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
#     outcome_regression_model,
#     glm_family
#   ))
#   variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
#   return(list("estimate" = estimate, "variance" = variance))
# }




##########
###### STC
##########


# stc <- function(centered_trial_AC, trial_BC, anchored, outcome_regression_model, covariate_names, outcome_family = gaussian) {
#   if (anchored) {
#     average_observed_BC <- trial_BC[ttt == "B", mean(Y_obs)] - trial_BC[ttt == "C", mean(Y_obs)]
#     stc_model <- glm(outcome_regression_model,
#                      family = outcome_family,
#                      data = centered_trial_AC)
#     simulated_AC <- stc_model$coefficients[["tttA"]]
#     # Equivalent to
#     # stc_model <- glm(Y_obs ~ ttt*X1, family = gaussian, data = trial_AC)
#     # estimate_AC <- stc_model$coefficients[["tttA"]] + stc_model$coefficients[["tttA:X1"]] * Ag_trial_BC[, sum(mean_X1*n)/sum(n)]
#     estimate_AB <- simulated_AC - average_observed_BC
#   } else {
#     outcome_regression_model <- paste0("Y_obs ~ ", paste0(covariate_names, collapse = " + "))
#     average_observed_B <- trial_BC[ttt == "B", mean(Y_obs)]
#     regression_model(trial_AC[ttt == "A", ], outcome_regression_model, outcome_family, "tttA")
#     
#     stc_model <- glm(outcome_regression_model,
#                      family = outcome_family,
#                      data = centered_trial_AC[ttt == "A",])
#     simulated_A <- stc_model$coefficients[["(Intercept)"]]
#     # Equivalent to
#     # ttt_A <- trial_AC[ttt == "A"]
#     # stc_model <- glm(Y_obs ~ X1, family = gaussian, data = ttt_A)
#     # simulated_A <- stc_model$coefficients[["(Intercept)"]] + stc_model$coefficients[["X1"]] * trial_BC[ttt == "B", mean(X1)]
#     estimate_AB <- simulated_A - average_observed_B
#   }
#   return(estimate_AB)
# }

# run_stc <- function(centered_trial_AC, trial_BC, anchored, outcome_regression_model, covariate_names, outcome_family) {
#   estimate <- stc(centered_trial_AC, trial_BC, anchored, outcome_regression_model, covariate_names, outcome_family)
#   boot_estimates <- lapply(1:N_BOOT_ITER, \(x) stc(centered_trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
#                                                    trial_BC[sample(1:.N, .N, replace = TRUE), .SD, by = ttt],
#                                                    anchored,
#                                                    outcome_regression_model,
#                                                    covariate_names,
#                                                    outcome_family))
#   variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
#   return(list("estimate" = estimate, "variance" = variance))
# }
