#######################
## ESTIMATORS DEFINTION
#######################

###############################
########## Unadjusted estimator
##############################

weighted.var <- function(x, w) {
  scaled_w <- w / sum(w)
  sum(scaled_w * (x - weighted.mean(x, scaled_w))^2)/(sum(scaled_w) - 1/sum(w))
}


# (Anchored) naive observed effect in pop_init
unadjusted_estimator <- function(trial_AC,
                                 trial_BC,
                                 anchored,
                                 glm_family) {
  if (anchored) { # ie two steps
    naive_conditional_model_AC <- glm(Y_obs ~ ttt, family = glm_family, data = trial_AC)
    naive_conditional_model_BC <- glm(Y_obs ~ ttt, family = glm_family, data = trial_BC)
    naive_AB <- naive_conditional_model_AC$coefficients[["tttA"]] - naive_conditional_model_BC$coefficients[["tttB"]]
    # vcov(naive_conditional_model_AC, complete = TRUE)["tttA", "tttA"] is equivalent to var(trial_AC[ttt == "A", Y_obs])/sum(trial_AC$ttt == "A") + var(trial_AC[ttt == "C", Y_obs])/sum(trial_AC$ttt == "C")
    var_AB <- vcov(naive_conditional_model_AC, complete = TRUE)["tttA", "tttA"] + vcov(naive_conditional_model_BC)["tttB", "tttB"]

  } else {
    both_trials <- rbind(trial_AC[ttt == "A", .(Y_obs, ttt)], trial_BC[ttt == "B", .(Y_obs, ttt)])
    both_trials[ ,ttt := relevel(ttt, ref = "B")]
    naive_conditional_model_AB <- glm(Y_obs ~ ttt, family = glm_family, data = both_trials)
    naive_AB <- naive_conditional_model_AB$coefficients[["tttA"]]
    var_AB <- vcov(naive_conditional_model_AB)["tttA", "tttA"]
  }
  return(list("estimate" = naive_AB, "variance" = var_AB))
}

# (Unanchored) unadjusted observed effect in pop_init
run_unadjusted_estimator <- function(trial_AC, trial_BC, anchored, glm_family) {
  result <- unadjusted_estimator(trial_AC, trial_BC, anchored, glm_family)
  return(list("estimate" = result$estimate, "variance" = result$variance))
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
mm <- function(df, mean_covariates, moments_2, var_covariates = NULL, are_binary_covariates = NULL) {
  centered_mean <- sweep(df, 2, mean_covariates, "-")

  if (moments_2) {
    if (is.null(are_binary_covariates)) stop("are_binary_covariates argument should be provided when balancing on second moments")
    centered_IPD2 <- sweep(df[, !are_binary_covariates, drop = FALSE]^2,
                           2,
                           mean_covariates[!are_binary_covariates] ^ 2 + var_covariates[!are_binary_covariates],
                           FUN = "-")
    centered_IPD <- cbind(centered_mean, centered_IPD2) |> data.matrix()
  } else {
    centered_IPD <- centered_mean
  }
  random_init <- rep(0, ncol(centered_IPD))

  params <- optim(random_init, fn = mm_obj_fun, gr = mm_grad_fun, method = "BFGS", X = centered_IPD)$par
  weights <- exp(centered_IPD %*% params)
  return(weights)
}

propensity_score <- function(trial_AC,
                             trial_BC,
                             covariate_names,
                             assignment_model,
                             anchored,
                             weight_estimation_method = c("max_likelihood", "moments_1", "moments_2"),
                             outcome_family,
                             are_binary_covariates = NULL) {
  if (anchored) {
    ######### Anchored
    df <- data.table::rbindlist(list(trial_AC, trial_BC), fill = TRUE, use.names = TRUE)
    df[, ttt := relevel(as.factor(ttt), ref = "B")]
    stopifnot(levels(df$ttt)[[1]] == "B")
    if (weight_estimation_method == "max_likelihood") {
      #### Propensity score weights
      df[, trial := relevel(as.factor(trial), ref = "AC")]
      stopifnot(levels(df$trial)[[1]] == "AC")
      trial_assignment_model <- paste0("trial ~ ", assignment_model)
      PS_BC_trial <- glm(trial_assignment_model, df, family = binomial(link = "logit"))$fitted.values
      trial_weights <- (df[["trial"]] == "BC") + (df[["trial"]] == "AC") * PS_BC_trial / (1 - PS_BC_trial)
    } else if (weight_estimation_method == "moments_1") {
      ##### MAIC weights
      trial_weights <- mm(data.matrix(trial_AC[, covariate_names, with = FALSE]),
                          colMeans(trial_BC[, covariate_names, with = FALSE]),
                          moments_2 = FALSE) |>
        c(rep(1, nrow(trial_BC)))
    } else if (weight_estimation_method == "moments_2") {
      trial_BC[, lapply(.SD, var), .SDcols = covariate_names]
      trial_weights <- mm(data.matrix(trial_AC[, covariate_names, with = FALSE]),
                          colMeans(trial_BC[, covariate_names, with = FALSE]),
                          moments_2 = TRUE,
                          var_covariates = trial_BC[, lapply(.SD, var), .SDcols = covariate_names] |> data.matrix(),
                          are_binary_covariates = are_binary_covariates
      ) |>
        c(rep(1, nrow(trial_BC)))
    } else {
      stop("No weight estimation method provided")
    }
    if (outcome_family$family == "binomial") {
      df$y_0 <- 1 - df$Y_obs
      outcome_model <- as.formula(cbind(Y_obs, y_0) ~ ttt + trial)
    } else {
      outcome_model <- as.formula(Y_obs ~ ttt + trial)
    }
  } else {
    ######## Non-anchored
    trial_A <- trial_AC[ttt == "A"]
    trial_B <- trial_BC[ttt == "B"]
    df <- data.table::rbindlist(list(trial_A, trial_B), use.names = TRUE)
    if (weight_estimation_method == "max_likelihood") {
      #### Propensity score weights
      df[, ttt := relevel(as.factor(ttt), ref = "A")]
      stopifnot(levels(df$ttt)[[1]] == "A")
      ttt_assignment_model <- paste0("ttt ~ ", assignment_model)
      PS_B_ttt <- glm(ttt_assignment_model, df, family = binomial(link = "logit"))$fitted.values
      trial_weights <- (PS_B_ttt / (1 - PS_B_ttt)) * (df[["ttt"]] == "A") + (df[["ttt"]] == "B")
    } else if (weight_estimation_method == "moments_1") {
      ##### MAIC weights
      trial_weights <- mm(data.matrix(trial_A[, covariate_names, with = FALSE]),
                          colMeans(trial_B[, covariate_names, with = FALSE]),
                          moments_2 = FALSE) |>
        c(rep(1, nrow(trial_B)))
    } else if (weight_estimation_method == "moments_2") {
      trial_weights <- mm(data.matrix(trial_A[, covariate_names, with = FALSE]),
                          colMeans(trial_B[, covariate_names, with = FALSE]),
                          moments_2 = TRUE,
                          var_covariates = trial_B[, lapply(.SD, var), .SDcols = covariate_names] |> data.matrix(),
                          are_binary_covariates = are_binary_covariates
      ) |>
        c(rep(1, nrow(trial_B)))
      if (any(is.infinite(trial_weights))) print("ouuuh Infinite weights")
    } else {
      stop("No weighting method provided")
    }
    df[, ttt := relevel(as.factor(ttt), ref = "B")]
    stopifnot(levels(df$ttt)[[1]] == "B")
    if (outcome_family$family == "binomial") {
      df$y_0 <- 1 - df$Y_obs
      outcome_model <- as.formula(cbind(Y_obs, y_0) ~ ttt)
    } else {
      outcome_model <- as.formula(Y_obs ~ ttt)
    }
  }
    df[, trial_weights := trial_weights]

  # ALl estimated in one step, but similar as doing it in two steps when not adjusted on any confoundings. The "anchored" comparison is performed by the " + trial" in the model
  # fitted_glm <- glm(outcome_model,
  #                   family = outcome_family,
  #                   data = df,
  #                   weights = trial_weights)
  estimate_AB <- tryCatch({
    design_glm <- survey::svydesign(id=~1, weights =~trial_weights, data = df)
    fitted_glm <- survey::svyglm(outcome_model, design = design_glm, family = outcome_family)
    fitted_glm$coefficients[["tttA"]]
  },
  warning = \(w) w,
  error = \(e) e
  )

  # variance

  # These two variances effectively only count one arm in the case of unanchored comparisons
  variance_BC <- df[trial == "BC", .(var = var(Y_obs)/.N), by = ttt][, var] |> sum()
  #TODO: reprendre d'ici, récupérer erreur quand problème avec la variance pondérée : peut être récupérer simplement le nombre de NA
  variance_AC <- df[trial == "AC"][, .(var = weighted.var(Y_obs, trial_weights)/sum(trial_weights)), by = ttt][, var] |> sum()
  # Estimate AC
  estimate_A_and_C <- df[trial == "AC"][, .(mean = weighted.mean(Y_obs, trial_weights)), by = ttt]
  if (anchored) estimate_AC <- estimate_A_and_C[ttt == "A", mean] - estimate_A_and_C[ttt == "C", mean] else estimate_AC <- estimate_A_and_C[ttt == "A", mean]
  return(list(estimate_AB = estimate_AB, estimate_AC = estimate_AC, variance_AC = variance_AC, variance_BC = variance_BC, df = df))
}

run_propensity_score <- function(trial_AC, trial_BC, covariate_names, assignment_model, anchored, weight_estimation_method, outcome_family, retrieve_ps_weights, are_binary_covariates = NULL) {
  results <- propensity_score(trial_AC,
                              trial_BC,
                              covariate_names,
                              assignment_model,
                              anchored,
                              weight_estimation_method,
                              outcome_family,
                              are_binary_covariates)
  estimate_AB <- results$estimate_AB
  variance_BC <- results$variance_BC
  # variance_AC <- results$variance_AC
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) propensity_score(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                                trial_BC,
                                                                covariate_names,
                                                                assignment_model,
                                                                anchored,
                                                                weight_estimation_method,
                                                                outcome_family,
                                                                are_binary_covariates))
  boot_estimates_AC <- sapply(boot_estimates, \(x) x$estimate_AC)
  boot_errors <- sapply(boot_estimates[!is.finite(boot_estimates_AC)], \(x) x$df, simplify = FALSE)
  boot_warnings <- sapply(boot_estimates[!is.finite(boot_estimates_AC)], \(x) x$estimate_AC, simplify = FALSE)

  # variance_AC <- Filter(is.numeric, boot_estimates_AC) |> unlist() |> var()
  variance_AC <- tryCatch(Filter(\(x) is.finite(x), boot_estimates_AC) |> var(),
                          error = function(e) {print("no valid estimation for AC variance"); return(NA)})

  # Should be close to consistent with bootstrap
  # boot_estimates_var_AC <- lapply(1:1000, \(x) propensity_score(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
  #                                                           trial_BC,
  #                                                           covariate_names,
  #                                                           assignment_model,
  #                                                           anchored,
  #                                                           weight_estimation_method,
  #                                                           outcome_family,
  #                                                           studying_populations)$variance_AC)
  # mean_variance_AC <- Filter(is.numeric, boot_estimates_var_AC) |> unlist() |> mean()

  variance_AB <- variance_BC + variance_AC # works only on a linear scale, so estimate output has to remain linear
  list_ps_results <- list("estimate" = estimate_AB, "variance" = variance_AB)
  if (length(boot_errors) > 0) list_ps_results$boot_errors <- boot_errors
  if (length(boot_errors) > 0) list_ps_results$boot_warnings <- boot_warnings
  if (retrieve_ps_weights) list_ps_results$df <- results$df
  return(list_ps_results)
}

############################
########## REGRESSION MODELS
############################
regression_model <- function(trial_AC,
                             trial_BC,
                             predictors_model,
                             covariate_names,
                             anchored,
                             full_ipd,
                             outcome_family) {
  # if (all(c("X1", "X2") %in% covariate_names) & bY_A_X1 == 1 & bY_A_X2 == 1 & ) browser()

    if (full_ipd) {
      # Classic IPD -->
      ### Classic regression model
      df_full_ipd <- data.table::rbindlist(list(trial_AC, trial_BC),
                                           fill = TRUE,
                                           use.names = TRUE)
      if (anchored) {
        # 2 differences between anchored/unanchored:
        # "two times" more data to estimate predictors' effect + 'trial' variable in the model
        predictors_model <- paste0(predictors_model, " + trial")
      } else {
        df_full_ipd <- df_full_ipd[ttt %in% c("A", "B"), ]
      }
      df_full_ipd[, ttt := relevel(as.factor(ttt), ref = "B")]
      stopifnot(levels(df_full_ipd$ttt)[[1]] == "B")
      mean_covariates_BC <- colMeans(trial_BC[, ..covariate_names])
      df_full_ipd_centered <- sweep(df_full_ipd[, ..covariate_names], 2, mean_covariates_BC, "-") |>
        cbind(df_full_ipd[, .(duree_rando_suivi_j60, Y_obs, ttt, trial)]) |> data.table::as.data.table()
      if (outcome_family == "survival") {
        outcome_regression_model <- as.formula(paste0("Surv(duree_rando_suivi_j60, Y_obs) ~ ", predictors_model))
        fitted_model <- coxph(outcome_regression_model, data = df_full_ipd_centered)
      } else {
        outcome_regression_model <- paste0("Y_obs ~ ", predictors_model)
        fitted_model <- glm(outcome_regression_model,
                          data = df_full_ipd_centered,
                          family = outcome_family)
      }
      estimate_AB <- fitted_model$coefficients[["tttA"]]
      variance_AB <- vcov(fitted_model)["tttA", "tttA"]

      # summary(fitted_model)
      # df_full_ipd_centered[, var(Y_obs)/.N, by = c("trial", "ttt")][, V1] |> sum()

    } else {
      #### STC
      if (anchored) {
        trial_AC_to_center <- copy(trial_AC)
        trial_AC_to_center[, ttt := relevel(factor(ttt), ref = "C")]
        stopifnot(levels(trial_AC_to_center$ttt)[[1]] == "C")
        trial_BC[, ttt := relevel(factor(ttt), ref = "C")]
        predictors_model <- paste0(predictors_model, " + trial")
      } else {
        trial_AC_to_center <- trial_AC[ttt == "A", ]
        trial_BC <- trial_BC[ttt == "B", ]
        predictors_model <- gsub("*ttt", "", predictors_model, fixed = TRUE)
      }
      mean_covariates_BC <- colMeans(trial_BC[, ..covariate_names])
      centered_trial_AC <- sweep(trial_AC_to_center[, ..covariate_names], 2, mean_covariates_BC, "-") |>
        cbind(trial_AC_to_center[, .(duree_rando_suivi_j60, Y_obs, ttt)])

      dtf <- data.table::rbindlist(list(
        centered_trial_AC,
        trial_BC
      ))
      dtf[, ttt := relevel(factor(ttt), ref = "B")]
      if (outcome_family == "survival") {
        outcome_regression_model <- as.formula(paste0("Surv(duree_rando_suivi_j60, Y_obs) ~ ", predictors_model))
        fitted_model <- coxph(outcome_regression_model, data = dtf)
        estimate_AB <- fitted_model$coefficients[["tttA"]]
        variance_AB <- vcov(fitted_model)["tttA", "tttA"]
        # if (anchored) {
        #   fitted_model_BC <- coxph(Surv(duree_rando_suivi_j60, Y_obs) ~ ttt, data = trial_BC)
        # } else {
        #   fitted_model_BC <- coxph(Surv(duree_rando_suivi_j60, Y_obs) ~ 1, data = trial_BC)
        # }
      } else {
        outcome_regression_model <- paste0("Y_obs ~ ", predictors_model)
        fitted_model_AC <- glm(outcome_regression_model,
                               data = centered_trial_AC,
                               family = outcome_family)
        if (anchored) {
          fitted_model_BC <- glm(Y_obs ~ ttt, data = trial_BC, family = outcome_family)
        } else {
          fitted_model_B <- glm(Y_obs ~ 1, data = trial_BC, family = outcome_family)

        }
      }
      if (anchored) {
        estimate_BC <- fitted_model_BC$coefficients[["tttB"]]
        estimate_AC <- fitted_model_AC$coefficients[["tttA"]]
        estimate_AB <- estimate_AC - estimate_BC

        variance_AC <- vcov(fitted_model_AC)["tttA", "tttA"]
        variance_BC <- vcov(fitted_model_BC)["tttB", "tttB"]
        variance_AB <- variance_AC + variance_BC
      } else {
        estimate_B <- fitted_model_B$coefficients[["(Intercept)"]]
        estimate_A <- fitted_model_AC$coefficients[["(Intercept)"]]
        estimate_AB <- estimate_A - estimate_B

        variance_AC <- vcov(fitted_model_AC)["(Intercept)", "(Intercept)"]
        variance_BC <- vcov(fitted_model_B)["(Intercept)", "(Intercept)"]
        variance_AB <- variance_AC + variance_BC
      }
    }
  return(list(estimate_AB = estimate_AB, variance_AB = variance_AB))
}

run_regression_model <- function(trial_AC, trial_BC, outcome_regression_model, covariate_names, full_ipd, anchored, outcome_family) {
  results <- regression_model(trial_AC, trial_BC, outcome_regression_model, covariate_names, anchored, full_ipd, outcome_family)
  # boot_estimates <- lapply(1:N_BOOT_ITER, \(x) {
  #   regression_model(trial_AC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
  #                    trial_BC[sample(1:.N, size = .N, replace = TRUE), .SD, by = ttt],
  #                    outcome_regression_model,
  #                    covariate_names,
  #                    anchored,
  #                    full_ipd,
  #                    outcome_family)
  # })
  # variance <- Filter(is.numeric, boot_estimates) |> unlist() |> var()
  return(list("estimate" = results$estimate_AB, "variance" = results$variance_AB))
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
