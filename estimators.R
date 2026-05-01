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

run_unadjusted_estimator <- function(trial_AC, trial_BC, anchored, glm_family) {
  result <- unadjusted_estimator(trial_AC, trial_BC, anchored, glm_family)
  return(list("estimate" = result$estimate, "variance" = result$variance))
}


###########################
########## PROPENSITY SCORE
###########################

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

  tryCatch.W.E <- function(expr) {
    W <- NULL
    w.handler <- function(w) { # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
         warning = W)
  }

  if (anchored) {
    # Two-step approach: separate weighted AC model and unweighted BC model.
    # Avoids conditioning on `trial` in a logistic model, which would estimate
    # a conditional OR (non-collapsible) rather than the targeted marginal OR.
    df_AC <- df[trial == "AC"]
    df_BC <- df[trial == "BC"]
    df_AC[, ttt := relevel(as.factor(ttt), ref = "C")]
    df_BC[, ttt := relevel(as.factor(ttt), ref = "C")]
    if (outcome_family$family == "binomial") {
      df_AC$y_0 <- 1 - df_AC$Y_obs
      df_BC$y_0 <- 1 - df_BC$Y_obs
      outcome_model_step <- as.formula(cbind(Y_obs, y_0) ~ ttt)
    } else {
      outcome_model_step <- as.formula(Y_obs ~ ttt)
    }
    design_AC <- survey::svydesign(id = ~1, weights = ~trial_weights, data = df_AC)
    fitted_BC_glm <- tryCatch(
      glm(outcome_model_step, data = df_BC, family = outcome_family),
      error = function(e) NULL
    )
    result_step <- tryCatch.W.E({
      fitted_AC <- survey::svyglm(outcome_model_step, design = design_AC, family = outcome_family)
      coef_AC <- fitted_AC$coefficients[["tttA"]]
      coef_BC <- if (!is.null(fitted_BC_glm)) fitted_BC_glm$coefficients[["tttB"]] else stop("BC model failed")
      list(AB = coef_AC - coef_BC, AC = coef_AC)
    })
    step_succeeded <- !inherits(result_step$value, "error")
    estimate_AC <- if (step_succeeded) result_step$value$AC else NA_real_
    estimate_AB <- list(
      value   = if (step_succeeded) result_step$value$AB else result_step$value,
      warning = result_step$warning
    )
    variance_BC <- if (outcome_family$family == "binomial") {
      # log-odds scale: use model-based variance to match the estimate scale
      if (!is.null(fitted_BC_glm))
        tryCatch(vcov(fitted_BC_glm)["tttB", "tttB"], error = function(e) NA_real_)
      else NA_real_
    } else {
      df_BC[, .(var = var(Y_obs) / .N), by = ttt][, var] |> sum()
    }
    variance_AC <- NA_real_   # computed via bootstrap in run_propensity_score
  } else {
    estimate_AB <- tryCatch.W.E({
      design_glm <- survey::svydesign(id=~1, weights =~trial_weights, data = df)
      fitted_glm <- survey::svyglm(outcome_model, design = design_glm, family = outcome_family)
      fitted_glm$coefficients[["tttA"]]
    })
    variance_BC <- df[trial == "BC", .(var = var(Y_obs)/.N), by = ttt][, var] |> sum()
    variance_AC <- df[trial == "AC"][, .(var = weighted.var(Y_obs, trial_weights)/sum(trial_weights)), by = ttt][, var] |> sum()
    estimate_AC <- df[trial == "AC"][, .(mean = weighted.mean(Y_obs, trial_weights)), by = ttt][ttt == "A", mean]
  }
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
  boot_estimates <- lapply(1:N_BOOT_ITER, \(x) propensity_score(trial_AC[sample(1:.N, .N, replace = TRUE), .SD, by = c("ttt")],
                                                                trial_BC,
                                                                covariate_names,
                                                                assignment_model,
                                                                anchored,
                                                                weight_estimation_method,
                                                                outcome_family,
                                                                are_binary_covariates))
  estimate_error <- if (!is.null(estimate_AB$warning)) {
    error_estimate = estimate_AB$warning
  } else if ("error" %in% class(estimate_AB$value)) {
    error_estimate = estimate_AB$value
    estimate_AB$value <- NA # estimate_AB$value is a real number when no error, ie when warning or nothing
  } else {
    error_estimate = NA
  }

  boot_estimates_AC <- sapply(boot_estimates, \(x) x$estimate_AC)
  boot_errors <- sapply(boot_estimates[!is.finite(boot_estimates_AC)], \(x) x$df, simplify = FALSE)
  boot_warnings <- sapply(boot_estimates[!is.finite(boot_estimates_AC)], \(x) x$estimate_AC, simplify = FALSE)

  variance_AC <- tryCatch(Filter(\(x) is.finite(x), boot_estimates_AC) |> var(),
                          error = function(e) {print("no valid estimation for AC variance"); return(NA)})


  variance_AB <- variance_BC + variance_AC # works only on a linear scale, so estimate output has to remain linear
  list_ps_results <- list("estimate" = estimate_AB$value, "variance" = variance_AB, "error_estimate" = error_estimate)
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
    if (full_ipd) {
      df_full_ipd <- data.table::rbindlist(list(trial_AC, trial_BC),
                                           fill = TRUE,
                                           use.names = TRUE)
      if (anchored) {
        predictors_model <- paste0(predictors_model, " + trial")
      } else {
        df_full_ipd <- df_full_ipd[ttt %in% c("A", "B"), ]
      }
      df_full_ipd[, ttt := relevel(as.factor(ttt), ref = "B")]
      stopifnot(levels(df_full_ipd$ttt)[[1]] == "B")
      mean_covariates_BC <- colMeans(trial_BC[, ..covariate_names])
      if (outcome_family$family == "survival") {
        df_full_ipd_centered <- sweep(df_full_ipd[, ..covariate_names], 2, mean_covariates_BC, "-") |>
          cbind(df_full_ipd[, .(duree_rando_suivi_j60, Y_obs, ttt, trial)]) |> data.table::as.data.table()
        outcome_regression_model <- as.formula(paste0("Surv(duree_rando_suivi_j60, Y_obs) ~ ", predictors_model))
        fitted_model <- coxph(outcome_regression_model, data = df_full_ipd_centered)
      } else {
        df_full_ipd_centered <- sweep(df_full_ipd[, ..covariate_names], 2, mean_covariates_BC, "-") |>
          cbind(df_full_ipd[, .(Y_obs, ttt, trial)]) |> data.table::as.data.table()
        outcome_regression_model <- paste0("Y_obs ~ ", predictors_model)
        fitted_model <- glm(outcome_regression_model,
                          data = df_full_ipd_centered,
                          family = outcome_family)
      }
      estimate_AB <- fitted_model$coefficients[["tttA"]]
      variance_AB <- vcov(fitted_model)["tttA", "tttA"]

    } else {
      #### STC
      if (anchored) {
        trial_AC_to_center <- copy(trial_AC)
        trial_AC_to_center[, ttt := relevel(factor(ttt), ref = "C")]
        stopifnot(levels(trial_AC_to_center$ttt)[[1]] == "C")
        trial_BC[, ttt := relevel(factor(ttt), ref = "C")]
      } else {
        trial_AC_to_center <- trial_AC[ttt == "A", ]
        trial_BC <- trial_BC[ttt == "B", ]
        predictors_model <- gsub("*ttt", "", predictors_model, fixed = TRUE)
      }
      mean_covariates_BC <- colMeans(trial_BC[, ..covariate_names])
      if (outcome_family$family == "survival") {
      centered_trial_AC <- sweep(trial_AC_to_center[, ..covariate_names], 2, mean_covariates_BC, "-") |>
        cbind(trial_AC_to_center[, .(duree_rando_suivi_j60, Y_obs, ttt)])
      } else {
        centered_trial_AC <- sweep(trial_AC_to_center[, ..covariate_names], 2, mean_covariates_BC, "-") |>
          cbind(trial_AC_to_center[, .(Y_obs, ttt)])
      }
      dtf <- data.table::rbindlist(list(
        centered_trial_AC,
        trial_BC[, names(centered_trial_AC), with = FALSE]
      ))
      dtf[, ttt := relevel(factor(ttt), ref = "B")]
      if (outcome_family$family == "survival") {
        outcome_regression_model <- as.formula(paste0("Surv(duree_rando_suivi_j60, Y_obs) ~ ", predictors_model))
        fitted_model <- coxph(outcome_regression_model, data = dtf)
        estimate_AB <- fitted_model$coefficients[["tttA"]]
        variance_AB <- vcov(fitted_model)["tttA", "tttA"]
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
  return(list("estimate" = results$estimate_AB, "variance" = results$variance_AB))
}

