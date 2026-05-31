library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
source("env_variables.R")

dir_experience_results <- file.path("results_simulations", DATE_EXPERIMENT)

list_files <- {
  dirs <- list.dirs(dir_experience_results, full.names = TRUE)
  res <- lapply(dirs, \(x) list.files(x, full.names = TRUE, pattern = "^experiment_results.*\\.RDS"))
  names(res) <- dirs
  Filter(f = \(x) length(x) != 0, res)
}

nested_list_results_df <- rapply(list_files, classes = "character", how = "replace", \(x) {
  sapply(x, readRDS, simplify = FALSE)
})
long_df_results <- lapply(nested_list_results_df, \(l) {
  lapply(l, rbindlist, idcol = "iteration", use.names = TRUE)
  }) |>
  lapply(rbindlist, idcol = "estimator_num", use.names = TRUE) |>
  lapply(\(df) {suppressWarnings(df$error_estimate <- NULL); return(df)}) |> # removing the column problems, to use it separately
  rbindlist(idcol = "population_parameters_num", use.names = TRUE)
long_df_results[, population_parameters_num := basename(population_parameters_num) |> as.integer()]
long_df_results[, estimator_num := gsub(pattern = ".*(?<=experiment_results_)(\\d+)(?=\\.RDS).*",
                                        replacement = "\\1", x = estimator_num, perl = TRUE)]

long_df_problems <- lapply(nested_list_results_df, \(l) {
  lapply(l, rbindlist, idcol = "iteration", use.names = TRUE)
}) |>
  lapply(rbindlist, idcol = "estimator_num", use.names = TRUE) |>
  lapply(\(df) {df$estimate <- NULL; df$variance <- NULL; return(df)}) |>
  rbindlist(idcol = "population_parameters_num", use.names = TRUE, fill = TRUE) |>
  dplyr::mutate(estimator_num = sub(estimator_num, pattern = ".*experiment_results_", replacement = "", perl = TRUE)) |>
  dplyr::mutate(population_parameters_num = basename(population_parameters_num) |> as.integer(),
                estimator_num = sub(estimator_num, pattern = "(?<=[0-9])\\.RDS", replacement = "", perl = TRUE) |> as.integer(),
                adjustment = case_match(adjustment,
                                        "unadjusted" ~ "Unadjusted",
                                        "regression" ~ "Regression",
                                        "iptw" ~ "IPTW"),
                model = case_match(model,
                                   "unadjusted" ~ "Unadjusted",
                                   "glm" ~ "GLM",
                                   "stc" ~ "STC",
                                   "ml" ~ "ML",
                                   "maic_1" ~ "MAIC_1",
                                   "maic_2" ~ "MAIC_2",
                                   .default = model),
                anchored = case_match(anchored,
                                      "anchored" ~ "Anchored",
                                      "unanchored" ~ "Unanchored",
                                      .default = anchored)) |>
  dplyr::rowwise() |>
  dplyr::mutate(error_estimate = purrr::pluck(error_estimate, "message", .default = NA)) |>
  dplyr::ungroup() |>
  dplyr::rename_with(.fn = stringr::str_to_title)


df_true_effects <- list.dirs(dir_experience_results) |>
  lapply(list.files, pattern = "average_outcome_df", full.names = TRUE) |>
  Filter(f = \(x) length(x) != 0) |>
  lapply(\(x) {
    name <- stringr::str_extract(pattern = "\\d+(?=\\/average_outcome_df\\.RDS)", x)
    value <- readRDS(x)
    alist <- list(value)
    names(alist) <- name
    return(alist)
    }) |>
  unlist(recursive = FALSE) |>
  rbindlist(idcol = "population_parameters_num", use.names = TRUE)

df_population_parameters <- readRDS(file = file.path(dir_experience_results, "df_population_parameters.RDS")) |>
  dplyr::mutate(across(where(is.list), as.character))
df_estimators_parameters <- readRDS(file = file.path(dir_experience_results, "df_estimators_parameters.RDS"))

long_df_results <- long_df_results |>
  dplyr::mutate(data = case_match(model,
                           c("maic_1", "maic_2", "stc") ~ "PAIC",
                           c("ml", "glm") ~ "IPD",
                           "unadjusted" ~ "AgD") |> factor(levels = c("AgD", "PAIC", "IPD")),
         outcome_type = case_match(adjustment,
                                   c("unadjusted", "iptw") ~ "marginal",
                                   "regression" ~ "conditional"),
         anchored = factor(anchored, levels = c("unanchored", "anchored")) |> stringr::str_to_title(),
         adjustment = ifelse(adjustment == "iptw", "IPTW", stringr::str_to_title(adjustment)),
         model = ifelse(model == "unadjusted", "Unadjusted", toupper(model)))

joined_results <- df_true_effects[trial == "BC", .(population_parameters_num = as.integer(population_parameters_num), outcome_type, true_effect = AB)][long_df_results, , on = c("population_parameters_num", "outcome_type")] |>
  dplyr::mutate(estimator_num = as.integer(estimator_num)) |>
  dplyr::rename_with(stringr::str_to_title) |>
  tibble::as_tibble() |>
  dplyr::inner_join(long_df_problems, by = c("Model", "Population_parameters_num", "Anchored", "Estimator_num", "Iteration", "Adjustment"))

get_bias <- function(obs, theo) {
  mean(obs - theo, na.rm = TRUE)
}
get_RMSE <- function(obs, theo) {
  sqrt(mean((obs - theo)**2, na.rm = TRUE))
}
get_VR <- function(obs, se_obs) {
  mean(se_obs, na.rm = TRUE) / sd(obs, na.rm = TRUE)
}
get_cov_95 <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975)*se
  lb <- coef - qnorm(0.975)*se
  covered <- theo < ub & theo > lb
  mean(covered, na.rm = TRUE)
}

get_number_na_estimate <- function(Estimate) sum(is.na(Estimate))

get_se_bias <- function(obs, theo) {
  n <- sum(!is.na(obs - theo))
  sd(obs - theo, na.rm = TRUE) / sqrt(n)
}

get_se_RMSE <- function(obs, theo) {
  sq_err <- (obs - theo)^2
  rmse_val <- sqrt(mean(sq_err, na.rm = TRUE))
  n <- sum(!is.na(sq_err))
  sqrt(var(sq_err, na.rm = TRUE) / n) / (2 * rmse_val)
}

get_se_VR <- function(obs, se_obs) {
  n <- sum(!is.na(obs) & !is.na(se_obs))
  mean_se <- mean(se_obs, na.rm = TRUE)
  sd_est  <- sd(obs, na.rm = TRUE)
  vr <- mean_se / sd_est
  # Delta method: Var(VR) ≈ VR² * (Var(se)/n/mean_se² + 1/(2*(n-1)))
  vr * sqrt(var(se_obs, na.rm = TRUE) / (n * mean_se^2) + 1 / (2 * (n - 1)))
}

get_se_cov_95 <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975) * se
  lb <- coef - qnorm(0.975) * se
  covered <- (theo < ub & theo > lb)
  p <- mean(covered, na.rm = TRUE)
  n <- sum(!is.na(covered))
  sqrt(p * (1 - p) / n)
}

correct_decision <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975)*se
  lb <- coef - qnorm(0.975)*se
  contains_0 <- ub >= 0 & lb <= 0
  sign_estimate_theo_identical <- (theo >= 0 & coef >= 0) | (theo <= 0 & coef <= 0)
  theo_is_null <- theo == 0
  correct_decision <- ifelse(theo_is_null, contains_0, sign_estimate_theo_identical)
  mean(correct_decision, na.rm = TRUE)
}

df_stats_wide <- joined_results |>
  group_by(Population_parameters_num, Estimator_num, Adjustment, Model, Anchored, Data) |>
  mutate(Estimate = ifelse(!is.na(Error_estimate), NA, Estimate),
         Variance = ifelse(!is.na(Error_estimate), NA, Variance)) |>
  summarize(bias = get_bias(Estimate, True_effect),
            se_bias = get_se_bias(Estimate, True_effect),
            rmse = get_RMSE(Estimate, True_effect),
            se_rmse = get_se_RMSE(Estimate, True_effect),
            vr = get_VR(Estimate, sqrt(Variance)),
            se_vr = get_se_VR(Estimate, sqrt(Variance)),
            cov_95 = get_cov_95(Estimate, sqrt(Variance), True_effect),
            se_cov_95 = get_se_cov_95(Estimate, sqrt(Variance), True_effect),
            number_na_estimate = get_number_na_estimate(Estimate),
            correct_decision = correct_decision(Estimate, Variance, True_effect),
            .groups = "drop")

group_cols <- c("Population_parameters_num", "Estimator_num", "Adjustment", "Model", "Anchored", "Data")

df_stats <- df_stats_wide |>
  select(-starts_with("se_")) |>
  pivot_longer(cols = c("bias", "rmse", "vr", "cov_95", "correct_decision", "number_na_estimate"),
               names_to = "indicator", values_to = "values") |>
  dplyr::left_join(
    df_stats_wide |>
      pivot_longer(cols = c("se_bias", "se_rmse", "se_vr", "se_cov_95"),
                   names_to = "indicator", names_prefix = "se_", values_to = "se_mc") |>
      dplyr::select(all_of(c(group_cols, "indicator", "se_mc"))),
    by = c(group_cols, "indicator")
  )

dir.create(file.path(dir_experience_results, "processed_results"))
saveRDS(df_stats, file.path(dir_experience_results, "processed_results", "df_stats.rds"))
saveRDS(joined_results, file.path(dir_experience_results, "processed_results", "joined_results.rds"))
saveRDS(long_df_problems, file.path(dir_experience_results, "processed_results", "long_df_problems.rds"))
