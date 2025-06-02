library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
source("env_variables.R")

dir_experience_results <- file.path("results_simulations", DATE_EXPERIMENT)

list_files <- lapply(list.dirs(dir_experience_results, full.names = TRUE), \(x) {
  list.files(x, full.names = TRUE, pattern = "^experiment_results.*\\.RDS")
}) |> Filter(f = \(x) length(x) != 0)

nested_list_results_df <- rapply(list_files, classes = "character", how = "replace", \(x) {
  sapply(x, readRDS, simplify = FALSE)
})
long_df_results <- lapply(nested_list_results_df, \(l) {
  lapply(l, rbindlist, idcol = "iteration", use.names = TRUE)
  }) |>
  lapply(rbindlist, idcol = "estimator_num", use.names = TRUE) |>
  lapply(\(df) {suppressWarnings(df$error_estimate <- NULL); return(df)}) |> # removing the column problems, to use it separately
  rbindlist(idcol = "population_parameters_num", use.names = TRUE)
long_df_results[, estimator_num := gsub(pattern = ".*(?<=experiment_results_)(\\d+)(?=\\.RDS).*",
                                        replacement = "\\1", x = estimator_num, perl = TRUE)]

long_df_problems <- lapply(nested_list_results_df, \(l) {
  lapply(l, rbindlist, idcol = "iteration", use.names = TRUE)
}) |>
  lapply(rbindlist, idcol = "estimator_num", use.names = TRUE) |>
  lapply(\(df) {df$estimate <- NULL; df$variance <- NULL; return(df)}) |>
  rbindlist(idcol = "population_parameters_num", use.names = TRUE, fill = TRUE) |>
  dplyr::mutate(estimator_num = sub(estimator_num, pattern = ".*experiment_results_", replacement = "", perl = TRUE)) |>
  dplyr::mutate(estimator_num = sub(estimator_num, pattern = "(?<=[0-9])\\.RDS", replacement = "", perl = TRUE) |> as.integer(),
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

combined_parameters <- merge(df_population_parameters[, .((.SD), key = 1)],
                             df_estimators_parameters[, .((.SD), key = 1)],
                             all.x = TRUE,
                             allow.cartesian = TRUE)

df_true_effects_merged <- df_true_effects[, population_parameters_num := as.integer(population_parameters_num)] |>
  merge(combined_parameters,
        all.x = TRUE,
        allow.cartesian = TRUE)

########################
# Aligning DGM and estimators models:
# finding out which scenarios are biased and which are not
########################


true_PF <- df_population_parameters |>
  select(population_parameters_num, matches("bY_X")) |>
  mutate(across(matches("bY_X"), \(x) x != 0)) |>
  rename_with(\(x) sub("bY_", "", x)) |>
  pivot_longer(-population_parameters_num, names_to = "PF") |>
  left_join(combined_parameters[, c("population_parameters_num", "estimator_num")],
            relationship = "many-to-many")

true_TEM <- df_population_parameters |>
  select(population_parameters_num, matches("bY_A_X")) |>
  mutate(across(matches("bY_A_X"), \(x) x != 0)) |>
  rename_with(\(x) sub("bY_A_", "", x)) |>
  pivot_longer(-population_parameters_num, names_to = "TEM") |>
  left_join(combined_parameters[, c("population_parameters_num", "estimator_num")],
            relationship = "many-to-many")

marginal_model <- df_estimators_parameters[, .(estimator_num, covariate_names)] |>
  unnest_longer(covariate_names, values_to = "variable") |>
  mutate(present = TRUE)

conditional_model <- df_estimators_parameters[, .(estimator_num,
                                                  TEM = stringr::str_extract_all(outcome_regression_model, "X\\d+(?=\\*)"),
                                                  PF = stringr::str_extract_all(outcome_regression_model, "X\\d+"))] |>
  mutate(present = TRUE)

prognostics_conditional <- unnest_longer(conditional_model[,c("estimator_num", "PF", "present")],
                                         col = "PF", values_to = "variable")
TEM_conditional <- unnest_longer(conditional_model[,c("estimator_num", "TEM", "present")], col = "TEM", values_to = "variable")

confounded_prognostics <- left_join(true_PF,
                                    prognostics_conditional,
                                    by = c("estimator_num", "PF" = "variable"),
                                    relationship = "many-to-many") |>
  mutate(present = replace_na(present, FALSE)) |>
  mutate(adequate = value == present,
         confounded = value & !present,
         imprecise = !value & present) |>
  group_by(population_parameters_num, estimator_num) |>
  summarize(adequate = all(adequate),
            confounded = any(confounded),
            imprecise = any(imprecise),
            confounder = "prognostic") |>
  ungroup() |>
  rowwise() |>
  mutate(string = list(c("adequate", "confounded", "imprecise")[c(adequate, confounded, imprecise)]))

confounded_tem <- left_join(true_TEM,
                            TEM_conditional,
                            by = c("estimator_num", "TEM" = "variable"),
                            relationship = "many-to-many") |>
  mutate(present = replace_na(present, FALSE)) |>
  mutate(adequate = value == present,
         confounded = value & !present,
         imprecise = !value & present) |>
  group_by(population_parameters_num, estimator_num) |>
  summarize(adequate = all(adequate),
            confounded = any(confounded),
            imprecise = any(imprecise),
            confounder = "tem") |>
  ungroup() |>
  rowwise() |>
  mutate(string = list(c("adequate", "confounded", "imprecise")[c(adequate, confounded, imprecise)]))

confounded_conditional <- bind_rows(confounded_prognostics,
                                    confounded_tem) |>
  group_by(population_parameters_num, estimator_num) |>
  summarize(adequate = all(adequate),
            confounded = any(confounded),
            imprecise = any(imprecise)) |>
  ungroup() |>
  rowwise() |>
  mutate(model_conditional = list(c("adequate", "confounded", "imprecise")[c(adequate, confounded, imprecise)]))

confounded_marginal <- full_join(true_TEM,
                                 true_PF,
                                 by = c("population_parameters_num", "estimator_num", "TEM" = "PF")) |>
  mutate(value = value.x | value.y) |>
  mutate(value = replace_na(value, FALSE)) |>
  select(-value.x, -value.y) |>
  left_join(marginal_model,
            by = c("estimator_num", "TEM" = "variable")) |>
  mutate(present = replace_na(present, FALSE)) |>
  mutate(adequate = value == present,
         confounded = value & !present,
         imprecise = !value & present) |>
  group_by(population_parameters_num, estimator_num) |>
  summarize(adequate = all(adequate),
            confounded = any(confounded),
            imprecise = any(imprecise)) |>
  ungroup() |>
  rowwise() |>
  mutate(model_marginal = list(c("adequate", "confounded", "imprecise")[c(adequate, confounded, imprecise)]))

nonnormal_distribution <- df_population_parameters |>
  select(c(population_parameters_num, matches("f_X"))) |>
  # TODO: pattern à mettre à jour pour récupérer toutes les façons de créer
  # des variables distribuées non normalement
  mutate(across(matches("f_X"), \(x) grepl("sample", x))) |>
  rowwise() |>
  mutate(any_nonnormal = any(c_across(f_X1:f_X4)),
         nonnormal_distribution = list(c("f_X1", "f_X2", "f_X3", "f_X4")[c(f_X1, f_X2, f_X3, f_X4)]))


long_df_results <- long_df_results |>
  mutate(data = case_match(model,
                           c("maic_1", "maic_2", "stc") ~ "PAIC",
                           c("ml", "glm") ~ "IPD",
                           "unadjusted" ~ "AgD") |> factor(levels = c("AgD", "PAIC", "IPD")),
         outcome_type = case_match(adjustment,
                                   c("unadjusted", "iptw") ~ "marginal",
                                   "regression" ~ "conditional"),
         anchored = factor(anchored, levels = c("unanchored", "anchored")) |> stringr::str_to_title(),
         adjustment = ifelse(adjustment == "iptw", "IPTW", stringr::str_to_title(adjustment)),
         model = ifelse(model == "unadjusted", "Unadjusted", toupper(model)))

joined_results <- df_true_effects[trial == "BC", .(population_parameters_num, outcome_type, true_effect = AB)][long_df_results, , on = c("population_parameters_num", "outcome_type")] |>
  mutate(across(c(population_parameters_num, estimator_num), as.integer)) |>
  rename_with(stringr::str_to_title) |>
  tibble::as_tibble() |>
  dplyr::inner_join(long_df_problems, by = c("Model", "Population_parameters_num", "Anchored", "Estimator_num", "Iteration", "Adjustment"))


classification_scenario <- left_join(combined_parameters,
                                     confounded_marginal[, c("estimator_num", "population_parameters_num", "model_marginal")],
                                     by = c("estimator_num", "population_parameters_num")) |>
  left_join(confounded_conditional[, c("estimator_num", "population_parameters_num", "model_conditional")],
            by = c("estimator_num", "population_parameters_num")) |>
  left_join(nonnormal_distribution[, c("population_parameters_num", "any_nonnormal")],
            by = "population_parameters_num")

# Long indicators
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

correct_decision <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975)*se
  lb <- coef - qnorm(0.975)*se
  contains_0 <- ub >= 0 & lb <= 0
  sign_estimate_theo_identical <- (theo >= 0 & coef >= 0) | (theo <= 0 & coef <= 0)
  theo_is_null <- theo == 0
  correct_decision <- ifelse(theo_is_null, contains_0, sign_estimate_theo_identical)
  mean(correct_decision, na.rm = TRUE)
}

df_stats <- joined_results |>
  group_by(Population_parameters_num, Estimator_num, Adjustment, Model, Anchored, Data) |>
  # Replacing NAs with warnings
  mutate(Estimate = ifelse(!is.na(Error_estimate), NA, Estimate),
          Variance = ifelse(!is.na(Error_estimate), NA, Variance)) |>
  summarize(bias = get_bias(Estimate, True_effect ),
            rmse = get_RMSE(Estimate, True_effect ),
            vr = get_VR(Estimate, sqrt(Variance)),
            cov_95 = get_cov_95(Estimate, sqrt(Variance), True_effect ),
            number_na_estimate = get_number_na_estimate(Estimate),
            correct_decision = correct_decision(Estimate, Variance, True_effect),
            .groups = "drop") |>
  pivot_longer(cols = c("bias", "rmse", "vr", "cov_95", "correct_decision", "number_na_estimate"),
               names_to = "indicator", values_to = "values")



dir.create(file.path(dir_experience_results, "processed_results"))
saveRDS(df_stats, file.path(dir_experience_results, "processed_results", "df_stats.rds"))
saveRDS(joined_results, file.path(dir_experience_results, "processed_results", "joined_results.rds"))
saveRDS(long_df_problems, file.path(dir_experience_results, "processed_results", "long_df_problems.rds"))

