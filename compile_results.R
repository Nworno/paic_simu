library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemr)
ggthemr::ggthemr("pale")
results_simulations <- readRDS("results_simulations_sacha_3/results_simulations.RDS")

long_results <- data.table::rbindlist(results_simulations, idcol ="n_iter")[
  , data := ifelse(model %in% c("maic", "stc"),
                  "PAIC",
                  ifelse(model %in% c("ml", "glm"),
                         "IPD",
                         ifelse(model == "unadjusted", "AgD", NA)))]

true_conditional <- long_results[name == "true" & model == "conditional", unique(estimate)]
true_marginal <- long_results[name == "true" & model == "marginal", unique(estimate)]

df_true_effects <- data.frame(name = c("unadjusted", "regression", "iptw"), true = c(true_conditional, true_conditional, true_marginal)) |> data.table::as.data.table()


# Long results overall
long_results |>
  filter(name != "true") |>
  ggplot(aes(x = estimate, y = name, color = data, fill = data)) +
  # geom_point(position = "jitter") +
  geom_violin() +
  facet_wrap(~anchored) +
  geom_vline(xintercept = true_conditional, linetype = "dashed")

# Long indicators
get_bias <- function(obs, theo) {
  mean(obs - theo, na.rm = FALSE)
}
get_RMSE <- function(obs, theo) {
  sqrt(mean((obs - theo)**2, na.rm = FALSE))
}
get_RV <- function(obs, se_obs) {
  mean(se_obs, na.rm = FALSE) / sd(obs, na.rm = FALSE)
}
get_cov_95 <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975)*se
  lb <- coef - qnorm(0.975)*se
  covered <- theo < ub & theo > lb
  mean(covered, na.rm = FALSE)
}

df_stats <- df_true_effects[long_results[name != "true"], on = "name"][, .(bias = get_bias(estimate, true),
                                                                           rmse = get_RMSE(estimate, true),
                                                                           rv = get_RV(estimate, sqrt(variance)),
                                                                           cov_95 = get_cov_95(estimate, sqrt(variance), true)),
                                                                       by = .(name, model, anchored, data)] |>
  tidyr::pivot_longer(cols = c("bias", "rmse", "rv", "cov_95"),
                      names_to = "indicator", values_to = "values")

df_stats |> ggplot(aes(x = values, color = data, shape = anchored)) +
  geom_point(aes(y = "test"), position = "jitter", size =3) +
  # geom_vline(aes(xintercept = intercept),
  #            data = data.frame(indicator = c("bias", "cov_95", "rmse", "rv"),
  #                              intercept = c(0, 0.95, 0, 1)),
  #            linetype = "dashed") +
  facet_grid(cols = vars(indicator), rows = vars(name), scales = "free", shrink = TRUE)


