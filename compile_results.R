library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggthemr)
library(patchwork)
library(data.table)
ggthemr::ggthemr("pale")
results_simulations <- readRDS("results_simulations_sacha_5/results_simulations.RDS")

long_results <- bind_rows(results_simulations, .id = "n_iter") |>
  mutate(data = ifelse(model %in% c("maic", "stc"),
                       "PAIC",
                       ifelse(
                         model %in% c("ml", "glm"),
                         "IPD",
                         ifelse(
                           model == "unadjusted", "AgD", NA)
                         )
                       ),
         data = factor(data, levels = c("AgD", "PAIC", "IPD")),
         anchored = factor(anchored, levels = c("unanchored", "anchored")),
         adjustment = ifelse(adjustment == "iptw", "IPTW", stringr::str_to_title(adjustment)),
         across(c(model, anchored), .fns = stringr::str_to_title)) |>
  rename_with(stringr::str_to_title)


true_conditional <- long_results |>
  filter(Adjustment == "True" & Model == "Conditional") |>
  distinct(Estimate) |>
  pull(Estimate)
true_marginal <- long_results |>
  filter(Adjustment == "True" & Model == "Marginal") |>
  distinct(Estimate) |>
  pull(Estimate)

df_true_effects <- data.frame(Adjustment = c("Unadjusted", "Regression", "IPTW"),
                              true = c(true_conditional, true_conditional, true_marginal))


# Variation of the Estimate
long_results |>
  filter(Adjustment != "True") |>
  ggplot(aes(x = Estimate, y = Adjustment, color = Data, fill = Data)) +
  geom_violin(alpha = 0.8) +
  facet_wrap(~Anchored) +
  geom_vline(aes(xintercept = true_conditional, linetype = type),
             show.legend = TRUE,
             data = data.frame(effect = true_conditional, type = "dashed")) +
  scale_linetype_manual(values = c("dashed"), labels = c(str_wrap("True conditional effect", width = 10))) +
  labs(y = NULL,
       title = "Treatment effect estimates")

# Variance
long_results |>
  filter(Adjustment != "True") |>
  ggplot(aes(x = Variance, y = Adjustment, color = Data, fill = Data)) +
  geom_violin() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~Anchored) +
  labs(y = NULL,
       x = "Standard error of the estimate",
       title = str_to_title("Treatment effect estimate standard error"))


# Long indicators
get_bias <- function(obs, theo) {
  mean(obs - theo, na.rm = FALSE)
}
get_RMSE <- function(obs, theo) {
  sqrt(mean((obs - theo)**2, na.rm = FALSE))
}
get_VR <- function(obs, se_obs) {
  mean(se_obs, na.rm = FALSE) / sd(obs, na.rm = FALSE)
}
get_cov_95 <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975)*se
  lb <- coef - qnorm(0.975)*se
  covered <- theo < ub & theo > lb
  mean(covered, na.rm = FALSE)
}

df_stats <- long_results |> filter(Adjustment != "True") |>
  left_join(df_true_effects, by = "Adjustment") |>
  group_by(Adjustment, Model, Anchored, Data) |>
  summarize(bias = get_bias(Estimate, true),
            rmse = get_RMSE(Estimate, true),
            vr = get_VR(Estimate, sqrt(Variance)),
            cov_95 = get_cov_95(Estimate, sqrt(Variance), true)) |>
  pivot_longer(cols = c("bias", "rmse", "vr", "cov_95"),
               names_to = "indicator", values_to = "values")


bias_plot <- ggplot(filter(df_stats, indicator == "bias"),
                    aes(x = values, color = Data, shape = Anchored)) +
  geom_point(aes(y = Adjustment), position = "jitter", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Bias")
rmse_plot <- ggplot(filter(df_stats, indicator == "rmse"),
                    aes(x = values, color = Data, shape = Anchored)) +
  geom_point(aes(y = Adjustment), position = "jitter", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "RMSE")
vr_plot <- ggplot(filter(df_stats, indicator == "vr"),
                    aes(x = values, color = Data, shape = Anchored)) +
  geom_point(aes(y = Adjustment), position = "jitter", size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "VR")
ci_plot <- ggplot(filter(df_stats, indicator == "cov_95"),
                    aes(x = values, color = Data, shape = Anchored)) +
  geom_point(aes(y = Adjustment), position = "jitter", size = 3) +
  geom_vline(xintercept = 0.95, linetype = "dashed") +
  labs(title = "95% coverage")

bias_plot + rmse_plot + vr_plot + ci_plot +
  plot_annotation(title = "Performance of the estimators") +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom',
        legend.box = "vertical",
        axis.title = element_blank())


