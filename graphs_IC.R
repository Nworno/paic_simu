# Exploring the problem of coverage and RMSE
library("ggplot2")
library("data.table")
library("dplyr")

ggthemr::ggthemr("pale")

dir_experiment <- file.path("results_simulations", DATE_EXPERIMENT)
df_estimators_parameters <- readRDS(file.path(dir_experiment, "df_estimators_parameters.RDS"))
df_population_parameters <- readRDS(file.path(dir_experiment, "df_population_parameters.RDS"))

estimator_num <- "3"
num_experiment <- "1"
path_experiment <- file.path(file.path(dir_experiment, num_experiment))
results_experiment <- readRDS(file.path(path_experiment, paste0("experiment_", estimator_num, ".RDS")))
true_results <- readRDS(file.path(path_experiment, "average_outcome_df.RDS"))

true_effect <- true_results[outcome_type == "conditional", AB]

data_long <- results_experiment |>
  rbindlist(use.names = TRUE) |>
  group_by(adjustment, model, anchored) |>
  mutate(mean_estimate = mean(estimate),
         # sd_estimate = sd(estimate),
         lb = estimate - qnorm(0.975)*sqrt(variance),
         ub = estimate + qnorm(0.975)*sqrt(variance),
         includes_true_effect = (lb <= true_effect) & (ub >= true_effect),
         includes_0 = (lb <= 0) & (ub >= 0),
         correct_decision = (true_effect == 0 & includes_0) |
           ((true_effect != 0) & (sign(true_effect) == sign(estimate)) & !includes_0),
  )

# data_long  |>
#   # group_by(adjustment, model, anchored) |>
#   ggplot() +
#   geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.01) +
#   geom_errorbarh(aes(xmin = lb, xmax = ub, y = anchored), height = 0.1, alpha = 0.1) +
#   facet_wrap(~adjustment + model) +
#   geom_jitter(aes(x = estimate, y = anchored), alpha = 0.1) +
#   geom_vline(xintercept = true_effect, linetype = "dashed")
# # geom_vline(xintercept = 0, linetype = "dashed", colour = "black")  +


data_long  |>
  ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = true_effect, ymin = true_effect, ymax = Inf), fill = "#A0D2AD", alpha = 0.02) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), fill = "grey", alpha = 0.01) +
  geom_point(aes(x = lb, y = ub, color = anchored), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
  facet_wrap(~adjustment + model)
  # geom_vline(xintercept = c(true_effect), linetype = "dashed") +
  # geom_vline(xintercept = c(0), colour = "black", linetype = "dashed") +
  # geom_hline(yintercept = true_effect, linetype = "dashed") +
  # geom_hline(yintercept = 0, colour = "black", linetype = "dashed")

