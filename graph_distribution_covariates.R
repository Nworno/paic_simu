library(ggplot2)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("flat")
theme_update(plot.background = element_blank())
theme_update(panel.background = element_blank())

###################################
## Drawing distributions covariates
###################################
library(data.table)
source("env_variables.R")
source("data_generation.R")

##################################
# Drawing covariates distributions
##################################

N_pop <- 10^5
# path_experiment <- file.path("results_simulations", "test")
# if (!dir.exists(path_experiment)) dir.create(path_experiment)
path_experiment <- file.path("results_simulations", DATE_EXPERIMENT)
df_population_parameters <- readRDS(file.path(path_experiment, "df_population_parameters.RDS"))

for (num_population in df_population_parameters$population_parameters_num) {
  print(num_population)
  path_results_experiments <- file.path(path_experiment, num_population)
  if (!dir.exists(path_results_experiments)) dir.create(path_results_experiments)
  list_simulation_parameters <- df_population_parameters[df_population_parameters$population_parameters_num == num_population, ] |>
    apply(2, \(col) ifelse(is.factor(col), as.character(col), col)) |>
    unlist(recursive = TRUE)
  populations <- creating_population(list_simulation_parameters)
  pop_init <- populations$pop_init

  all_individuals <- pop_init |>
    data.table::melt(measure.vars = patterns("X[0-9]+"), value.name = "variable_value", number = as.numerical) |>
    # Because they are the only variables used for now
    dplyr::filter(variable %in% c("X1", "X2")) |>
    data.table::melt(measure.vars = c("A", "B", "C"), value.name = "Y_obs", variable.name = "ttt", number = as.numerical) |>
    dplyr::mutate(trial = factor(trial, levels = c("AC", "BC"), labels = c("AC (IPD)", "BC (AgD)")))

  plot_propensity_distribution <- all_individuals |>
    ggplot() +
    geom_density(aes(prob_imbalanced_trial, fill = trial), alpha = 0.3) +
    geom_density(aes(prob_imbalanced_trial), color = "black") + # Both trials together
    labs(x = NULL, y = NULL, title = "Propensity distributions")
  saveRDS(plot_propensity_distribution, file.path(path_results_experiments, "propensity_distribution.RDS"))
  ggsave(file.path(path_results_experiments, "propensity_distribution.png"),
         plot = plot_propensity_distribution, width = 10, height = 5)

  plot_covariates_distribution <- all_individuals |>
    ggplot() +
    geom_density(aes(variable_value, fill = trial), alpha = 0.4) +
    facet_wrap(facets = "variable", scales = "free") +
    labs(x = NULL, y = NULL, title = "Covariates distributions") +
    theme(strip.text = element_text(size = 12))
  saveRDS(plot_covariates_distribution, file.path(path_results_experiments, "covariates_distribution.RDS"))
  ggsave(file.path(path_results_experiments, "covariates_distribution.png"),
         plot = plot_covariates_distribution, width = 10, height = 5)

  if (list_simulation_parameters$outcome_distribution == "normal") {
    plot_outcome_distribution <- all_individuals |>
      ggplot() +
      # geom_density(aes(Y_obs, fill = ttt), alpha = 0.4) +
      geom_violin(aes(ttt, Y_obs, fill = ttt), alpha = 0.4) +
      geom_pointrange(aes(y = mean_Y_obs, x = ttt, ymin = low, ymax = up), color = "black", size = 1,
                      data = all_individuals |>
                        dplyr::group_by(trial, ttt) |>
                        dplyr::summarise(mean_Y_obs = mean(Y_obs), sd_Y_obs = sd(Y_obs), low = mean_Y_obs - sd_Y_obs, up = mean_Y_obs + sd_Y_obs)) +
      facet_wrap(facets = "trial", ncol = 2, scales = "free_x") +
      labs(x = NULL, y = NULL, title = "Outcome distribution")

  } else if (list_simulation_parameters$outcome_distribution == "binomial") {
    # plot the distribution of the outcome as a barplot, with proportions of the outcome as stack bars for each treatment group
    plot_outcome_distribution <- all_individuals |>
      dplyr::mutate(Y_obs = ifelse(Y_obs > 0, 1, 0), fill = ttt) |>
      dplyr::group_by(trial, ttt) |>
      dplyr::summarize("0" = 1L - mean(Y_obs), "1" = mean(Y_obs)) |>
      tidyr::pivot_longer(cols = c("0", "1"), names_to = "prop_Y_obs", values_to = "value") |>
      # dplyr::count(Y_obs) |>
      # dplyr::summarize(Y_obs = dplyr::count(Y_obs), .by = c("trial", "ttt")) |>
      ggplot() +
      geom_bar(aes(y = value, x = ttt, fill = prop_Y_obs), stat = "identity", alpha = 0.4) +
      # geom_bar(aes(Y_obs, position = "dodge", alpha = 0.4) +
      facet_wrap(facets = "trial", ncol = 1, scales = "free") +
      # geom_bar(aes(Y_obs, fill = trial), alpha = 0.4) +
      labs(x = NULL, y = NULL, title = "Outcome distribution")
  }
  # print(plot_outcome_distribution)
  print(num_population)
  print(path_results_experiments)
  saveRDS(plot_outcome_distribution, file.path(path_results_experiments, "outcomes_distribution.RDS"))
  ggsave(file.path(path_results_experiments, "outcomes_distribution.png"),
         plot = plot_outcome_distribution, width = 10, height = 5)
}


#### Graph distributions weighting
# num_population = 6
# path_results_experiments <- file.path("studying_weighting_fn", date_experiment, num_population)
# results_experiment <- readRDS(file.path(path_results_experiments, "experiment_3.RDS"))
#
# all_iterations <- lapply(results_experiment, data.table::rbindlist, fill = TRUE, idcol = "trial") |>
#   data.table::rbindlist(fill = TRUE, idcol = "iteration") |>
#   dplyr::select(iteration, trial, id, ttt, X1, X2, prob_w_trial_AC, prob_w_trial_BC, maic_w, ps_w) |>
#   tidyr::pivot_longer(
#     cols = tidyselect::matches("X[0-9]"),
#     names_to = "variable",
#     values_to = "variable_value") |>
#   dplyr::mutate(unweighted = 1) |>
#   tidyr::pivot_longer(cols = c("maic_w", "ps_w", "unweighted"), names_to = "weight_name", values_to = "weight_value") |>
#   tidyr::replace_na(list(weight_value = 1)) |>
#   dplyr::mutate(e = weight_value / (weight_value + 1))
#
#
# all_iterations |>
#   dplyr::filter(ttt %in% c("A", "B")) |>
#   dplyr::filter(variable == "X2") |>
#   ggplot() +
#   geom_density(aes(variable_value, fill = ttt, weight = weight_value), alpha = 0.4, bw = "nrd") +
#   # geom_histogram(aes(variable_value, fill = ttt, weight = weight_value), alpha = 0.4, position = "dodge") +
#   facet_wrap("weight_name") +
#   labs(x = NULL, y = NULL) +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         strip.text = element_text(size = 12))
#
#
# ### Plots weights themselves
# all_iterations |>
#   dplyr::filter(ttt %in% c("A")) |>
#   dplyr::filter(variable == "X2") |>
#   dplyr::filter(weight_name != "unweighted") |>
#   dplyr::mutate(weight_name = factor(weight_name, labels = c("maic_w" = "MAIC", "ps_w" = "IPTW"))) |>
#   ggplot() +
#   geom_point(aes(variable_value, e), alpha = 0.05, color = "black") +
#   geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   facet_grid(~weight_name) +
#   # geom_smooth(aes(variable_value, e), color = "red", method = "lm", se = FALSE, alpha = 1, fullrange = TRUE) +
#   ylim(c(0, 1)) +
#   geom_density(aes(variable_value, weight = weight_value, color = "weighted", fill = "weighted"),
#                alpha = 0.4,
#                bw = "nrd") +
#   geom_density(aes(variable_value, fill = ttt, color = ttt),
#                alpha = 0.4,
#                bw = "nrd",
#                inherit.aes = FALSE,
#                data = all_iterations[all_iterations$ttt %in% c("A", "B") &
#                                        all_iterations$variable == "X2" &
#                                        all_iterations$weight_name != "unweighted", ] |>
#                  dplyr::mutate(weight_name = factor(weight_name, labels = c("maic_w" = "MAIC", "ps_w" = "IPTW")))) +
#   scale_fill_manual(values = c("weighted" = "red", "A" = "#3498db", "B" = "#2ecc71"),
#                     name = "Trial",
#                     labels = c("weighted" = "AC weighted", "A" = "AC unweighted", "B" = "BC")) +
#   scale_color_manual(values = c("weighted" = "red", "A" = "#3498db", "B" = "#2ecc71"),
#                     name = "Trial",
#                     labels = c("weighted" = "AC weighted", "A" = "AC unweighted", "B" = "BC")) +
#   labs(y = "Propensity score", x = "X2 value") +
#   theme(axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 15),
#         strip.text = element_text(size = 15),
#         axis.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         legend.position = "bottom", legend.box = "vertical")
#
#
#
# all_iterations |>
#   dplyr::filter(ttt %in% c("A")) |>
#   dplyr::filter(variable == "X2") |>
#   dplyr::filter(weight_name != "unweighted") |>
#   ggplot() +
#   geom_point(aes(variable_value, weight_value, color = weight_name), alpha = 0.1) +
#   facet_grid(~weight_name) +
#   # geom_smooth(aes(variable_value, e, color = weight_name), alpha = 1) +
#   geom_density(aes(variable_value, fill = ttt),
#                alpha = 0.4,
#                bw = "nrd",
#                inherit.aes = FALSE,
#                data = all_iterations[all_iterations$ttt %in% c("A", "B") &
#                                        all_iterations$variable == "X2" &
#                                        all_iterations$weight_name != "unweighted", ])




# one_iteration <- results_experiment[[100]] |>
#   data.table::rbindlist(fill = TRUE) |>
#   dplyr::select(id, ttt, X1, X2, prob_w_trial_AC, prob_w_trial_BC, maic_w, ps_w) |>
#   tidyr::pivot_longer(cols = tidyselect::matches("X[0-9]"), names_to = "variable", values_to = "variable_value") |>
#   tidyr::pivot_longer(cols = c("maic_w", "ps_w"), names_to = "weight_name", values_to = "weight_value") |>
#   dplyr::mutate(weighted_value = variable_value * weight_value) |>
#   tidyr::pivot_longer(cols = c("variable_value", "weighted_value"),
#                       names_to = "weight_type",
#                       values_to = "variable_value")
#
#
#
#
# one_iteration |>
#   # dplyr::select(id, ttt, X1, X2, w) |>
#   # dplyr::mutate(across(tidyselect::matches("X[0-9]+"), as.double)) |>
#   # data.table::melt(measure.vars = patterns("X[0-9]+"), value.name = "variable_value", number = as.numerical) |>
#   ggplot() +
#   geom_density(aes(variable_value, fill = weight_type), alpha = 0.4) +
#   facet_grid(ttt ~ weight_name, scales = "free") +
#   labs(x = NULL, y = NULL) +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         strip.text = element_text(size = 12))

