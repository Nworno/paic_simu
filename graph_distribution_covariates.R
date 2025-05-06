library(ggplot2)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("flat")
# theme_update(plot.background = element_blank())
# theme_update(panel.background = element_blank())

###################################
## Drawing distributions covariates
###################################
library(data.table)
library(patchwork)
source("env_variables.R")
source("data_generation.R")

GRAPH_PUBLICATION <- TRUE

N_pop <- 10^5
path_experiment <- file.path("results_simulations", DATE_EXPERIMENT)
if (!dir.exists(file.path(path_experiment, "plots_publication"))) {
  dir.create(file.path(path_experiment, "plots_publication"))
}
df_population_parameters <- readRDS(file.path(path_experiment, "df_population_parameters.RDS"))

###############################################################
# Drawing covariates distributions, one scenario after the other
###############################################################

scenario_of_interest <- 1:6
plots_distribution_covariates <- list()
for (num_population in scenario_of_interest) {
  print(num_population)
  path_results_experiments <- file.path(path_experiment, num_population)
  if (!dir.exists(path_results_experiments)) dir.create(path_results_experiments)
  list_simulation_parameters <- df_population_parameters[df_population_parameters$population_parameters_num == num_population, ] |>
    apply(2, \(col) ifelse(is.factor(col), as.character(col), col)) |>
    unlist(recursive = TRUE)
  populations <- creating_population(list_simulation_parameters)
  pop_init <- populations$pop_init
  average_outcome_df <- populations$average_outcome_df
  average_outcome_df[, trial := factor(trial, levels = c("AC", "BC"), labels = c("AC (IPD)", "BC (AgD)"))]

  all_individuals <- pop_init |>
    data.table::melt(measure.vars = patterns("X[0-9]+"), value.name = "variable_value", number = as.numerical) |>
    data.table::melt(measure.vars = c("A", "B", "C"), value.name = "Y_obs", variable.name = "ttt", number = as.numerical) |>
    dplyr::mutate(trial = factor(trial, levels = c("AC", "BC"), labels = c("AC (IPD)", "BC (AgD)")))

  plot_propensity_distribution <- all_individuals |>
    ggplot() +
    geom_density(aes(prob_BC, fill = trial), alpha = 0.3) +
    labs(x = NULL, y = NULL, title = "Propensity distributions")
  saveRDS(plot_propensity_distribution, file.path(path_results_experiments, "propensity_distribution.RDS"))
  ggsave(file.path(path_results_experiments, "propensity_distribution.png"),
         plot = plot_propensity_distribution, width = 10, height = 5)

  plot_covariates_distribution <- all_individuals |>
    ggplot() +
    geom_density(aes(variable_value, fill = trial), alpha = 0.4) +
    facet_grid(rows = "variable", scales = "free") +
    labs(x = NULL, y = NULL, title = "Covariates distributions") +
    theme(strip.text = element_text(size = 12))
  saveRDS(plot_covariates_distribution, file.path(path_results_experiments, "covariates_distribution.RDS"))
  ggsave(file.path(path_results_experiments, "covariates_distribution.png"),
         plot = plot_covariates_distribution, width = 10, height = 5)


  plots_distribution_covariates[[num_population]] <- all_individuals[variable == "X1", ] |>
    ggplot() +
    geom_density(aes(variable_value, fill = trial, color = trial), alpha = 0.7) +
    labs(x = NULL, y = NULL, fill = "Trial", color = "Trial", title = paste0("DGM-", num_population)) +
    scale_fill_manual(values = c("AC (IPD)" = "#2ecc71", "BC (AgD)" = "#f1c40f"),
                      labels = c("AC (IPD)" = expression(italic(a) *  "(IPD)"),
                                 "BC (AgD)" = expression(italic(b) * "(AgD)"))
                      ) +
    scale_color_manual(values = c("AC (IPD)" = "#2ecc71", "BC (AgD)" = "#f1c40f"),
                      labels = c("AC (IPD)" = expression(italic(a) * "(IPD)"),
                                 "BC (AgD)" = expression(italic(b) * "(AgD)"))
                      ) +
    theme(
      title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "bottom",
      panel.background = element_blank(),
    )

  saveRDS(plot_covariates_distribution, file.path(path_results_experiments, "covariates_distribution.RDS"))
  # ggsave(file.path(path_results_experiments, "covariates_distribution.png"),
  #        plot = plot_covariates_distribution, width = 10, height = 5)

  df_outcome <- readRDS(file.path(path_results_experiments, "average_outcome_df.RDS"))
  if (list_simulation_parameters$outcome_distribution == "normal") {

    # Rajouter tableau df_outcome

    plot_outcome_distribution <- all_individuals |>
      ggplot() +
      # geom_density(aes(Y_obs, fill = ttt), alpha = 0.4) +
      geom_violin(aes(ttt, Y_obs, fill = ttt), alpha = 0.4) +
      geom_pointrange(aes(y = mean_Y_obs, x = ttt, ymin = low, ymax = up), color = "black", size = 1,
                      data = all_individuals |>
                        dplyr::group_by(trial, ttt) |>
                        dplyr::summarise(mean_Y_obs = mean(Y_obs), sd_Y_obs = sd(Y_obs), low = mean_Y_obs - sd_Y_obs, up = mean_Y_obs + sd_Y_obs)) +
      geom_point(aes(ttt, Y_obs, color = "conditional"), size = 3, data = average_outcome_df[outcome_type == "conditional", ] |>
                   melt(measure.vars = c("A", "B", "C"), value.name = "Y_obs", variable.name = "ttt")) +
      facet_wrap(facets = "trial", ncol = 2, scales = "free_x") +
      scale_color_manual(values = c("conditional" = "red"), name = NULL, labels = c("conditional" = "Conditional Effect")) +
      labs(x = NULL, y = NULL, title = "Marginal outcome distribution")

  } else if (list_simulation_parameters$outcome_distribution == "binomial") {
    # plot the distribution of the outcome as a barplot, with proportions of the outcome as stack bars for each treatment group
    # browser()
    marginal_outcome_distribution <- all_individuals |>
      dplyr::mutate(Y_obs = ifelse(Y_obs > 0, 1, 0), fill = ttt) |>
      dplyr::group_by(trial, ttt) |>
      dplyr::summarize("0" = 1L - mean(Y_obs), "1" = mean(Y_obs)) |>
      tidyr::pivot_longer(cols = c("0", "1"), names_to = "prop_Y_obs", values_to = "value")

    conditional_outcome_distribution <- average_outcome_df[outcome_type == "conditional", ] |>
      melt(measure.vars = c("A", "B", "C"), id.vars = c("trial"), variable.name = "ttt") |>
      data.table:::DT(, .(value_0 = 1 - value, trial, ttt, value)) |>
      melt(measure.vars = c("value", "value_0"), variable.name = "prop_Y_obs") |>
      data.table:::DT(, .(prop_Y_obs = factor(prop_Y_obs, levels = c("value_0", "value"), labels = c("0", "1")), trial, ttt, value))

    plot_outcome_distribution <- rbindlist(list("marginal" = marginal_outcome_distribution,
                                                "conditional" = conditional_outcome_distribution),
                                           use.names = TRUE,
                                           idcol = "distribution",
                                           ) |>
      ggplot() +
      geom_bar(aes(y = value, x = ttt, fill = prop_Y_obs), stat = "identity", alpha = 0.7) +
      scale_fill_manual(values = c("a" = "#65ADC2", "b" = "233B43"),
                        labels = c("a" = expression(italic(a)),
                                   "b" = expression(italic(b))
                        )) +
      # geom_bar(aes(Y_obs, position = "dodge", alpha = 0.4) +
      facet_grid("distribution ~ trial") +
      # geom_bar(aes(Y_obs, fill = trial), alpha = 0.4) +
      labs(x = NULL, y = NULL, title = paste0("DGM-", num_population))

  }
  # print(plot_outcome_distribution)
  print(num_population)
  print(path_results_experiments)
  saveRDS(plot_outcome_distribution, file.path(path_results_experiments, "outcomes_distribution.RDS"))
  # ggsave(file.path(path_results_experiments, "outcomes_distribution.png"),
  #        plot = plot_outcome_distribution, width = 10, height = 5)
}
g <- plots_distribution_covariates[[1]] +
  plots_distribution_covariates[[2]] +
  plots_distribution_covariates[[3]] +
  plots_distribution_covariates[[4]] +
  plots_distribution_covariates[[5]] +
  plots_distribution_covariates[[6]] +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  patchwork::plot_annotation(title  = "Covariates (X1 and X2) distributions in a and b trials")
ggplot2::ggsave(file.path(path_experiment, "plots_publication", "all_covariates_distribution.pdf"), g)


# plot_weighting <- function(df, weight_column) {
#   plot_distribution <- df |>
#     dplyr::select(-X3, -X4) |>
#     tidyr::pivot_longer(cols = c("X1", "X2"), names_to = "variable", values_to = "values") |>
#     ggplot() +
#     geom_density(aes(values, fill = trial, col = trial, weight = {{ weight_column }}), alpha = 0.5) +
#     facet_wrap(facets = "variable")
#   return(plot_distribution)
# }
#
# list_dfs <- readRDS(file.path("results_simulations/20241016_184147", "1", paste0("experiment_dfs_", "3", ".RDS")))
# for (num_population in df_population_parameters$population_parameters_num) {
#   for (num_estimator in df_estimators_parameters$estimator_num) {
#     list_dfs <- readRDS(file.path(path_experiment, num_population, paste0("experiment_dfs_", num_estimator, ".RDS")))
#     sample_list_dfs <- list_dfs[sample(1:length(list_dfs), min(length(list_dfs), 3), replace = FALSE)]
#
#     list_weighting_plots <- sapply(sample_list_dfs,
#                                    \(iteration) {
#                                      # sapply(iteration,
#                                      # \(sublist) {
#                                      non_weighted <- plot_weighting(iteration[[1]], weight_column = NULL)
#                                      ml <- plot_weighting(iteration[[1]], weight_column = ml)
#                                      maic_1 <- plot_weighting(iteration[[1]], weight_column = maic_1)
#                                      maic_2 <- plot_weighting(iteration[[1]], weight_column = maic_2)
#                                      return(list("non_weighted" = non_weighted,
#                                                  "ml" = ml,
#                                                  "maic_1" = maic_1,
#                                                  "maic_2" = maic_2))
#                                    },
#                                    USE.NAMES = TRUE,
#                                    simplify = FALSE)
#     saveRDS(list_weighting_plots,
#             file = file.path(file.path(path_experiment, num_population, paste0("sample_weighting_plots", num_estimator, ".RDS"))))
#   }
# }