library(ggplot2)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("flat")

###################################
## Drawing distributions covariates
###################################
library(data.table)
library(patchwork)
source("env_variables.R")
N_pop <- 10^6
source("data_generation.R")

GRAPH_PUBLICATION <- TRUE

path_experiment <- file.path("results_simulations", DATE_EXPERIMENT)
if (!dir.exists(file.path(path_experiment, "plots_publication"))) {
  dir.create(file.path(path_experiment, "plots_publication"))
}
df_population_parameters <- readRDS(file.path(path_experiment, "df_population_parameters.RDS"))

###############################################################
# Drawing covariates distributions, one scenario after the other
###############################################################

scenario_of_interest <- as.character(sort(unique(df_population_parameters$population_parameters_num)))

# Scenarios sharing identical covariate distributions (same f_X1/X2 and trial assignment)
covariate_groups <- list(c("1", "7", "9"), c("2", "8", "10"))
scenario_group_label    <- setNames(paste0("DGM-", scenario_of_interest), scenario_of_interest)
scenario_representative <- setNames(scenario_of_interest, scenario_of_interest)
for (grp in covariate_groups) {
  label <- paste0("DGM-", paste(grp, collapse = ", "))
  for (s in grp) {
    scenario_group_label[s]    <- label
    scenario_representative[s] <- grp[1]
  }
}
unique_scenarios <- scenario_of_interest[sapply(scenario_of_interest, \(s) scenario_representative[s] == s)]

plots_distribution_covariates <- list()
for (num_population in scenario_of_interest) {
  if (scenario_representative[num_population] != num_population) next
  print(num_population)
  path_results_experiments <- file.path(path_experiment, num_population)
  if (!dir.exists(path_results_experiments)) dir.create(path_results_experiments)
  list_simulation_parameters <- df_population_parameters[df_population_parameters$population_parameters_num == num_population, ] |>
    apply(2, \(col) ifelse(is.factor(col), as.character(col), col)) |>
    unlist(recursive = TRUE)
  populations <- creating_population(list_simulation_parameters, N_pop)
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

  plot_covariates_distribution <- all_individuals |>
    ggplot() +
    geom_density(aes(variable_value, fill = trial), alpha = 0.4) +
    facet_grid(rows = "variable", scales = "free") +
    labs(x = NULL, y = NULL, title = "Covariates distributions") +
    theme(strip.text = element_text(size = 12))

  plots_distribution_covariates[[num_population]] <- all_individuals[variable == "X1", ] |>
    ggplot() +
    geom_density(aes(variable_value, fill = trial, color = trial), alpha = 0.7) +
    labs(x = NULL, y = NULL, fill = "Trial", color = "Trial", title = scenario_group_label[num_population]) +
    scale_fill_manual(values = c("AC (IPD)" = "#2ecc71", "BC (AgD)" = "#f1c40f"),
                      labels = c("AC (IPD)" = expression(italic(a) *  "(IPD)"),
                                 "BC (AgD)" = expression(italic(b) * "(AgD)"))
                      ) +
    scale_color_manual(values = c("AC (IPD)" = "#2ecc71", "BC (AgD)" = "#f1c40f"),
                      labels = c("AC (IPD)" = expression(italic(a) * "(IPD)"),
                                 "BC (AgD)" = expression(italic(b) * "(AgD)"))
                      ) +
    theme(
      title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.position = "bottom",
      panel.background = element_blank(),
    )

  df_outcome <- readRDS(file.path(path_results_experiments, "average_outcome_df.RDS"))
  if (list_simulation_parameters$outcome_distribution == "normal") {

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
    marginal_outcome_distribution <- all_individuals |>
      dplyr::mutate(Y_obs = ifelse(Y_obs > 0, 1, 0), fill = ttt) |>
      dplyr::group_by(trial, ttt) |>
      dplyr::summarize("0" = 1L - mean(Y_obs), "1" = mean(Y_obs)) |>
      tidyr::pivot_longer(cols = c("0", "1"), names_to = "prop_Y_obs", values_to = "value")

    conditional_outcome_distribution <- average_outcome_df[outcome_type == "conditional", ] |>
      melt(measure.vars = c("A", "B", "C"), id.vars = c("trial"), variable.name = "ttt") |>
      (\(dt) dt[, .(value_0 = 1 - value, trial, ttt, value)])() |>
      melt(measure.vars = c("value", "value_0"), variable.name = "prop_Y_obs") |>
      (\(dt) dt[, .(prop_Y_obs = factor(prop_Y_obs, levels = c("value_0", "value"), labels = c("0", "1")), trial, ttt, value)])()

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
      facet_grid("distribution ~ trial") +
      labs(x = NULL, y = NULL, title = paste0("DGM-", num_population))

  }
  print(num_population)
  print(path_results_experiments)
}

n_panels    <- length(plots_distribution_covariates)
plot_height <- ceiling(n_panels / 2) * 4
g <- patchwork::wrap_plots(plots_distribution_covariates[unique_scenarios]) +
  plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom") &
  patchwork::plot_annotation(title = "Covariates (X1 and X2) distributions in a and b trials")
ggplot2::ggsave(plot = g,
                filename = file.path(path_experiment, "plots_publication", "all_covariates_distribution.pdf"),
                width = 10, height = plot_height)
