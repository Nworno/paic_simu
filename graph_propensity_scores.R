library(ggplot2)
library(patchwork)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("fresh")
source("env_variables.R")
dir_results <- file.path("results_simulations", DATE_EXPERIMENT)
print(DATE_EXPERIMENT)

if (!dir.exists(file.path(dir_results, "plots_publication"))) dir.create(file.path(dir_results, "plots_publication"))


df_pop_params <- readRDS(file.path(dir_results, "df_population_parameters.RDS"))
list_population_parameters <- as.character(sort(unique(df_pop_params$population_parameters_num)))

# Scenarios sharing identical PS distributions (same covariate distribution and trial assignment)
covariate_groups <- list(c("1", "7", "9"), c("2", "8", "10"))
scenario_group_label    <- setNames(paste0("DGM-", list_population_parameters), list_population_parameters)
scenario_representative <- setNames(list_population_parameters, list_population_parameters)
for (grp in covariate_groups) {
  label <- paste0("DGM-", paste(grp, collapse = ", "))
  for (s in grp) {
    scenario_group_label[s]    <- label
    scenario_representative[s] <- grp[1]
  }
}
unique_scenarios <- list_population_parameters[sapply(list_population_parameters, \(s) scenario_representative[s] == s)]

# Plot propensity distribution --------
list_plots <- list()
for (config in 1:2) {
  if (config == 1) {
    var <- c("X1", "X2")
    anchored <- "unanchored"
  } else {
    var <- "X1"
    anchored <- "anchored"
  }
  for (population_parameters_num in list_population_parameters) {
    if (scenario_representative[population_parameters_num] != population_parameters_num) next
    experiment_dfs <- readRDS(file.path(dir_results,
                                        population_parameters_num,
                                        switch(anchored,
                                               anchored = "experiment_dfs_1.RDS",
                                               unanchored = "experiment_dfs_3.RDS")))

    combined_data <- lapply(experiment_dfs, \(x) x[[switch(anchored, unanchored = 1, anchored = 2)]]) |> # 1 is unanchored, 2 is anchored
      data.table::rbindlist()
    if (nrow(combined_data) == 0 || !all(c(var, "prob_BC", "trial", "ttt", "ml", "maic_1", "maic_2") %in% names(combined_data))) {
      message("Skipping DGM ", population_parameters_num, " (", anchored, "): missing data or columns")
      next
    }
    base_plot <- combined_data |>
      dplyr::select(all_of(var), true_PS = prob_BC, trial, ttt, ml, maic_1, maic_2) |>
      dplyr::mutate(unweighted = 1,
                    true_logit_ps = log(true_PS / (1 - true_PS))) |>
      tidyr::pivot_longer(cols = c("ml", "maic_1", "maic_2", "unweighted"),
                          names_to = "weight_type", values_to = "weights") |>
      dplyr::mutate(ps = weights / (weights + 1),
                    trial = dplyr::case_match(trial,
                                              "AC" ~ "a",
                                              "BC" ~ "b"),
                    weight_type = dplyr::case_match(weight_type,
                                                    "ml" ~ "PSW",
                                                    "maic_1" ~ "MAIC-1",
                                                    "maic_2" ~ "MAIC-2",
                                                    "unweighted" ~ "Unweighted") |>
                      factor(levels = c("Unweighted", "PSW", "MAIC-1", "MAIC-2")))
    if (length(var) > 1)  {
      base_plot <- tidyr::pivot_longer(base_plot, cols = all_of(var), names_to = "X_name", values_to = "X_value")
      var_col <- "X_value"
      alpha_points <- 50000/nrow(base_plot)
    } else {
      var_col <- var
      alpha_points <- 10000/nrow(base_plot)
    }
    list_plots[[population_parameters_num]] <- ggplot(base_plot) +
      geom_point(
        data = base_plot |> dplyr::filter(trial == "a" & weight_type != "Unweighted"),
        aes(x = .data[[var_col]], y = ps, color = "PS"),
        # alpha = alpha_points
        alpha = 0.01
      ) +
      geom_density(aes(x = .data[[var_col]], fill = trial, weight = weights, color = trial), alpha = 0.7) +
      (if (length(var) > 1) facet_grid(X_name ~ weight_type) else facet_wrap(~weight_type)) +
      labs(x = var,
           y = "Propensity Score",
           title = scenario_group_label[population_parameters_num],
           fill = "Trial",
           color = NULL) +
      scale_y_continuous(name = "Propensity score", breaks = c(0, 1)) +
      scale_color_manual(values = c("PS" = "black"), labels = c("PS" = expression("PS in " * italic(a) * " (IPD) trial"))) +
      scale_fill_manual(values = c("a" = "#2ecc71", "b" = "#f1c40f"),
                        labels = c("a" = expression(italic(a) *" (IPD)"),
                                   "b" = expression(italic(b) * " (AgD)"))
      ) +
      guides(color = guide_legend(override.aes = list(alpha = 1, color = 'black'))) +
      theme(text = element_text(size = 16),
            strip.text = element_text(size = 16)) +
      expand_limits(y = c(0, 1)) +
      theme(legend.position = "bottom",
            legend.box = "vertical") +
      (if (length(var) > 1) labs(x = "X value") else labs(x = var))
  }

  valid_plots <- Filter(Negate(is.null), list_plots[unique_scenarios])
  if (length(valid_plots) == 0) next
  n_panels    <- length(valid_plots)
  plot_height <- ceiling(n_panels / 2) * 5
  plot_propensity_distributions <- patchwork::wrap_plots(valid_plots) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(plot.title = element_text(size = 20),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text = element_text(size = 16))
  ggsave(filename = file.path(dir_results, "plots_publication", paste0("PS_distributions_all_", var_col, "_", anchored, ".jpeg")),
         width = 15, height = plot_height,
         plot = plot_propensity_distributions)

}


