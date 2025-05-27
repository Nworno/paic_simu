library(ggplot2)
library(patchwork)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("fresh")
source("env_variables.R")
dir_results <- file.path("results_simulations", DATE_EXPERIMENT)
print(DATE_EXPERIMENT)

if (!dir.exists(file.path(dir_results, "plots_publication"))) dir.create(file.path(dir_results, "plots_publication"))


list_population_parameters <- as.character(1:8)
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
    experiment_dfs <- readRDS(file.path(dir_results,
                                        population_parameters_num,
                                        switch(anchored,
                                               anchored = "experiment_dfs_1.RDS",
                                               unanchored = "experiment_dfs_3.RDS")))

    base_plot <- lapply(experiment_dfs, \(x) x[[switch(anchored, unanchored = 1, anchored = 2)]]) |> # 1 is unanchored, 2 is anchored
      data.table::rbindlist() |>
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
      print("oh")
      base_plot <- tidyr::pivot_longer(base_plot, cols = all_of(var), names_to = "X_name", values_to = "X_value")
      var_col <- "X_value"
      alpha_points <- 50000/nrow(base_plot)
    } else {
      print("yeah")
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
           title = paste0("DGM-", population_parameters_num),
           fill = "Trial",
           color = NULL) +
      scale_y_continuous(name = "Propensity score", breaks = c(0, 1)) +
      scale_color_manual(values = c("PS" = "black"), labels = c("PS" = expression("PS in " * italic(a) * " (IPD) trial"))) +
      scale_fill_manual(values = c("a" = "#2ecc71", "b" = "#f1c40f"),
                        labels = c("a" = expression(italic(a) *" (IPD)"),
                                   "b" = expression(italic(b * "(AgD)")))
      ) +
      guides(color = guide_legend(override.aes = list(alpha = 1, color = 'black'))) +
      theme(text = element_text(size = 16),
            strip.text = element_text(size = 16)) +
      expand_limits(y = c(0, 1)) +
      theme(legend.position = "bottom",
            legend.box = "vertical") +
      (if (length(var) > 1) labs(X = "X value") else labs(X = var))
  }

  plot_propensity_distributions <- list_plots[["1"]] +
    list_plots[["2"]] +
    list_plots[["3"]] +
    list_plots[["4"]] +
    list_plots[["5"]] +
    list_plots[["6"]] +
    list_plots[["7"]] +
    list_plots[["8"]] +
    plot_layout(ncol = 2, guides = "collect") &
    theme(plot.title = element_text(size = 20),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text = element_text(size = 16))
  ggsave(filename = file.path(dir_results, "plots_publication", paste0("PS_distributions_all_", var_col, "_", anchored, ".jpeg")), # for some reasons, pdf is much larger for that one
         width = 15, height = 20,
         plot = plot_propensity_distributions)

}


