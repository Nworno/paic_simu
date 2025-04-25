library(ggplot2)
library(magrittr, include.only = "%$%")
library(patchwork)
library(fontawesome)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("fresh")
source("env_variables.R")
DATE_EXPERIMENT <- "20241016_184147" # Results simulations
# DATE_EXPERIMENT <- "20250401_161025" # Propensity score distribution
dir_results <- file.path("results_simulations", DATE_EXPERIMENT)
print(DATE_EXPERIMENT)
joined_results <- readRDS(file.path(dir_results, "processed_results", "joined_results.rds"))
df_stats <- readRDS(file.path(dir_results, "processed_results", "df_stats.rds"))

if (!dir.exists(file.path(dir_results, "plots_publication"))) dir.create(file.path(dir_results, "plots_publication"))

df_estimators_parameters <- readRDS(file.path("~/git_repos/paic_simu/results_simulations/", DATE_EXPERIMENT, "/df_estimators_parameters.RDS")) |>
  dplyr::rename(Estimator_num = estimator_num)

df_population_parameters <- readRDS(file.path("~/git_repos/paic_simu/results_simulations/", DATE_EXPERIMENT, "/df_population_parameters.RDS"))

na_counts <- df_stats |>
  dplyr::filter(indicator == "number_na_estimate") |>
  dplyr::rename(missing_count = values) |>
  dplyr::filter(missing_count > 0) |>
  dplyr::select(Population_parameters_num, Estimator_num, Adjustment, Model, Anchored, Data, missing_count)

df_stats <- df_stats |>
  dplyr::left_join(na_counts,
                   by = c("Population_parameters_num", "Estimator_num", "Adjustment", "Model", "Anchored", "Data"))




# Treatment effect, only weighting methods, only complete models --------
for (population_parameters_num in c("1", "2", "3", "4")) {

  df_plot <- df_stats |>
    dplyr::filter(Population_parameters_num == population_parameters_num,
                  ((Estimator_num == 3 & Anchored == "Unanchored") | (Estimator_num == 1 & Anchored == "Anchored")), # plotting complete models only here
    ) |>
    dplyr::filter(!indicator %in% c("correct_decision", "number_na_estimate") &
                    `Adjustment` %in% c("Unadjusted", "IPTW")) |>
    dplyr::mutate(
      objective = ifelse(indicator %in% c("bias", "rmse"), 0,
                         ifelse(indicator == "vr", 1,
                                ifelse(indicator == "cov_95", 0.95, NA))),
      Model = dplyr::case_match(Model,
                                "MAIC_1" ~ "MAIC-1",
                                "MAIC_2" ~ "MAIC-2",
                                "ML" ~ "PSW",
                                "Unadjusted" ~ "Unadjusted"),
      indicator = dplyr::case_match(indicator,
                                    "bias" ~ "Bias",
                                    "rmse" ~ "RMSE",
                                    "vr" ~ "VR",
                                    "cov_95" ~ "95% coverage") |> factor(levels = c("Bias", "RMSE", "VR", "95% coverage"))
    )
  plot_results_DGM <- ggplot(df_plot) +
    facet_wrap(~indicator, scales = "free_x") +
    geom_point(aes(y = Model, x = values, color = Anchored, shape = Anchored), size = 7, alpha = 0.7) +
    geom_text(data = dplyr::filter(df_plot, !is.na(missing_count)),
              aes(x = values, y = Model, label = paste0("*")),
              vjust = 0, hjust = 0, size = 8, color = "darkred") +
    scale_shape_manual(breaks = c("Anchored", "Unanchored"), values = c(17, 16)) +
    geom_vline(aes(xintercept = objective), linetype = "dashed", color = "black") +
    labs(x = NULL, y = NULL, color = NULL, shape = NULL) +
    theme_bw() +
    theme(text = element_text(size = 16),
          strip.text = element_text(size = 16),
          plot.caption = element_text(size = 12, hjust = 0))

  if (any(!is.na(df_plot$missing_count))) warning(paste0("Missing values for population_parameters_num: ", population_parameters_num))

   ggsave(filename = file.path(dir_results, "plots_publication", paste0("performance_estimators", population_parameters_num, ".pdf")),
         width = 15, height = 7,
         plot = plot_results_DGM)
}




# Plots exploring confounding with missing variables --------
df_bias <- df_stats |>
  dplyr::filter(Model %in% c("MAIC_1",
                             # "MAIC_2",
                             "ML"),
                df_stats$Population_parameters_num %in% c("1", "2", "3", "4"),
                indicator == "bias") |>
  dplyr::mutate(
    Covariates = dplyr::case_match(Estimator_num,
                                   1 ~ "X1",
                                   2 ~ "X2",
                                   3 ~ "X1, X2") |> factor(levels = c("X1, X2", "X2", "X1")),
    values = abs(values), # In order to align the biases across all DGM
    objective = 0,
    Model = dplyr::case_match(Model,
                              "MAIC_1" ~ "MAIC-1",
                              "MAIC_2" ~ "MAIC-2",
                              "ML" ~ "PSW",
                              "Unadjusted" ~ "Unadjusted"),
    is_complete = ifelse(Anchored == "Anchored" & Estimator_num %in% c(1, 3), "complete",
                         ifelse(Anchored == "Unanchored" & Estimator_num == 3, "complete",
                                "incomplete")),
    indicator = "Bias",
    Covariates_offset = as.numeric(factor(Covariates)) +
      ifelse(Model == "MAIC-1", -0.05, 0.05))
fixed_scale_x <- with(df_bias, c(min(abs(values), na.rm = TRUE), max(abs(values), na.rm = TRUE)))

list_plots <- list()
for (population_parameters_num in c("1", "2", "3", "4")) {
  df_plot <- df_bias |>
    dplyr::filter(Population_parameters_num == population_parameters_num)
  covariate_labels <- levels(factor(df_plot$Covariates))

  list_plots[[population_parameters_num]] <- ggplot(df_plot) +
    geom_point(
      aes(
        x = values,
        y = Covariates_offset,
        shape = Anchored,
        fill = Model,
        colour = is_complete
      ),
      size = 7,
      alpha = 0.6,
      stroke = 1.5
    ) +
    # Dummy point just for the red triangle in legend
    geom_point(
      data = data.frame(x = NA_real_, Covariates_offset = NA_real_, is_complete = "incomplete"),
      aes(x = x, y = Covariates_offset, colour = is_complete),
      shape = 24,
      size = 7,
      stroke = 1.5,
      fill = "white"
    ) +
    geom_text(data = dplyr::filter(df_plot, !is.na(missing_count)),
              aes(x = values, y = Covariates_offset, label = paste0("*")),
              vjust = 0, hjust = 0, size = 9, color = "darkred") +
    scale_shape_manual(
      values = c("Anchored" = 24, "Unanchored" = 21)  # Use triangle for unanchored
    ) +
    scale_fill_manual(
      values = c("MAIC-1" = "#109B37",
                 # "MAIC-2" = "#168E7F",
                 "PSW" = "#C29365"),
      name = "Model"
    ) +
    scale_colour_manual(
      values = c("complete" = "black", "incomplete" = "red"),
      name = "Confounding",
      labels = c("incomplete" = "Incomplete model")  # Custom label
    ) +
    geom_vline(
      aes(xintercept = objective),
      linetype = "dashed",
      color = "black"
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, colour = "black")),  # Show model fill with consistent appearance
      colour = guide_legend(override.aes = list(shape = 22, fill = "white"))  # Show triangle with red border
    ) +
    # Relabel y-axis
    scale_y_continuous(
      breaks = 1:length(covariate_labels),
      labels = covariate_labels
    ) +
    scale_x_continuous(limits = fixed_scale_x, expand = c(0, 0.1)) + #TODO: reprendre d'ici, comprendre pourquoi ne fonctionne pas
    labs(title = paste0("DGM-", population_parameters_num), x = NULL, y = "Covariates", shape = "Anchoring") +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 16))

  if (any(!is.na(df_plot$missing_count))) warning(paste0("Missing values for population_parameters_num: ", population_parameters_num))
}

plot_bias_confounding <- list_plots[["1"]] +
  list_plots[["2"]] +
  list_plots[["3"]] +
  list_plots[["4"]] +
  plot_layout(ncol = 2,
              guides = "collect") &
  theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.box = "vertical")
ggsave(filename = file.path(dir_results, "plots_publication", paste0("bias_confounding.pdf")),
       width = 15, height = 12,
       plot = plot_bias_confounding)



# Plot propensity distribution --------
list_plots <- list()
for (population_parameters_num in c("1", "2", "3", "4")) {
  # var <- c("X1", "X2")
  # anchored <- "unanchored"
  var <- "X1"
  anchored <- "anchored"
  experiment_dfs <- readRDS(file.path(dir_results,
                                      population_parameters_num,
                                      switch(anchored,
                                             anchored = "experiment_dfs_1.RDS",
                                             unanchored = "experiment_dfs_3.RDS")))

  base_plot <- lapply(experiment_dfs, \(x) x[[switch(anchored, unanchored = 1, anchored = 2)]]) |> # 1 is unanchored, 2 is anchored
    data.table::rbindlist() |>
    dplyr::select(all_of(var), true_PS = prob_BC, trial, ttt, ml, maic_1, maic_2) |>
    dplyr::mutate(unweighted = 1) |>
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
    alpha_points <- 10000/nrow(base_plot)
  } else {
    var_col <- var
    alpha_points <- 5000/nrow(base_plot)
  }
  list_plots[[population_parameters_num]] <- ggplot(base_plot) +
    geom_point(
      data = base_plot |> dplyr::filter(trial == "a" & weight_type != "Unweighted"),
      aes(x = .data[[var_col]], y = ps, color = "PS"),
      alpha = alpha_points
    ) +
    geom_density(aes(x = .data[[var_col]], fill = trial, weight = weights, color = trial), alpha = 0.7) +
    (if (length(var) > 1) facet_grid(X_name ~ weight_type) else facet_wrap(~weight_type)) +
    labs(x = var,
         y = "Propensity Score",
         title = paste0("DGM-", population_parameters_num),
         fill = "Trial",
         color = NULL) +
    scale_color_manual(values = c("PS" = "black"), labels = c("PS" = expression("PS in " * italic(a) * " (IPD) trial"))) +
    scale_fill_manual(values = c("a" = "#2ecc71", "b" = "#f1c40f"),
                      labels = c("a" = expression(italic(a) *" (IPD)"),
                                 "b" = expression(italic(b * "(AgD)")))
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, color = 'black'))) +
    theme(text = element_text(size = 16),
          strip.text = element_text(size = 16)) +
    ylim(0, 1) +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    (if (length(var) > 1) labs(X = "X value") else labs(X = var))
}

plot_propensity_distributions <- list_plots[["1"]] +
  list_plots[["2"]] +
  list_plots[["3"]] +
  list_plots[["4"]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.box = "vertical")
ggsave(filename = file.path(dir_results, "plots_publication", paste0("PS_distributions_all_", var_col, "_", anchored, ".jpeg")), # for some reasons, pdf is much larger for that one
       width = 15, height = 12,
       plot = plot_propensity_distributions)


