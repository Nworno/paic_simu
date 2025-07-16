library(ggplot2)
library(magrittr, include.only = "%$%")
library(patchwork)
library(fontawesome)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("fresh")
source("env_variables.R")
# DATE_EXPERIMENT <- "20241016_184147" # Results simulations
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


list_population_parameters <- as.character(1:8)

# Treatment effect, only weighting methods, only complete models --------
for (population_parameters_num in list_population_parameters) {

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

# Graph with pair of treatments
pair_list_population_parameters <- list(
  "1" = c("1", "2"),
  "2" = c("3", "4"),
  "3" = c("5", "6"),
  "4" = c("7", "8")
)
for (pair_population_parameters_num in pair_list_population_parameters) {

  df_plot <- df_stats |>
    dplyr::filter(Population_parameters_num %in% pair_population_parameters_num,
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
  if (any(!is.na(df_plot$missing_count))) warning(paste0("Missing values for population_parameters_num: ", pair_population_parameters_num))
  list_plots <- lapply(pair_population_parameters_num, function(elem) {
     sub_df_plot <- df_plot |> dplyr::filter(Population_parameters_num == elem)
      ggplot(sub_df_plot) +
      facet_wrap(~indicator, scales = "free_x") +
      geom_point(aes(y = Model, x = values, color = Anchored, shape = Anchored), size = 7, alpha = 0.7) +
      geom_text(data = dplyr::filter(sub_df_plot, !is.na(missing_count)),
                aes(x = values, y = Model, label = paste0("*")),
                vjust = 0, hjust = 0, size = 8, color = "darkred") +
      scale_shape_manual(breaks = c("Anchored", "Unanchored"), values = c(17, 16)) +
      geom_vline(aes(xintercept = objective), linetype = "dashed", color = "black") +
      labs(x = NULL, y = NULL, color = NULL, shape = NULL) +
      theme_bw() +
      theme(text = element_text(size = 16),
            strip.text = element_text(size = 16),
            plot.caption = element_text(size = 12, hjust = 0)) +
      labs(title = paste0("DGM-", elem))
  })

  list_plots[[1]] + list_plots[[2]] +
    plot_layout(ncol = 1, guides = "collect") &
    theme(plot.title = element_text(size = 20),
          legend.position = "bottom",
            legend.box = "vertical")
    ggsave(filename = file.path(dir_results, "plots_publication", paste0("performance_estimators_pair", paste(pair_population_parameters_num, collapse = "_"),  ".pdf")),width = 12, height = 15)
}



# Plots exploring confounding with missing variables --------
df_bias <- df_stats |>
  dplyr::filter(Model %in% c("MAIC_1",
                             # "MAIC_2",
                             "ML"),
                df_stats$Population_parameters_num %in% list_population_parameters,
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
for (population_parameters_num in list_population_parameters) {
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
    labs(title = paste0("DGM-", population_parameters_num), x = "Bias", y = "Covariates", shape = "Anchoring") +
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
  list_plots[["5"]] +
  list_plots[["6"]] +
  list_plots[["7"]] +
  list_plots[["8"]] +
  plot_layout(ncol = 2,
              guides = "collect") &
  theme(plot.title = element_text(size = 20),
        legend.position = "bottom",
        legend.box = "vertical")
ggsave(filename = file.path(dir_results, "plots_publication", paste0("bias_confounding.pdf")),
       width = 15, height = 16,
       plot = plot_bias_confounding)


