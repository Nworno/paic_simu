library(ggplot2)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("flat")
# Graphs distributions confoundings across scenario

date_experiment <- "20240109_122850"
n_populations <- list.dirs(path = file.path("results_simulations", date_experiment),
                           recursive = FALSE, 
                           full.names = FALSE) |> length()
for (num_population in 1:n_populations) {
  path_results_experiments <- file.path("results_simulations", date_experiment, num_population)
  pop_init <- readRDS(file.path(path_results_experiments, "pop_init.RDS"))
  N_RCT <- 2*10^5
  selected_individuals_AC <- pop_init[sample(id, N_RCT, replace = TRUE, prob = prob_w_trial_AC)][
    , ttt := rep_len(c("A", "C"), length.out = .N)]
  selected_individuals_BC <- pop_init[sample(id, N_RCT, replace = TRUE, prob = prob_w_trial_BC)][
    , ttt := rep_len(c("B", "C"), length.out = .N)]
  
  covariates_distribution <- data.table::rbindlist(list("AC" = selected_individuals_AC, "BC" = selected_individuals_BC),
                                                   idcol = "trial") |> 
    dplyr::mutate(across(tidyselect::matches("X[0-9]+"), as.double)) |> 
    data.table::melt(measure.vars = patterns("X[0-9]+"), value.name = "variable_value", number = as.numerical) |> 
    # Because they are the only variables used for now
    dplyr::filter(variable %in% c("X1", "X2")) |> 
    ggplot() +
    geom_density(aes(variable_value, fill = trial), alpha = 0.4) +
    facet_wrap(facets = "variable", scales = "free") +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_blank(),
           axis.ticks = element_blank(),
           strip.text = element_text(size = 12))
  saveRDS(covariates_distribution, file.path(path_results_experiments, "covariates_distribution.RDS"))
  ggsave(file.path(path_results_experiments, "covariates_distribution.png"),
         plot = covariates_distribution, width = 10, height = 5)
}
