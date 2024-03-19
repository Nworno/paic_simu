library(ggplot2)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("flat")
theme_update(plot.background = element_blank())
theme_update(panel.background = element_blank())

###################################
## Drawing distributions covariates
###################################
library(data.table)
options(mc.cores = 1)
source("estimators.R")
source("env_variables.R")

# df_population_parameters <- list(
#   prop_X1 = 0.5, # Variable binaire, prevalence dans la population
#   bT_X1 = c(3),   # Effet de la variable binaire sur la probabilité d'être dans l'essai AC
#   bT_X2 = c(5),   # Effet de la variable continue X2...
#   bT_X3 = 0,     # Idem, mais inutile pour le moment
#   bT_X4 = 0,     # Idem, mais inutile pour le moment
#   bY_X1 = 1.5,   # Effet de X1 sur l'outcome
#   bY_X2 = 0.5,   # Effet de X2 sur l'outcome
#   bY_X3 = 0,     # Effet de X3 sur l'outcome (inutile pour le moment)
#   bY_X4 = 0.5,     # Effet de X4 sur l'outcome
#   bY_A_X1 = c(0), # Interaction A et X1 dans le modèle outcome
#   bY_A_X2 = c(0, 1), # Interaction A et X2 dans le modèle outcome
#   bY_A_X3 = c(0), # Interaction A et X3 dans le modèle outcome
#   bY_A_X4 = c(0), # Interaction A et X4 dans le modèle outcome
#   binary_marker = c(bquote(rbinom(N_pop, 1, 0.5))), # Utilisé pour la variable bimodale
#   f_X1 = c(bquote(rbinom(N_pop, 1, 0.5))), # distribution de X1 (revoir car il faudrait utiliser prop_X1)
#   f_X2 = c(bquote(rnorm(N_pop, 0.5, 1)),
#            bquote(binary_marker * rnorm(N_pop, -2, 1) + (1 - binary_marker) * rnorm(N_pop, 2, 1))), # distribution de X2
#   f_X3 = c(bquote(0)), # inutile pour le moment
#   # f_X4 = c(bquote(binary_marker * rnorm(N_pop, -1.5, 1) + (1 - binary_marker) * rnorm(N_pop, 1.5, 1))),
#   f_X4 = c(bquote(0)), # inutile pour le moment
#   bY_A = 1.5,  # Effet de A par rapport à C
#   bY_B = 1.5,  # Effet de B par rapport à C
#   bY_C = 0,    # Pas d'effet de C sur l'outcome
#   AC_trial_model = c(bquote(X1 * bT_X1 + X2 * bT_X2 + X3 * bT_X3 + X4 * bT_X4)), # Modèle d'attribution de l'essai AC
#   BC_trial_model = c(bquote(0)),  # Modèle d'attribution de l'essai BC
#   outcome_generation_formula =  c(bquote(
#     bY_X1*X1 + bY_X2 * X2 + bY_X3 * X3 + bY_X4 * X4 + (bY_A + bY_A_X1*X1 + bY_A_X2*X2 + bY_A_X3*X3 + bY_A_X4*X4) * A +  bY_B*B + bY_C*C
#   ))
# )  |> 
#   expand.grid() |>
#   as.data.table()
# df_population_parameters[, population_parameters_num := 1:.N]


# df_population_parameters <- df_population_parameters[population_parameters_num == 6,]

######### Models
##################
## Note David: we could have
## one binary variable X1
## one continuous variable X2
## one parameter to switch the binary variable between prognostic only (bYA_X1 = 0) or effect modifier (bYA_X1 != 0)
## one parameter to switch the continuous variable between prognostic only (bYA_X2Ò = 0) or effect modifier (bYA_X2 != 0)
## one parameter to switch the continuous variable distribution: normal (symmetrical) or lognormal (asymmetrical)
## Overall, 8 scenarios here
##
## Then, run all estimators. For estimators taking into accounts covariates, run three estimations:
## - only X1
## - only X2
## - both X1 and X2
## For estimators taking into accounts moments, run with (to be discussed):
## - first moment only
## - first and second moments

##############################################
########### Creating an overarching population
##############################################

## Génère une data.frame de 10^6 ou 7 lignes
creating_population <- function(list_simulation_parameters) {
  attach(list_simulation_parameters)
  binary_marker <- rbinom(N_pop, 1, 0.5)
  pop_init <- data.table(
    id = 1:N_pop,
    X1 = eval(f_X1),
    X2 = eval(f_X2),
    X3 = rlnorm(N_pop, 0, 1), # not used
    X4 = binary_marker * rnorm(N_pop, -1.5, 1) + (1 - binary_marker) * rnorm(N_pop, 1.5, 1) # not used
  ) |>
    setkey("id")
  print(pop_init)
  
  trial_assignement_prob <- function(trial_assignment_model, df) {
    predicted <- with(df, eval(trial_assignment_model))
    return(plogis(predicted))
  }
  
  predict_outcome <- function(outcome_model, df) {
    with(df, eval(outcome_model))
  }
  
  stopifnot(imbalanced_trial %in% c("AC", "BC"))
  AC_trial_model <- switch(imbalanced_trial, "AC" = imbalanced_trial_model, "BC" = balanced_trial_model)
  BC_trial_model <- switch(imbalanced_trial, "BC" = imbalanced_trial_model, "AC" = balanced_trial_model)
  
  pop_init$prob_w_trial_AC <- trial_assignement_prob(AC_trial_model, df = pop_init)
  pop_init$prob_w_trial_BC <- trial_assignement_prob(BC_trial_model, df = pop_init)
  
  covariate_names <- c("X1", "X2", "X3", "X4")
  df_outcomes <- sapply(list(A = pop_init[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                             B = pop_init[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                             C = pop_init[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
                        predict_outcome,
                        outcome_model = outcome_generation_formula,
                        simplify = FALSE) |>
    c("id" = list(1:nrow(pop_init))) |>
    as.data.table() |>
    melt(id.vars = c("id"), variable.name = "ttt", value.name = "Y_theo")
  df_outcomes[, Y_obs := Y_theo + rnorm(n = length(Y_theo), mean = 0, sd = 1)]
  
  average_pop_init <- pop_init[, lapply(.SD, mean), .SDcols = covariate_names]
  
  average_conditional_outcome <- sapply(list("A" = average_pop_init[, .(A = 1L, B = 0L, C = 0L, (.SD)), .SDcols = covariate_names],
                                             "B" = average_pop_init[, .(A = 0L, B = 1L, C = 0L, (.SD)), .SDcols = covariate_names],
                                             "C" = average_pop_init[, .(A = 0L, B = 0L, C = 1L, (.SD)), .SDcols = covariate_names]),
                                        predict_outcome,
                                        outcome_model = outcome_generation_formula,
                                        simplify = FALSE) |>
    as.data.table() |>
    melt(measure.vars = c("A", "B", "C"), variable.name = "ttt", value.name = "outcome")
  marginal_outcome <- df_outcomes[, .(outcome = mean(Y_obs)), by = c("ttt")]
  average_outcome_df <- rbindlist(
    list("conditional" = average_conditional_outcome,
         "marginal" = marginal_outcome),
    use.names = TRUE,
    idcol = "outcome_type") |>
    dcast(outcome_type ~ ttt, value.var = "outcome")
  average_outcome_df[, AB := A - B]
  # population_variance <- df_outcomes[, .(var_Y_obs = var(Y_obs)), by = c("ttt")] |> 
  #   dcast(. ~ ttt, value.var = "var_Y_obs") |> 
  #   dplyr::rename(var_population = `.`) |> 
  #   dplyr::mutate(var_AB = A + B)
  # 
  # df_outcomes |> ggplot() + geom_violin(aes(x = Y_obs, y = ttt))
  # diff_AB <- df_outcomes[ttt == "A", Y_obs] - df_outcomes[ttt == "B", Y_obs]
  # mean(diff_AB)
  # mean((diff_AB - mean(diff_AB))^2)
  # 
  # average_outcome_df <- merge(average_outcome_df, population_variance, by = "ttt")
  # browser()
  
  # Correcting theoretical marginal effect, so that it is set to 0 when there is actually
  # no difference between theoretical conditional and marginal, as it should be
  # Useful to quantify estimators alpha and beta nominal risk level 
  if (average_outcome_df[outcome_type == "conditional", AB] == 0) average_outcome_df[, AB := 0]
  detach(list_simulation_parameters)
  return(list(
    "pop_init" = pop_init,
    "df_outcomes" = df_outcomes,
    "average_outcome_df" = average_outcome_df
  ))
}


##################################
# Drawing covariates distributions 
##################################

N_pop <- 10^5
path_experiment <- file.path("results_simulations", DATE_EXPERIMENT)
df_population_parameters <- readRDS(file.path(path_experiment, "df_population_parameters.RDS"))

for (num_population in df_population_parameters$population_parameters_num) {
  path_results_experiments <- file.path(path_experiment, num_population)
  list_simulation_parameters <- df_population_parameters[df_population_parameters$population_parameters_num == num_population, ] |>
    apply(2, \(col) ifelse(is.factor(col), as.character(col), col)) |> 
    unlist(recursive = TRUE)
  populations <- creating_population(list_simulation_parameters)
  pop_init <- populations$pop_init
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


#### Graph distributions weighting
num_population = 6
path_results_experiments <- file.path("studying_weighting_fn", date_experiment, num_population)
results_experiment <- readRDS(file.path(path_results_experiments, "experiment_3.RDS"))

all_iterations <- lapply(results_experiment, data.table::rbindlist, fill = TRUE, idcol = "trial") |> 
  data.table::rbindlist(fill = TRUE, idcol = "iteration") |>
  dplyr::select(iteration, trial, id, ttt, X1, X2, prob_w_trial_AC, prob_w_trial_BC, maic_w, ps_w) |>
  tidyr::pivot_longer(
    cols = tidyselect::matches("X[0-9]"),
    names_to = "variable",
    values_to = "variable_value") |> 
  dplyr::mutate(unweighted = 1) |> 
  tidyr::pivot_longer(cols = c("maic_w", "ps_w", "unweighted"), names_to = "weight_name", values_to = "weight_value") |> 
  tidyr::replace_na(list(weight_value = 1)) |> 
  dplyr::mutate(e = weight_value / (weight_value + 1))
  

all_iterations |> 
  dplyr::filter(ttt %in% c("A", "B")) |> 
  dplyr::filter(variable == "X2") |> 
  ggplot() +
  geom_density(aes(variable_value, fill = ttt, weight = weight_value), alpha = 0.4, bw = "nrd") +
  # geom_histogram(aes(variable_value, fill = ttt, weight = weight_value), alpha = 0.4, position = "dodge") +
  facet_wrap("weight_name") +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 12))

  
### Plots weights themselves
all_iterations |> 
  dplyr::filter(ttt %in% c("A")) |> 
  dplyr::filter(variable == "X2") |> 
  dplyr::filter(weight_name != "unweighted") |> 
  dplyr::mutate(weight_name = factor(weight_name, labels = c("maic_w" = "MAIC", "ps_w" = "IPTW"))) |> 
  ggplot() +
  geom_point(aes(variable_value, e), alpha = 0.05, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(~weight_name) +
  # geom_smooth(aes(variable_value, e), color = "red", method = "lm", se = FALSE, alpha = 1, fullrange = TRUE) +
  ylim(c(0, 1)) +
  geom_density(aes(variable_value, weight = weight_value, color = "weighted", fill = "weighted"),
               alpha = 0.4,
               bw = "nrd") +
  geom_density(aes(variable_value, fill = ttt, color = ttt),
               alpha = 0.4,
               bw = "nrd",
               inherit.aes = FALSE, 
               data = all_iterations[all_iterations$ttt %in% c("A", "B") &
                                       all_iterations$variable == "X2" &
                                       all_iterations$weight_name != "unweighted", ] |> 
                 dplyr::mutate(weight_name = factor(weight_name, labels = c("maic_w" = "MAIC", "ps_w" = "IPTW")))) +
  scale_fill_manual(values = c("weighted" = "red", "A" = "#3498db", "B" = "#2ecc71"),
                    name = "Trial", 
                    labels = c("weighted" = "AC weighted", "A" = "AC unweighted", "B" = "BC")) +
  scale_color_manual(values = c("weighted" = "red", "A" = "#3498db", "B" = "#2ecc71"),
                    name = "Trial", 
                    labels = c("weighted" = "AC weighted", "A" = "AC unweighted", "B" = "BC")) +
  labs(y = "Propensity score", x = "X2 value") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "bottom", legend.box = "vertical")



all_iterations |> 
  dplyr::filter(ttt %in% c("A")) |> 
  dplyr::filter(variable == "X2") |> 
  dplyr::filter(weight_name != "unweighted") |> 
  ggplot() +
  geom_point(aes(variable_value, weight_value, color = weight_name), alpha = 0.1) +
  facet_grid(~weight_name) +
  # geom_smooth(aes(variable_value, e, color = weight_name), alpha = 1) + 
  geom_density(aes(variable_value, fill = ttt),
               alpha = 0.4,
               bw = "nrd",
               inherit.aes = FALSE, 
               data = all_iterations[all_iterations$ttt %in% c("A", "B") &
                                       all_iterations$variable == "X2" &
                                       all_iterations$weight_name != "unweighted", ])



  
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

