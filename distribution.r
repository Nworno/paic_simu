library(data.table)

options(mc.cores = 1)
source("estimators.R")
N_pop <- 10^6
N_BOOT_ITER <- 10
# Parameters
####################

df_population_parameters <- list(
  prop_X1 = 0.5, # Variable binaire, prevalence dans la population
  bT_X1 = c(3),   # Effet de la variable binaire sur la probabilité d'être dans l'essai AC
  bT_X2 = c(5),   # Effet de la variable continue X2...
  bT_X3 = 0,     # Idem, mais inutile pour le moment
  bT_X4 = 0,     # Idem, mais inutile pour le moment
  bY_X1 = 1.5,   # Effet de X1 sur l'outcome
  bY_X2 = 0.5,   # Effet de X2 sur l'outcome
  bY_X3 = 0,     # Effet de X3 sur l'outcome (inutile pour le moment)
  bY_X4 = 0.5,     # Effet de X4 sur l'outcome
  bY_A_X1 = c(0), # Interaction A et X1 dans le modèle outcome
  bY_A_X2 = c(0, 1), # Interaction A et X2 dans le modèle outcome
  bY_A_X3 = c(0), # Interaction A et X3 dans le modèle outcome
  bY_A_X4 = c(0), # Interaction A et X4 dans le modèle outcome
  binary_marker = c(bquote(rbinom(N_pop, 1, 0.5))), # Utilisé pour la variable bimodale
  f_X1 = c(bquote(rbinom(N_pop, 1, 0.5))), # distribution de X1 (revoir car il faudrait utiliser prop_X1)
  f_X2 = c(bquote(rnorm(N_pop, 0.5, 1)),
           bquote(binary_marker * rnorm(N_pop, -2, 1) + (1 - binary_marker) * rnorm(N_pop, 2, 1))), # distribution de X2
  f_X3 = c(bquote(0)), # inutile pour le moment
  # f_X4 = c(bquote(binary_marker * rnorm(N_pop, -1.5, 1) + (1 - binary_marker) * rnorm(N_pop, 1.5, 1))),
  f_X4 = c(bquote(0)), # inutile pour le moment
  bY_A = 1.5,  # Effet de A par rapport à C
  bY_B = 1.5,  # Effet de B par rapport à C
  bY_C = 0,    # Pas d'effet de C sur l'outcome
  AC_trial_model = c(bquote(X1 * bT_X1 + X2 * bT_X2 + X3 * bT_X3 + X4 * bT_X4)), # Modèle d'attribution de l'essai AC
  BC_trial_model = c(bquote(0)),  # Modèle d'attribution de l'essai BC
  outcome_generation_formula =  c(bquote(
    bY_X1*X1 + bY_X2 * X2 + bY_X3 * X3 + bY_X4 * X4 + (bY_A + bY_A_X1*X1 + bY_A_X2*X2 + bY_A_X3*X3 + bY_A_X4*X4) * A +  bY_B*B + bY_C*C
  ))
)  |>
  expand.grid() |>
  as.data.table()
df_population_parameters[, population_parameters_num := 1:.N]
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
  ) |>
    setkey("id")


  trial_assignement_prob <- function(trial_assignment_model, df) {
    predicted <- with(df, eval(trial_assignment_model))
    return(plogis(predicted))
  }

  predict_outcome <- function(outcome_model, df) {
    with(df, eval(outcome_model))
  }

  pop_init[, prob_w_trial_AC := trial_assignement_prob(AC_trial_model, df = pop_init)]
  pop_init[, prob_w_trial_BC := trial_assignement_prob(BC_trial_model, df = pop_init)]

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

  return(list(
    "pop_init" = pop_init,
    "df_outcomes" = df_outcomes,
    "average_outcome_df" = average_outcome_df
  ))
}


cat("Scenario ", 1, "\n")
pop <- creating_population(unlist(df_population_parameters[1, ]))
# pop <- creating_population(unlist(df_population_parameters[30, ]))
pop_init <- pop$pop_init
df_outcomes <- pop$df_outcomes

selected_individuals_AC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_AC)][
  , ttt := rep_len(c("A", "C"), length.out = .N)]
# selected_outcomes_AC <- df_outcomes[selected_individuals_AC, on = c("id", "ttt")][
# , c("id", "ttt", "Y_obs")]
# trial_AC <- pop_init[selected_outcomes_AC, on = "id"][, ttt := factor(ttt, levels = c("C", "A"))]
trial_AC <- selected_individuals_AC

selected_individuals_BC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_BC)][
  , ttt := rep_len(c("B", "C"), length.out = .N)]
# selected_outcomes_BC <- df_outcomes[selected_individuals_BC, on = c("id", "ttt")][
# , c("id", "ttt", "Y_obs")]
# trial_BC <- pop_init[selected_outcomes_BC, on = "id"][, ttt := factor(ttt, levels = c("C", "B"))]
trial_BC <- selected_individuals_BC
trial_AC$trial <- "AC"
trial_BC$trial <- "BC"
trials <- rbind(trial_AC, trial_BC)
trials$scenario <- 1
trials$bT_X2 <- paste0("bT_X2=", formatC(unlist(df_population_parameters[1, "bT_X2"]), digits= 1, format = "f"))
trials$f_X2 <- paste0("f_X2=", df_population_parameters[1, "f_X2_char"])

library(ggplot2)
ggplot(data = pop$pop_init, aes(x = X2)) + geom_density()
ggplot(data = pop_init, aes(x = X2)) + geom_density()
ggplot(data = trials, aes(x = X2, color = trial)) + geom_density()


for (i in 2:nrow(df_population_parameters)) {
  cat("Scenario ", i, "\n")
  pop <- creating_population(unlist(df_population_parameters[i, ]))
  pop_init <- pop$pop_init
  df_outcomes <- pop$df_outcomes

  selected_individuals_AC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_AC)][
    , ttt := rep_len(c("A", "C"), length.out = .N)]
  # selected_outcomes_AC <- df_outcomes[selected_individuals_AC, on = c("id", "ttt")][
  # , c("id", "ttt", "Y_obs")]
  # trial_AC <- pop_init[selected_outcomes_AC, on = "id"][, ttt := factor(ttt, levels = c("C", "A"))]
  trial_AC <- selected_individuals_AC

  selected_individuals_BC <- pop_init[sample(id, N_RCT, replace = FALSE, prob = prob_w_trial_BC)][
    , ttt := rep_len(c("B", "C"), length.out = .N)]
  # selected_outcomes_BC <- df_outcomes[selected_individuals_BC, on = c("id", "ttt")][
  # , c("id", "ttt", "Y_obs")]
  # trial_BC <- pop_init[selected_outcomes_BC, on = "id"][, ttt := factor(ttt, levels = c("C", "B"))]
  trial_BC <- selected_individuals_BC

  trial_AC$trial <- "AC"
  trial_BC$trial <- "BC"
  tmp <- rbind(trial_AC, trial_BC)
  tmp$scenario <- i
  tmp$bT_X2 <- paste0("bT_X2=", formatC(unlist(df_population_parameters[i, "bT_X2"]), digits= 1, format = "f"))
  tmp$f_X2 <- paste0("f_X2=", df_population_parameters[i, "f_X2_char"])
  trials <- rbind(trials, tmp)
}

trials$f_X2f <- factor(trials$f_X2, paste0("f_X2=", unique(df_population_parameters$f_X2_char)))
trials$bT_X2f <- factor(trials$bT_X2, paste0("bT_X2=", unique(formatC(unlist(df_population_parameters$bT_X2), digits= 1, format = "f"))))

dummy <- trials %>%
  group_by(bT_X2f, f_X2f, trial) %>%
  summarize(mean = mean(X2), median = median(X2), sd = sd(X2))
dummy$x <- ifelse(dummy$trial == "AC", -4, 3)

trials$facet <- interaction(trials$bT_X2f, trials$f_X2f, drop = TRUE)
plots <- levels(trials$facet)
coolplots <- plots[grepl("rnormbimod\\(", plots) & grepl("bT_X2=-?[1-2]", plots)]
arnaudplots <- plots[(grepl("rnorm\\(", plots) | grepl("rlnorm\\(", plots)) & grepl("bT_X2=0.2", plots)]

print(
  ggplot(data = trials, aes(x = X2, color = trial)) +
    geom_density() +
    facet_grid(bT_X2f~f_X2f) +
    geom_vline(data = dummy, aes(xintercept = mean, color = trial)) +
    geom_vline(data = dummy, aes(xintercept = median, color = trial), linetype = "dotted") +
    geom_text(data = dummy, aes(x = x, y = -Inf, label = round(sd, 2)), hjust = -0.1, vjust = -1) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 0.7)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              data = ~ subset(., facet %in% coolplots),
              colour = "red", fill = NA, inherit.aes = FALSE) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              data = ~ subset(., facet %in% arnaudplots),
              colour = "blue", fill = NA, inherit.aes = FALSE)
)
ggsave("density.pdf", width = 8, height = 10)
print(
  ggplot(data = trials, aes(x = X2, fill = trial, color = trial)) + geom_histogram(bins = 50, alpha = 0.1) + facet_grid(bT_X2f~f_X2f) +
    geom_vline(data = dummy, aes(xintercept = mean, color = trial)) +
    geom_vline(data = dummy, aes(xintercept = median, color = trial), linetype = "dotted") +
    geom_text(data = dummy, aes(x = x, y = -Inf, label = round(sd, 2)), hjust = -0.1, vjust = -1) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1000)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              data = ~ subset(., facet %in% coolplots),
              colour = "red", fill = NA, inherit.aes = FALSE) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              data = ~ subset(., facet %in% arnaudplots),
              colour = "blue", fill = NA, inherit.aes = FALSE)
)
ggsave("histogram.pdf", width = 8, height = 10)
