library(data.table)
library(dplyr)
logit <- function(x) log(x/(1 - x))

N_pop <- 10^6
N_RCT <- 200

prop_X1 <- 0.5
bT_X1 <- 2

bY_X1 <- 1.5
bYA_X1 <- 1.2

bY_A <- 1.5
bY_B <- 1.5
bY_C <- 0

trial_assignment_model_AC <- expr(X1 * bT_X1)
trial_assignment_model_BC <- expr(0) # required that prob of trial assignment independent of baseline characteristics of overarching population for the ATT to represent the marginal effect in the initial (overall) population
outcome_model <- expr(bYA_X1*X1*A + bY_X1*X1 + bY_A*A + bY_B*B + bY_C*C)

# Creating an overarching population
pop_init <- data.table(
  id = 1:N_pop,
  X1 = rbinom(N_pop, 1, prop_X1)
) |>
  setkey("id")

trial_assignement_prob <- function(df, trial_assignment_model, deviates = FALSE) {
  predicted <- with(df, eval(trial_assignment_model))
  # variance of a logistic function is set ((s^2\pi^2)/3), with s=1 when using normal law -->
  # making an assumption about bias and residual heterogeneity
  return(plogis(predicted))
}

predict_outcome <- function(df, outcome_model, deviates = FALSE) {
  predicted <- with(df, eval(outcome_model))
  if (deviates) return(predicted + rnorm(n = length(predicted), mean = 0, sd = 1)) else return(predicted)
}

pop_init$prob_AC <- trial_assignement_prob(pop_init, trial_assignment_model_AC)
pop_init$prob_BC <- trial_assignement_prob(pop_init, trial_assignment_model_BC)


pop_init$Y_theo_A <- copy(pop_init)[, c("A", "B", "C") := .(1L, 0L, 0L)] |>
  predict_outcome(outcome_model)
pop_init$Y_theo_B <- copy(pop_init)[, c("A", "B", "C") := .(0L, 1L, 0L)] |>
  predict_outcome(outcome_model)
pop_init$Y_theo_C <- copy(pop_init)[, c("A", "B", "C") := .(0L, 0L, 1L)] |>
  predict_outcome(outcome_model)

ttt_names <- c("A", "B", "C")
cols_obs <- paste0("Y_obs_", ttt_names)
cols_theo <- paste0("Y_theo_", ttt_names)
# Assigning "observed" values for all the individual in the population
pop_init[, eval(expr(cols_obs)) := lapply(.SD, \(x) x + rnorm(n = length(x))), .SD = cols_theo]
# equivalent
# pop_init[, sub("theo", "obs", cols_theo, fixed = TRUE) := lapply(.SD, \(x) x + rnorm(n = length(x))), .SD = cols_theo]


pop_init_long <- pop_init |> tidyr::pivot_longer(
  cols = starts_with("Y_"),
  names_to = c("measure_type", "ttt"),
  names_prefix = "Y_",
  names_sep = "_"
) |>
  tidyr::pivot_wider(names_from = "measure_type", values_from = "value") |>
  # C ref levels for subsequent models
  dplyr::mutate(ttt = factor(ttt, levels = c("C", "A", "B"))) |>
  data.table::setDT()
# may prefer alternative straightforward implementation using data.table
# pop_init_long <- melt(pop_init,
#                       id.vars = "id",
#                       measure.vars = patterns("theo", "obs"),
#                       value.name = c("theo", "obs"),
#                       variable.name = "Y",
#                       variable.factor = FALSE
# )


# Drawing trials
# Could do all the draws at once, and transform to long df with bind_rows, and then do one massive join

trial_AC <- pop_init[
  sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_AC), # because key(pop_init) == id
  .(id, ttt = factor(rep_len(c("A", "C"), length.out = N_RCT), levels = c("C","A", "B")))
][pop_init_long[,c("id", "X1", "ttt", "prob_AC", "theo", "obs")], on = .(id, ttt), nomatch = NULL]
trial_BC <- pop_init[
  sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_BC), # because key(pop_init) == id
  .(id, ttt = factor(rep_len(c("B", "C"), length.out = N_RCT), levels = c("C","A", "B")))
][pop_init_long[,c("id", "X1", "ttt", "prob_BC", "theo", "obs")], on = .(id, ttt), nomatch = NULL]
trials_combined <- data.table::rbindlist(list("AC" = trial_AC, "BC" = trial_BC), fill = TRUE, idcol = "trial")[,trial:= factor(trial, levels = c("AC", "BC"))]


# Marginal true effect in pop_init
true_marginal_AB <- mean(pop_init$Y_theo_A) - mean(pop_init$Y_theo_B)

# (Average) Conditional true effect in pop_init
true_conditional_AB <- eval(outcome_model, envir = list(X1 = prop_X1, A = 1, B = 0, C = 0)) -
  eval(outcome_model, envir = list(X1 = prop_X1, A = 0, B = 1, C = 0))


##############################################################
########## REGRESSION BASED OUTCOME MODEL (both treatment IPD)
##############################################################

# Anchored observed conditional effect with IPD (in population similar to BC trial)
model_AC <- glm(formula = obs ~ X1*ttt, family = "gaussian", data = trial_AC)
model_BC <- glm(formula = obs ~ X1*ttt, family = "gaussian", data = trial_BC)
obs_conditional_AC <- model_AC$coefficients[["tttA"]] +
  model_AC$coefficients[["X1:tttA"]]*mean(trial_BC$X1)
obs_conditional_BC <- model_BC$coefficients[["tttB"]] +
  model_BC$coefficients[["X1:tttB"]]*mean(trial_BC$X1)

obs_conditional_AB <- obs_conditional_AC - obs_conditional_BC
# Same but with mixed model (ie random effect meta analysis on individual patients)
#TODO: reprendre d'ici
lme4::lmer(formula = obs ~ X1*ttt + (ttt|trial), family = "gaussian", data = trials_combined)
model_AC <- glm(formula = obs ~ X1*ttt, family = "gaussian", data = trials_combined)
model_BC <- glm(formula = obs ~ X1*ttt, family = "gaussian", data = trial_BC)
obs_conditional_AC <- model_AC$coefficients[["tttA"]] +
  model_AC$coefficients[["X1:tttA"]]*mean(trial_BC$X1)
obs_conditional_BC <- model_BC$coefficients[["tttB"]] +
  model_BC$coefficients[["X1:tttB"]]*mean(trial_BC$X1)



# Unanchored observed conditional effect with IPD
# observed effect under population with distribution similar to what is observed in the BC trial
model_unanchored <- glm(formula = obs ~ X1*ttt, family = "gaussian", data = trials_combined)
observed_conditional_A <- model_unanchored$coefficients[["tttA"]] + model_unanchored$coefficients[["X1:tttA"]]*mean(trial_BC$X1)
observed_conditional_B <- model_unanchored$coefficients[["tttB"]] + model_unanchored$coefficients[["X1:tttB"]]*mean(trial_BC$X1)
unanchored_conditional_AB <- observed_conditional_A - observed_conditional_B




#################################################
########### PROPENSITY SCORE (both treatment IPD)
#################################################

# (Anchored) observed marginal effect (distribution of X1 similar to the observed BC trial)
trials_combined$PS_trial <- glm(trial ~ X1, trials_combined, family = binomial(link = "logit"))$fitted.values

# 'ATC'/'ATT' weights, ie weighting so AC trial gets reweighted following BC trial distribution
trials_combined[, AT_BC_w := (trial == "BC") + (trial == "AC") * PS_trial / (1 - PS_trial)]

# Anchored
obs_marginal_AC <- trials_combined[trial == "AC" & ttt == "A", sum(obs * AT_BC_w) / sum(AT_BC_w)] - trials_combined[trial == "AC" & ttt == "C", sum(obs * AT_BC_w) / sum(AT_BC_w)]
obs_marginal_BC <- trials_combined[trial == "BC" & ttt == "B", mean(obs)] - trials_combined[trial == "BC" & ttt == "C", mean(obs)]

obs_anchored_marginal_AB <- obs_marginal_AC - obs_marginal_BC

# Unanchored
obs_unanchored_marginal_AB <- trials_combined[trial == "AC" & ttt == "A", sum(obs * AT_BC_w) / sum(AT_BC_w)] -
  trials_combined[trial == "BC" & ttt == "B", mean(obs)]


##############################################
########## Naive unadjusted conditional effect
##############################################

# (Anchored) unadjusted observed conditional effect in pop_init
unadjusted_conditional_model_AC <- glm(obs ~ ttt, family = "gaussian", data = trial_AC)
unadjusted_conditional_model_BC <- glm(obs ~ ttt, family = "gaussian", data = trial_BC)

unadjusted_conditional_AB <- unadjusted_conditional_model_AC$coefficients[["tttA"]] - unadjusted_conditional_model_BC$coefficients[["tttB"]]

# (Anchored) unadjusted observed marginal effect in pop_init
unadjusted_marginal_AC <- trial_AC[ttt == "A", mean(obs)] - trial_AC[ttt == "C", mean(obs)]
unadjusted_marginal_BC <- trial_BC[ttt == "B", mean(obs)] - trial_BC[ttt == "C", mean(obs)]

unadjusted_marginal_AB <- unadjusted_marginal_AC - unadjusted_marginal_BC

### MAIC
Ag_trials_BC <- trial_BC[, .("mean_X1" = mean(X1), "sd_X1" = sd(X1),
                             "mean_theo" = mean(theo), "sd_theo" = sd(theo),
                             "mean_obs" = mean(obs), "sd_obs" = sd(obs)), by = ttt]



# Unanchored analysis
unadjusted_relative_ttt_effect <- model_AC$coefficients[[X1]]


# ] |> merge(pop_init_long[,c("id", "ttt", "theo", "obs")], by = c("id", "ttt"), all.x = TRUE, all.y = FALSE)



# observed_AC <- pop_init[trial_AC, .(id, X1, A, B, C, Y_A_theo, Y_C_theo)][
#   , Y_obs := predict_outcome(.SD, outcome_model, deviates = TRUE)
# ]
#
# observed_BC <- pop_init[trial_BC, .(id, X1, A, B, C, Y_B_theo, Y_C_theo)][
#   , Y_obs := predict_outcome(.SD, outcome_model, deviates = TRUE)
# ]
#
# cols_ttt <-
# pop_init[trial_BC, .(id, X1, A, B, C, Y_A_theo, Y_B_theo, Y_C_theo)][
#   , paste0("Y_", c("A", "B", "C"), "_obs") := lapply(.SD, rnorm(cols_ttt), )
# ]

###########################################################################
# TESTS -------------------------------------------------------------------
###########################################################################

microbenchmark::microbenchmark(plogis(-100:100))

# Microbenchmarking different ways to multiply columns --------------------
## Ranked from fastest to slowest
microbenchmark::microbenchmark(
  mapply(FUN = "*", source_pop[, .(V1, V2, V3)], c(coef_V1, coef_V2, coef_V3))
)
microbenchmark::microbenchmark(
  source_pop[, .(odds_trial = as.matrix(.SD[, c(V1, V2, V3)]) %*% c(coef_V1, coef_V2, coef_V3))]
)
microbenchmark::microbenchmark(
  sweep(source_pop[, .(V1, V2, V3)], MARGIN = 2, STATS = c(coef_V1, coef_V2, coef_V3), FUN = "*")
)


# Benchmarking ways to do matrix multiplications --------------------------

microbenchmark::microbenchmark(
  source_pop [, .(odds_trial = trial_odds_model(V1, V2, V3))]
)

microbenchmark::microbenchmark(
  as.matrix(source_pop[, .(V1, V2, V3)]) %*% c(coef_V1, coef_V2, coef_V3)
)


AC <- source_pop[sample(.N, 1000, prob = proba), ]
BC <- source_pop[sample(.N, 1000, prob = 1- proba), ]


source_pop[, sapply(X = .SD, is.character), .SDcols = c("V1", "V2")]



# old ---------------------------------------------------------------------

library(simstudy)
N_subject <- 100000


# standard DGM found in the genCorData package
# anchored setting
coef_age <- -0.01
coef_female <- -1
coef_visits <- 2
coef_Tr <- 1
coef_EM_female <- 0.1
coef_trial <- 0.2 # only used if Trial = 1 and treatment = 1, thus represent the relative difference between the two 'active' treatments
trial_odds_model <- function(age, female, visits) {
  coef_age*age + coef_female*female + coef_visits*visits
}

predict_outcome <- function(age, female, visits, trial, Tr) {
  #with binary effect modifiers
  coef_age*age + coef_female*female + coef_visits*visits +
    Tr*(coef_Tr + coef_EM_female*female) + # Main treatment effect + modifiers
    Tr*trial*coef_trial #TODO should this term also include an interaction term?
}
C <- matrix(c(1, 0.7, 0.2, 0.7, 1, 0.8, 0.2, 0.8, 1), nrow = 3)
non_correlated <- diag(4)
DGM <- tibble::tribble(
  ~varname, ~formula                        , ~variance, ~distri, ~link,
  "age"   ,	"50"	                           , 2       ,	"normal",	"identity",
  "female",	"-2  + age * 0.02"              , 0       ,	"binary",	"logit",
  "visits",	"1.5  - 0.02 * age + 0.5 * female", 0       ,	"poisson",	"log",
) %>% as.data.table()
source_pop <- simstudy::genData(n = 100000, dtDefs = DGM)
summary(source_pop)
source_pop <- genCorData(N_subject,
                         mu = c(0.4, -0.8, 0.3, 0.1),
                         sigma = c(0.1, 0.2, 0.3, 0.1),
                         cnames = c("V1", "V2", "V3", "Tr"),
                         corMatrix = non_correlated)

source_pop[, odds_trial := trial_odds_model(age, female, visits)]
source_pop[, trial := rbinom(length(odds_trial), 1, logistic(odds_trial))]
source_pop[, outcome := predict_outcome(age, female, visits, trial, Tr)]
# addCorData(idname = "cid", mu = c(0, 0), sigma = c(2, 0.2), rho = -0.2,
#            corstr = "cs", cnames = c("a0", "a1")) %>%

# Trick for faster sampling of many subsets, from https://stackoverflow.com/questions/5458271/fast-sampling-in-r
# -- I needed to divide a 235,000,000 population into random sets of 150,000 each. At first I tried sampling the sets individually, but that would've taken over a day, so I sampled only once: population <- sample(population, length(population)), then chunked. --> Divided by 50 time of the code
# source_pop <- source_pop[, .(.SD, proba = logistic(rowSums(.SD))), .SDcols = c("V1", "V2", "V3")]
# source_pop <- source_pop[, .(.SD, proba = logistic(sum(.SD))), .SDcols = c(".SD.V1", ".SD.V2", ".SD.V3"), by = .EACHI]



