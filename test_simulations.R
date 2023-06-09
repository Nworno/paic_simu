
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



