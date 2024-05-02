library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggthemr)
library(patchwork)
library(data.table)
library(ggrepel)

source("processing_results.R")
# View(classification_scenario)
path_results <- file.path("results_simulations", DATE_EXPERIMENT, "processed_results")

get_bias <- function(obs, theo) {
  mean(obs - theo, na.rm = FALSE)
}
get_RMSE <- function(obs, theo) {
  sqrt(mean((obs - theo)**2, na.rm = FALSE))
}
get_VR <- function(obs, se_obs) {
  mean(se_obs, na.rm = FALSE) / sd(obs, na.rm = FALSE)
}
get_cov_95 <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975)*se
  lb <- coef - qnorm(0.975)*se
  covered <- theo < ub & theo > lb
  mean(covered, na.rm = FALSE)
}

correct_decision <- function(coef, se, theo) {
  ub <- coef + qnorm(0.975)*se
  lb <- coef - qnorm(0.975)*se
  contains_0 <- ub >= 0 & lb <= 0
  sign_estimate_theo_identical <- (theo >= 0 & coef >= 0) | (theo <= 0 & coef <= 0)
  theo_is_null <- theo == 0
  correct_decision <- ifelse(theo_is_null, contains_0, sign_estimate_theo_identical)
  mean(correct_decision, na.rm = FALSE)
}


df_stats <- joined_results |>
  group_by(Population_parameters_num, Estimator_num, Adjustment, Model, Anchored, Data) |>
  summarize(bias = get_bias(Estimate, True_effect ),
            rmse = get_RMSE(Estimate, True_effect ),
            vr = get_VR(Estimate, sqrt(Variance)),
            cov_95 = get_cov_95(Estimate, sqrt(Variance), True_effect ), 
            correct_decision = correct_decision(Estimate, Variance, True_effect),
            .groups = "drop") |>
  pivot_longer(cols = c("bias", "rmse", "vr", "cov_95", "correct_decision"),
               names_to = "indicator", values_to = "values") 

classification_scenario$scenario_num <- 1:nrow(classification_scenario)
df_stats$scenario_num <- df_stats$Estimator_num + (df_stats$Population_parameters_num - 1) * 3 ## récupéré de "processing_results.R", j'espère que c'est bien ça
df_stats <- merge(df_stats, classification_scenario[, c("scenario_num", "covariate_names", "outcome_regression_model", "any_nonnormal", "bY_A_X1", "bY_A_X2")], by = "scenario_num", all.x = TRUE, all.y = TRUE) # récupération des variables incluses dans le modèle de propension ou le modèle de l'outcome
df_stats$covariate_names <- sapply(df_stats$covariate_names, paste, collapse = ", ")

df_stats$method <- NA
table(df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "MAIC" & df_stats$Anchored == "Unanchored", "MAIC unanchored", df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "MAIC" & df_stats$Anchored == "Anchored", "MAIC anchored", df_stats$method)
table(df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "STC" & df_stats$Anchored == "Unanchored", "STC unanchored", df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "STC" & df_stats$Anchored == "Anchored", "STC anchored", df_stats$method)
table(df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "Unadjusted" & df_stats$Anchored == "Unanchored", "Unadjusted unanchored", df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "Unadjusted" & df_stats$Anchored == "Anchored", "Unadjusted anchored", df_stats$method)
table(df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "ML" & df_stats$Anchored == "Unanchored", "Direct, IPTW unanchored", df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "ML" & df_stats$Anchored == "Anchored", "Direct, IPTW anchored", df_stats$method)
table(df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "GLM" & df_stats$Anchored == "Unanchored", "Direct, Regr unanchored", df_stats$method)
df_stats$method <- ifelse(df_stats$Model == "GLM" & df_stats$Anchored == "Anchored", "Direct, Regr anchored", df_stats$method)
table(df_stats$method)
# View(df_stats[is.na(df_stats$method), ])

df_stats$methodnoanch <- sub(" anchored| unanchored", "", df_stats$method)

df_stats$methodX <- NA
df_stats$methodX[df_stats$Model == "MAIC"] <- with(df_stats[df_stats$Model == "MAIC", ], paste0(method, " : ", covariate_names))
df_stats$methodX[df_stats$Model == "STC"] <- with(df_stats[df_stats$Model == "STC", ], paste0(method, " : ", outcome_regression_model))
df_stats$methodX[df_stats$Model == "ML"] <- with(df_stats[df_stats$Model == "ML", ], paste0(method, " : ", covariate_names))
df_stats$methodX[df_stats$Model == "GLM"] <- with(df_stats[df_stats$Model == "GLM", ], paste0(method, " : ", outcome_regression_model))
df_stats$methodX <- sub(" \\+ ", ", ", gsub("\\*ttt", "", sub("Y_obs ~ ", "", df_stats$methodX)))
df_stats$methodX[is.na(df_stats$methodX)] <- df_stats$method[is.na(df_stats$methodX)]
table(df_stats$methodX)

df_stats$methodXnoanch <- sub(" anchored| unanchored", "", df_stats$methodX)

ordre <- sort(unique(df_stats$method))
ordreX <- sort(unique(df_stats$methodX))
ordrenoanch <- sort(unique(df_stats$methodnoanch))
ordreXnoanch <- sort(unique(df_stats$methodXnoanch))

df_stats$method <- factor(df_stats$method, ordre)
df_stats$methodX <- factor(df_stats$methodX, ordreX)
df_stats$methodnoanch <- factor(df_stats$methodnoanch, ordrenoanch)
df_stats$methodXnoanch <- factor(df_stats$methodXnoanch, ordreXnoanch)

df_stats$nonnormal <- ifelse(df_stats$any_nonnormal, "X2 bimodal", "X2 norm")
df_stats$nonnormal <- factor(df_stats$nonnormal, c("X2 norm", "X2 bimodal"))

df_stats$effectmod <- ifelse(df_stats$bY_A_X1 == 0 & df_stats$bY_A_X2 == 0, 1, NA)
df_stats$effectmod <- ifelse(df_stats$bY_A_X1 > 0 & df_stats$bY_A_X2 == 0, 2, df_stats$effectmod)
df_stats$effectmod <- ifelse(df_stats$bY_A_X1 == 0 & df_stats$bY_A_X2 > 0, 3, df_stats$effectmod)
df_stats$effectmod <- ifelse(df_stats$bY_A_X1 > 0 & df_stats$bY_A_X2 > 0, 4, df_stats$effectmod)
df_stats$effectmod <- factor(df_stats$effectmod, 1:4, c("no TEM", "eff. mod. X1", "X2 TEM", "eff. mod. X1X2"))

df_stats$X <- sub("(^.+ : )(.+$)", "\\2", df_stats$methodX)

df_stats <- df_stats |> dplyr::mutate(X = paste0("Adjusted on ", X)) |> 
  dplyr::filter(!grepl("Unadjusted", X))

## Bias
ggplot(data = df_stats[df_stats$indicator == "bias", ], aes(x = values, y = method, color = Anchored, shape = X)) +
  geom_point(size = 3) + 
  geom_vline(xintercept = 0) +
  facet_grid(methodnoanch ~ nonnormal + effectmod, scales = "free_y") +
  scale_color_manual(values = alpha(c("#3262ab", "#de6757")), name = NULL, guide = FALSE) +
  scale_shape(name = "Confounding") +
  labs(x = "Bias", y = NULL) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(path_results, "bias.pdf"), width = 15, height = 7)
## Pourquoi 3 lignes ici ?
df_stats[df_stats$methodX == "Unadjusted anchored" & df_stats$Anchored == "Anchored" & df_stats$indicator == "bias" & df_stats$Population_parameters_num == 1, ]

## rmse
ggplot(data = df_stats[df_stats$indicator == "rmse", ], aes(x = values, y = method, color = Anchored, shape = X)) +
  geom_point(size = 3) + 
  geom_vline(xintercept = 0) +
  facet_grid(methodnoanch ~ nonnormal + effectmod, scales = "free_y") +
  scale_color_manual(values = alpha(c("#3262ab", "#de6757")), name = NULL, guide = FALSE) +
  scale_shape(name = "Confounding") +
  labs(x = "RMSE", y = NULL) +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(path_results, "rmse.pdf"), width = 15, height = 7)

## vr
ggplot(data = df_stats[df_stats$indicator == "vr", ], aes(x = values, y = method, color = Anchored, fill = Anchored, shape = X)) +
  geom_point() + # pour unajdusted, plusieurs résultats pour une mm combinaison (pas trouver pourquoi)
  geom_vline(xintercept = 1) +
  facet_grid(methodnoanch ~ Population_parameters_num + nonnormal + effectmod, scales = "free_y") +
  scale_color_manual(values = alpha(c("#3262ab", "#de6757"))) +
  theme(axis.text.x = element_text(size = 6))
  labs(x = "VR", y = NULL) + 
ggsave(file.path(path_results, "vr.pdf"), width = 15, height = 7)

## cov_95
ggplot(data = df_stats[df_stats$indicator == "cov_95", ], aes(x = values, y = method, color = Anchored, fill = Anchored, shape = X)) +
  geom_point() + # pour unajdusted, plusieurs résultats pour une mm combinaison (pas trouver pourquoi)
  geom_vline(xintercept = 0.95) +
  facet_grid(methodnoanch ~ Population_parameters_num + nonnormal + effectmod, scales = "free_y") +
  scale_color_manual(values = alpha(c("#3262ab", "#de6757"))) +
  theme(axis.text.x = element_text(size = 6))
  labs(x = "95% coverage", y = NULL) + 
ggsave(file.path(path_results, "cov_95.pdf"), width = 15, height = 7)

## Correct decision
ggplot(data = df_stats[df_stats$indicator == "correct_decision", ], aes(x = values, y = method, color = Anchored, fill = Anchored, shape = X)) +
  geom_point() + # pour unajdusted, plusieurs résultats pour une mm combinaison (pas trouver pourquoi)
  geom_vline(xintercept = 1) +
  facet_grid(methodnoanch ~ Population_parameters_num + nonnormal + effectmod, scales = "free_y") +
  scale_color_manual(values = alpha(c("#3262ab", "#de6757"))) +
  theme(axis.text.x = element_text(size = 6))
ggsave(file.path(path_results, "correct_decision.pdf"), width = 15, height = 7)

