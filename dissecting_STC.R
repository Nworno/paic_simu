# Proof that STC is merely estimating the effect using half of the data as compared to a standard regression model with 
# access to full IPD.

# Simulate some data
set.seed(123)
n <- 1000

# Individual data
x1 <- rnorm(n)
x2 <- rnorm(n)
# y <- 2*x1 + 3*x2 + rnorm(n)
# df <- data.frame(y, x1, x2)
# (average_effect_IndD <- mean(df$y))
ttt <- rep(c("A", "C"), n/2)
y <- 1 + 2*x1 + 3*x2 + 3*(ttt == "A") + x1*(ttt == "A") + rnorm(n)
df <- data.frame(y, x1, x2, ttt)
(average_effect_IndD <- mean(df$y[df$ttt == "A"]) - mean(df$y[df$ttt == "C"]))

# Aggregated data
x1 <- rnorm(n, 1)
x2 <- rnorm(n, 1)
# y <- 1 + 2*x1 + 3*x2 + rnorm(n)
# df2 <- data.frame(y, x1, x2)
# (average_effect_AgD <- mean(df2$y))
ttt <- rep(c("B", "C"), n/2)
y <- 1 + 2*x1 + 3*x2 + 2*(ttt == "B") + rnorm(n)
df2 <- data.frame(y, x1, x2, ttt = as.factor(ttt))
(average_effect_AgD <- mean(df2$y[df2$ttt == "B"]) - mean(df2$y[df2$ttt == "C"]))


# Unadjusted effect
mean(df$y) - mean(df2$y)

# Bucher ITC 
average_effect_IndD - average_effect_AgD

# Fit the standard regression model and estimate the ttt effect
# fitted_model <- lm(y ~ x1 + x2 + ttt, data = df)
# fitted_model$coefficients[["ttt1"]]
# fitted_model_2 <- lm(y ~ x1 + x2 + ttt, data = df2)
# fitted_model_2$coefficients[["ttt1"]]
# relative_ttt_effect <- fitted_model$coefficients[["ttt1"]] - fitted_model_2$coefficients[["ttt1"]]

# Using FULL IPD DATA: that is fitting a
## Unanchored 
fitted_model <- lm(y ~ x1 + x2, data = df)
fitted_model$coefficients[["(Intercept)"]]
fitted_model_2 <- lm(y ~ x1 + x2, data = df2)
fitted_model_2$coefficients[["(Intercept)"]]
(relative_ttt_effect <- fitted_model$coefficients[["(Intercept)"]] - fitted_model_2$coefficients[["(Intercept)"]])
## Unanchored one step 
fitted_model <- lm(y ~ x1 + x2 + ttt, data = dplyr::bind_rows("A" = df, "B" = df2, .id = "ttt"))
fitted_model$coefficients[["tttB"]]

## Anchored 
df_anchored <- dplyr::bind_rows("AC" = df, "BC" = df2, .id = "trial")
df_anchored$ttt <- relevel(as.factor(df_anchored$ttt), ref = "B")
df_anchored$x1 <- df_anchored$x1 - mean(df2$x1)
fitted_model <- lm(y ~ trial + x1:ttt + ttt, data = df_anchored)
fitted_model$coefficients

# Estimate the ttt effect by predicting the average effect in AgD with ttt and the average effect in AgD without ttt
predicted_average_effect_in_AgD_w_ttt <- predict(fitted_model, lapply(df2, mean)) |> mean()
predicted_average_effect_AgD <- predict(fitted_model_2, df2) |> mean()
(relative_ttt_effect <- predicted_average_effect_in_AgD_w_ttt - predicted_average_effect_AgD)
(relative_ttt_effect <- predicted_average_effect_in_AgD_w_ttt - average_effect_AgD)

# Estimate the ttt effect by fitting a model in the AgD while setting the coefficients of the AgD model to the coefficients of the IPD model
coefficients_IPD <- coef(fitted_model)
intercept <- coefficients_IPD[["(Intercept)"]]
coef_x1 <- coefficients_IPD[["x1"]]
coef_x2 <- coefficients_IPD[["x2"]]
fitted_model_2 <- lm(y ~ offset(coef_x1 * x1) + offset(coef_x2 * x2), data = df2)
(ttt_effect <- fitted_model_2$coefficients[["(Intercept)"]] - intercept)

# Estimate the ttt effect by predicting the ttt effect when combining both IPD and AgD, but fixing the coefficients to
# their estimation in the IPD model
new_model <- lm(y ~ offset(coef_x1*x1) + offset(coef_x2 * x2) + ttt, data = dplyr::bind_rows("A" = df, "B" = df2, .id = "ttt"))
new_model$coefficients[["tttB"]] 


# Fit the STC model as described in the princeps paper, that is by first centering the variables to the 
df_centered <- sweep(df, STATS = colMeans(df2), FUN = "-", MARGIN = 2)
df_centered[, "y"] <- df[, "y"]
fitted_model_stc <- lm(y ~ x1 + x2, data = df_centered)
(relative_effect <- coef(fitted_model_stc)[["(Intercept)"]] - mean(df2$y))
