library(shiny)
library(data.table)
library(ggplot2)
source("data_generation.R")


# Define UI
ui <- fluidPage(
  titlePanel("Covariate Distribution in Trials"),
  sidebarLayout(
    sidebarPanel(
      actionButton("runSimulation", "Run Simulation"),
      selectInput("outcome_distribution", "Outcome Distribution", choices = c("normal", "binomial")),
      selectInput("with_replace", "Sampling with replacement", choices = c(TRUE, FALSE)),
      selectInput("imbalanced_trial", "Imbalanced Trial", choices = c("AC", "BC")),
      numericInput("bY_A", "bY_A", value = 1.5),
      numericInput("bY_B", "bY_B", value = 1.5),
      numericInput("bY_C", "bY_C", value = 0),
      numericInput("bT_X1", "bT_X1", value = 0),
      numericInput("bT_X2", "bT_X2", value = 0 ),
      numericInput("bY_X1", "bY_X1", value = 0),
      numericInput("bY_X2", "bY_X2", value = 0),
      numericInput("bY_A_X1", "bY_A_X1", value = 0),
      numericInput("bY_A_X2", "bY_A_X2", value = 0),
      textInput("f_X1", "f_X1", value = "rnorm(N_pop, 0, 1)", placeholder ="rnorm(N_pop, 0, 1)"),
      textInput("f_X2", "f_X2", value = "rnorm(N_pop, 0, 1)", placeholder ="rnorm(N_pop, 0, 1)"),
      actionButton("runSimulation", "Run Simulation")
    ),
    mainPanel(
      plotOutput("plot_propensity_distribution"),
      plotOutput("plot_covariates_distribution"),
      plotOutput("plot_outcome_distribution"),
      tableOutput("marginal_effect_table")
    )
  )
)


# Define server logic
server <- function(input, output) {

  theme_set(theme_bw())
  theme_update(text = element_text(size = 20))

  bindEvent(input$runSimulation,
            x = observe({
              N_RCT <- 5*10^4 # Example value, adjust based on your simulation setup
              N_pop <- 5*10^5
              print(N_pop)
              print(N_RCT)
              list_simulation_parameters <- list(
                N_pop = 10^5,
                N_RCT = 10000,
                bT_X1 = input$bT_X1,   # Effet de la variable sur la probabilité d'être dans l'essai BC
                bT_X2 = input$bT_X2,   # Effet de la variable continue X2...
                bY_X1 = input$bY_X1,   # Effet de X1 sur l'outcome
                bY_X2 = input$bY_X2,   # Effet de X2 sur l'outcome
                bY_A_X1 = input$bY_A_X1, # Interaction A et X1 dans le modèle outcome
                bY_A_X2 = input$bY_A_X2, # Interaction A et X2 dans le modèle outcome
                binary_marker = bquote(rbinom(N_pop, 1, 0.5)), # Utilisé pour la variable bimodale
                f_X1 = parse(text = input$f_X1),
                f_X2 = parse(text = input$f_X2), # distribution de X2
                bY_A = input$bY_A,  # Effet de A par rapport à C
                bY_B = input$bY_B,  # Effet de B par rapport à C
                bY_C = input$bY_C,    # Pas d'effet de C sur l'outcome
                imbalanced_trial = c(input$imbalanced_trial),
                # imbalanced_trial = c("AC"),
                imbalanced_trial_model = bquote(X1 * bT_X1 + X2 * bT_X2 ), # Modèle d'attribution de l'essai AC
                balanced_trial_model = bquote(0),  # Modèle d'attribution de l'essai BC
                outcome_distribution = input$outcome_distribution,
                outcome_generation_formula =  bquote(
                  bY_X1*X1 + bY_X2 * X2 +  (bY_A + bY_A_X1*X1 + bY_A_X2*X2) * A +  bY_B*B + bY_C*C
                )
              )

              populations <- creating_population(list_simulation_parameters)
              pop_init <- populations$pop_init


              # Should subset the trials here based on whether anchored or not
              stopifnot(levels(pop_init$ttt)[[1]] == "C")
              all_individuals <- pop_init |>
                data.table::melt(measure.vars = patterns("X[0-9]+"), value.name = "variable_value", number = as.numerical) |>
                # Because they are the only variables used for now
                dplyr::filter(variable %in% c("X1", "X2")) |>
                data.table::melt(measure.vars = c("A", "B", "C"), value.name = "Y_obs", variable.name = "ttt", number = as.numerical) |>
                dplyr::mutate(trial = factor(trial, levels = c("AC", "BC"), labels = c("AC (IPD)", "BC (AgD)")))

              marginal_effect <- all_individuals |>
                dplyr::select(trial, ttt, Y_obs) |>
                dplyr::group_by(trial, ttt) |>
                dplyr::summarize(Y_obs = mean(Y_obs)) |>
                tidyr::pivot_wider(names_from = "ttt", values_from = "Y_obs") |>
                dplyr::mutate(AB = A - B)

              plot_covariates_distribution <- all_individuals |>
                ggplot() +
                geom_density(aes(variable_value, fill = trial), alpha = 0.4) +
                facet_wrap(facets = "variable") +
                labs(x = NULL, y = NULL, title = "Covariates distributions") +
                theme(strip.text = element_text(size = 12))

              plot_propensity_distribution <- all_individuals |>
                ggplot() +
                geom_density(aes(prob_imbalanced_trial, fill = trial), alpha = 0.3) +
                geom_density(aes(prob_imbalanced_trial), color = "black") + # Both trials together
                labs(x = NULL, y = NULL, title = "Propensity distributions")

              if (list_simulation_parameters$outcome_distribution == "normal") {
                plot_outcome_distribution <- all_individuals |>
                  ggplot() +
                  # geom_density(aes(Y_obs, fill = ttt), alpha = 0.4) +
                  geom_violin(aes(ttt, Y_obs, fill = ttt), alpha = 0.4) +
                  geom_pointrange(aes(y = mean_Y_obs, x = ttt, ymin = low, ymax = up), color = "black", size = 1,
                                  data = all_individuals |>
                                    dplyr::group_by(trial, ttt) |>
                                    dplyr::summarise(mean_Y_obs = mean(Y_obs), sd_Y_obs = sd(Y_obs), low = mean_Y_obs - sd_Y_obs, up = mean_Y_obs + sd_Y_obs)) +
                  facet_wrap(facets = "trial", ncol = 2, scales = "free_x") +
                  labs(x = NULL, y = NULL, title = "Outcome distribution")

              } else if (list_simulation_parameters$outcome_distribution == "binomial") {
                # plot the distribution of the outcome as a barplot, with proportions of the outcome as stack bars for each treatment group
                plot_outcome_distribution <- all_individuals |>
                  dplyr::mutate(Y_obs = ifelse(Y_obs > 0, 1, 0), fill = ttt) |>
                  dplyr::group_by(trial, ttt) |>
                  dplyr::summarize("0" = 1L - mean(Y_obs), "1" = mean(Y_obs)) |>
                  tidyr::pivot_longer(cols = c("0", "1"), names_to = "prop_Y_obs", values_to = "value") |>
                  # dplyr::count(Y_obs) |>
                  # dplyr::summarize(Y_obs = dplyr::count(Y_obs), .by = c("trial", "ttt")) |>
                  ggplot() +
                  geom_bar(aes(y = value, x = ttt, fill = prop_Y_obs), stat = "identity", alpha = 0.4) +
                  # geom_bar(aes(Y_obs, position = "dodge", alpha = 0.4) +
                  facet_wrap(facets = "trial", ncol = 1, scales = "free_x") +
                  # geom_bar(aes(Y_obs, fill = trial), alpha = 0.4) +
                  labs(x = NULL, y = NULL, title = "Outcome distribution")
              }
              output$plot_covariates_distribution <- renderPlot(plot_covariates_distribution)
              output$plot_propensity_distribution <- renderPlot(plot_propensity_distribution)
              output$plot_outcome_distribution <- renderPlot(plot_outcome_distribution)
              output$marginal_effect_table <- renderTable(marginal_effect)

            })
  )
}

# Run the app
shinyApp(ui = ui, server = server)
