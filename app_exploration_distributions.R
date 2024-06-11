library(shiny)
library(data.table)
library(ggplot2)
source("data_generation.R")


# Define UI
ui <- fluidPage(
  titlePanel("Covariate Distribution in Trials"),
  sidebarLayout(
    sidebarPanel(
      selectInput("imbalanced_trial", "Imbalanced Trial", choices = c("AC", "BC")),
      numericInput("bT_X1", "bT_X1", value = 10),
      numericInput("bT_X2", "bT_X2", value = -10),
      numericInput("bY_X1", "bY_X1", value = 0),
      numericInput("bY_X2", "bY_X2", value = 0),
      numericInput("bY_A_X1", "bY_A_X1", value = 0),
      numericInput("bY_A_X2", "bY_A_X1", value = 0),
      selectInput("imbalanced_trial", "Imbalanced Trial", choices = c("AC", "BC")),
      textInput("f_X1", "f_X1", value = "rnorm(N_pop, 0, 1)", placeholder ="rnorm(N_pop, 0, 1)"),
      textInput("f_X2", "f_X2", value = "rnorm(N_pop, 0, 1)", placeholder ="rnorm(N_pop, 0, 1)"),
      # selectInput("genFuncX1", "Generating Function for X1",
      #             choices = c(
      #               "bimod" = "bimod",
      #               "rnorm" = "rnorm", "runif" = "runif", "rbinom" = "rbinom")),
      # conditionalPanel(
      #   condition = "input.genFuncX1 == 'bimod'",
      #   numericInput("mean1X1", "Mean for X1 (rnorm)", value = -2),
      #   numericInput("mean2X1", "Mean2 for X1 (rnorm)", value = 2),
      #   numericInput("sdX1", "SD for X1 (rnorm)", value = 1, min = 0),
      #   numericInput("sd2X1", "SD2 for X1 (rnorm)", value = 1, min = 0)
      # ),
      # conditionalPanel(
      #   condition = "input.genFuncX1 == 'rbinom'",
      #   sliderInput("probX1", "Probability for X1 (rbinom)", min = 0, max = 1, value = 0.5)
      # ),
      # conditionalPanel(
      #   condition = "input.genFuncX1 == 'rnorm'",
      #   numericInput("meanX1", "Mean for X1 (rnorm)", value = 0),
      #   numericInput("sdX1", "SD for X1 (rnorm)", value = 1)
      # ),
      # conditionalPanel(
      #   condition = "input.genFuncX1 == 'runif'",
      #   numericInput("minX1", "Min for X1 (runif)", value = 0),
      #   numericInput("maxX1", "Max for X1 (runif)", value = 1)
      # ),
      # selectInput("genFuncX2", "Generating Function for X2",
      #             choices = c("rnorm" = "rnorm")),
      # conditionalPanel(
      #   condition = "input.genFuncX2 == 'rnorm'",
      #   numericInput("mean1X2", "Mean for X2 (rnorm)", value = -2),
      #   numericInput("sd1X2", "SD for X2 (rnorm)", value = 1, min = 0)
      # ),
      # conditionalPanel(
      #   condition = "input.genFuncX1 == 'rbinom'",
      #   sliderInput("probX1", "Probability for X1 (rbinom)", min = 0, max = 1, value = 0.5)
      # ),
      # conditionalPanel(
      #   condition = "input.genFuncX1 == 'rnorm'",
      #   numericInput("meanX1", "Mean for X1 (rnorm)", value = 0),
      #   numericInput("sdX1", "SD for X1 (rnorm)", value = 1)
      # ),
      # conditionalPanel(
      #   condition = "input.genFuncX1 == 'runif'",
      #   numericInput("minX1", "Min for X1 (runif)", value = 0),
      #   numericInput("maxX1", "Max for X1 (runif)", value = 1)
      # ),
      actionButton("runSimulation", "Run Simulation")
    ),
    mainPanel(
      plotOutput("plot_covariates_distribution"),
      plotOutput("plot_outcome_distribution")
    )
  )
)

# Define server logic
server <- function(input, output) {

  bindEvent(input$runSimulation,
            x = observe({
              N_RCT <- 10000 # Example value, adjust based on your simulation setup
              N_pop <- 10^5
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
                bY_A = 1.5,  # Effet de A par rapport à C
                bY_B = 1.5,  # Effet de B par rapport à C
                bY_C = 0,    # Pas d'effet de C sur l'outcome
                imbalanced_trial = c(input$imbalanced_trial),
                # imbalanced_trial = c("AC"),
                imbalanced_trial_model = bquote(X1 * bT_X1 + X2 * bT_X2 ), # Modèle d'attribution de l'essai AC
                balanced_trial_model = bquote(0),  # Modèle d'attribution de l'essai BC
                outcome_distribution = "normal",
                outcome_generation_formula =  bquote(
                  bY_X1*X1 + bY_X2 * X2 +  (bY_A + bY_A_X1*X1 + bY_A_X2*X2) * A +  bY_B*B + bY_C*C
                )
              )

              populations <- creating_population(list_simulation_parameters)
              pop_init <- populations$pop_init
              selected_individuals_AC <- pop_init[sample(id, N_RCT, replace = TRUE, prob = prob_w_trial_AC)][
                , ttt := rep_len(c("A", "C"), length.out = .N)]
              selected_individuals_BC <- pop_init[sample(id, N_RCT, replace = TRUE, prob = prob_w_trial_BC)][
                , ttt := rep_len(c("B", "C"), length.out = .N)]

              all_individuals <- data.table::rbindlist(list("AC" = selected_individuals_AC, "BC" = selected_individuals_BC),
                                                       idcol = "trial") |>
                dplyr::mutate(across(tidyselect::matches("X[0-9]+"), as.double)) |>
                data.table::melt(measure.vars = patterns("X[0-9]+"), value.name = "variable_value", number = as.numerical) |>
                # Because they are the only variables used for now
                dplyr::filter(variable %in% c("X1", "X2")) |>
                dplyr::mutate(trial = ifelse(trial == "AC", "AC (IPD)", trial))


              plot_covariates_distribution <- all_individuals |>
                ggplot() +
                geom_density(aes(variable_value, fill = trial), alpha = 0.4) +
                facet_wrap(facets = "variable") +
                labs(x = NULL, y = NULL, title = "Covariates distributions") +
                theme(strip.text = element_text(size = 12))


              if (list_simulation_parameters$outcome_distribution == "normal") {
                plot_outcome_distribution <- all_individuals |>
                  ggplot() +
                  # geom_density(aes(Y_obs, fill = ttt), alpha = 0.4) +
                  geom_violin(aes(ttt, Y_obs, fill = ttt), alpha = 0.4) +
                  geom_pointrange(aes(y = mean_Y_obs, x = ttt, ymin = low, ymax = up), color = "black", size = 1,
                                  data = all_individuals |>
                                    dplyr::group_by(trial, ttt) |>
                                    dplyr::summarise(mean_Y_obs = mean(Y_obs), sd_Y_obs = sd(Y_obs), low = mean_Y_obs - sd_Y_obs, up = mean_Y_obs + sd_Y_obs)) +
                  facet_wrap(facets = "trial", ncol = 2) +
                  labs(x = NULL, y = NULL, title = "Outcome distribution")
                print(plot_outcome_distribution)

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
                  facet_wrap(facets = "trial", ncol = 1) +
                  # geom_bar(aes(Y_obs, fill = trial), alpha = 0.4) +
                  labs(x = NULL, y = NULL, title = "Outcome distribution")
              }
              output$plot_covariates_distribution <- renderPlot(plot_covariates_distribution)
              output$plot_outcome_distribution <- renderPlot(plot_outcome_distribution)



              # binary_marker = rbinom(N_pop, 1, 0.5)
              # Dynamically generate X1 based on user input
              # generateX1 <- switch(input$genFuncX1,
              #                      "rbinom" = rbinom(N_pop, 1, input$probX1),
              #                      "rnorm" = rnorm(N_pop, input$meanX1, input$sdX1),
              #                      "runif" = runif(N_pop, input$minX1, input$maxX1),
              #                      "bimod" = binary_marker * rnorm(N_pop, input$mean1X1, input$sdX1) +
              #                        (1 - binary_marker) * rnorm(N_pop, input$mean2X1, input$sd2X1)
              # )
              # generateX2 <- switch(input$genFuncX2,
              #                      # "rbinom" = rbinom(N_pop, 1, input$probX1),
              #                      "rnorm" = rnorm(N_pop, input$mean1X2, input$sd1X2),
              #                      # "runif" = runif(N_pop, input$minX1, input$maxX1),
    #                      # "bimod" = binary_marker * rnorm(N_pop, input$mean1X1, input$sdX1) +
    #                      #   (1 - binary_marker) * rnorm(N_pop, input$mean2X1, input$sd2X1)
    # )
    # # Create population with the dynamic X1
    # pop_init <- data.table(
    #   id = 1:N_pop,
    #   bT_X1 = input$bT_X1,   # Effet de la variable binaire sur la probabilité d'être dans l'essai AC
    #   bT_X2 = input$bT_X2,   # Effet de la variable continue X2...
    #   bT_X3 = 0,     # Idem, mais inutile pour le moment
    #   bT_X4 = 0,     # Idem, mais inutile pour le moment
    #   X1 = generateX1,
    #   X2 = generateX2
    #   # Extend for other variables as needed
    # )
    #
    # # AC_trial_model = c(bquote(X1 * bT_X1 + X2 * bT_X2 + X3 * bT_X3 + X4 * bT_X4)) # Modèle d'attribution de l'essai AC
    # imbalanced_trial <- input$imbalanced_trial
    # imbalanced_model <- bquote(X1 * bT_X1 + X2 * bT_X2) # Modèle d'attribution de l'essai AC
    # if (imbalanced_trial == "AC") {
    #   BC_trial_model = bquote(0)  # Modèle d'attribution de l'essai BC
    #   AC_trial_model = imbalanced_model
    # } else {
    #   BC_trial_model = imbalanced_model
    #   AC_trial_model = bquote(0)
    # }
    #
    #
    # trial_assignement_prob <- function(trial_assignment_model, df) {
    #   predicted <- with(df, eval(trial_assignment_model))
    #   return(plogis(predicted))
    # }
    #
    # pop_init[, prob_w_trial_AC := trial_assignement_prob(AC_trial_model, df = pop_init)]
    # pop_init[, prob_w_trial_BC := trial_assignement_prob(BC_trial_model, df = pop_init)]
    #
    #
    # # Draw trials based on your code
    # selected_individuals_AC <- pop_init[sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_w_trial_AC)][
    #   , ttt := rep_len(c("A", "C"), length.out = .N)]
    # trial_AC <- selected_individuals_AC
    #
    # selected_individuals_BC <- pop_init[sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_w_trial_BC)][
    #   , ttt := rep_len(c("B", "C"), length.out = .N)]
    # trial_BC <- selected_individuals_BC
    #
    # df_trials <- rbindlist(list("AC" = trial_AC, "BC" = trial_BC), idcol = "trials")
              # return(df_trials)
            })
  )

  # output$plotX1 <- renderPlot({
  #   ggplot(simulationOutput(), aes(x = X1, color = trials)) +
  #     # geom_histogram(binwidth = 0.1) +
  #     geom_density() +
  #     ggtitle("Covariate Distribution in trials") +
  #     xlab("X1") +
  #     ylab("Frequency")
  # })
  #
  # output$plotX2 <- renderPlot({
  #   ggplot(simulationOutput(), aes(x = X2, color = trials)) +
  #     # geom_histogram(binwidth = 0.1) +
  #     geom_density() +
  #     ggtitle("Covariate Distribution in trials") +
  #     xlab("X2") +
  #     ylab("Frequency")
  # })

}

# Run the app
shinyApp(ui = ui, server = server)
