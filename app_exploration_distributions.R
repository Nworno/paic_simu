library(shiny)
library(data.table)
library(ggplot2)


N_RCT <- 10000 # Example value, adjust based on your simulation setup
N_pop <- 10^5

# Define UI
ui <- fluidPage(
  titlePanel("Covariate Distribution in Trials"),
  sidebarLayout(
    sidebarPanel(
      selectInput("imbalanced_trial", "Imbalanced Trial", choices = c("AC", "BC")),
      numericInput("bT_X1", "bT_X1", value = 1), 
      numericInput("bT_X2", "bT_X2", value = 1), 
      selectInput("genFuncX1", "Generating Function for X1",
                  choices = c(
                    "bimod" = "bimod",
                    "rnorm" = "rnorm", "runif" = "runif", "rbinom" = "rbinom")),
      conditionalPanel(
        condition = "input.genFuncX1 == 'bimod'",
        numericInput("mean1X1", "Mean for X1 (rnorm)", value = -2),
        numericInput("mean2X1", "Mean2 for X1 (rnorm)", value = 2),
        numericInput("sdX1", "SD for X1 (rnorm)", value = 1, min = 0),
        numericInput("sd2X1", "SD2 for X1 (rnorm)", value = 1, min = 0)
      ),
      conditionalPanel(
        condition = "input.genFuncX1 == 'rbinom'",
        sliderInput("probX1", "Probability for X1 (rbinom)", min = 0, max = 1, value = 0.5)
      ),
      conditionalPanel(
        condition = "input.genFuncX1 == 'rnorm'",
        numericInput("meanX1", "Mean for X1 (rnorm)", value = 0),
        numericInput("sdX1", "SD for X1 (rnorm)", value = 1)
      ),
      conditionalPanel(
        condition = "input.genFuncX1 == 'runif'",
        numericInput("minX1", "Min for X1 (runif)", value = 0),
        numericInput("maxX1", "Max for X1 (runif)", value = 1)
      ),
      selectInput("genFuncX2", "Generating Function for X2",
                  choices = c("rnorm" = "rnorm")),
      conditionalPanel(
        condition = "input.genFuncX2 == 'rnorm'",
        numericInput("mean1X2", "Mean for X2 (rnorm)", value = -2),
        numericInput("sd1X2", "SD for X2 (rnorm)", value = 1, min = 0)
      ),
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
      plotOutput("plotX1"),
      plotOutput("plotX2")
    )
  )
)

# Define server logic
server <- function(input, output) {
  simulationOutput <- eventReactive(input$runSimulation, {
    
    binary_marker = rbinom(N_pop, 1, 0.5)
    # Dynamically generate X1 based on user input
    generateX1 <- switch(input$genFuncX1,
                         "rbinom" = rbinom(N_pop, 1, input$probX1),
                         "rnorm" = rnorm(N_pop, input$meanX1, input$sdX1),
                         "runif" = runif(N_pop, input$minX1, input$maxX1), 
                         "bimod" = binary_marker * rnorm(N_pop, input$mean1X1, input$sdX1) +
                           (1 - binary_marker) * rnorm(N_pop, input$mean2X1, input$sd2X1)
    )
    generateX2 <- switch(input$genFuncX2,
                         # "rbinom" = rbinom(N_pop, 1, input$probX1),
                         "rnorm" = rnorm(N_pop, input$mean1X2, input$sd1X2),
                         # "runif" = runif(N_pop, input$minX1, input$maxX1), 
                         # "bimod" = binary_marker * rnorm(N_pop, input$mean1X1, input$sdX1) +
                         #   (1 - binary_marker) * rnorm(N_pop, input$mean2X1, input$sd2X1)
    )
    
    
    # Create population with the dynamic X1
    pop_init <- data.table(
      id = 1:N_pop,
      bT_X1 = input$bT_X1,   # Effet de la variable binaire sur la probabilité d'être dans l'essai AC
      bT_X2 = input$bT_X2,   # Effet de la variable continue X2...
      bT_X3 = 0,     # Idem, mais inutile pour le moment
      bT_X4 = 0,     # Idem, mais inutile pour le moment
      X1 = generateX1, 
      X2 = generateX2
      # Extend for other variables as needed
    )
    
    # AC_trial_model = c(bquote(X1 * bT_X1 + X2 * bT_X2 + X3 * bT_X3 + X4 * bT_X4)) # Modèle d'attribution de l'essai AC
    imbalanced_trial <- input$imbalanced_trial
    imbalanced_model <- bquote(X1 * bT_X1 + X2 * bT_X2) # Modèle d'attribution de l'essai AC
    if (imbalanced_trial == "AC") {
      BC_trial_model = bquote(0)  # Modèle d'attribution de l'essai BC
      AC_trial_model = imbalanced_model
    } else {
      BC_trial_model = imbalanced_model
      AC_trial_model = bquote(0)
    }
    
      
    trial_assignement_prob <- function(trial_assignment_model, df) {
      predicted <- with(df, eval(trial_assignment_model))
      return(plogis(predicted))
    }
    
    pop_init[, prob_w_trial_AC := trial_assignement_prob(AC_trial_model, df = pop_init)]
    pop_init[, prob_w_trial_BC := trial_assignement_prob(BC_trial_model, df = pop_init)]
    
    
    # Draw trials based on your code
    selected_individuals_AC <- pop_init[sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_w_trial_AC)][
      , ttt := rep_len(c("A", "C"), length.out = .N)]
    trial_AC <- selected_individuals_AC
    
    selected_individuals_BC <- pop_init[sample(pop_init$id, N_RCT, replace = FALSE, prob = pop_init$prob_w_trial_BC)][
      , ttt := rep_len(c("B", "C"), length.out = .N)]
    trial_BC <- selected_individuals_BC
    
    df_trials <- rbindlist(list("AC" = trial_AC, "BC" = trial_BC), idcol = "trials")
    return(df_trials)
  })
  
  output$plotX1 <- renderPlot({
    ggplot(simulationOutput(), aes(x = X1, color = trials)) + 
      # geom_histogram(binwidth = 0.1) + 
      geom_density() + 
      ggtitle("Covariate Distribution in trials") +
      xlab("X1") +
      ylab("Frequency")
  })
  
  output$plotX2 <- renderPlot({
    ggplot(simulationOutput(), aes(x = X2, color = trials)) + 
      # geom_histogram(binwidth = 0.1) + 
      geom_density() + 
      ggtitle("Covariate Distribution in trials") +
      xlab("X2") +
      ylab("Frequency")
  })
  
}

# Run the app
shinyApp(ui = ui, server = server)
