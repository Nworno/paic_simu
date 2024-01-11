library(shiny)
library(ggplot2)
library(data.table)
library(dplyr)
library(shinydashboard)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("flat")
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("fresh")
theme_set(theme_bw() + theme(text = element_text(size = 16)))

ui <- dashboardPage(
  dashboardHeader(title = "PAIC Results Exploration"),
  dashboardSidebar(collapsed = TRUE),
  dashboardBody(
    selectInput("num_experiment", "Number of Experiment", choices = 1:8, selected = 2),
    fluidRow(
      column(8, 
             plotOutput("resultsPlot1"),
             plotOutput("resultsPlot2"),
             plotOutput("resultsPlot3"),
      ),
      column(4,
             plotOutput("treatment_effect"),
             plotOutput("covariates_distribution_Output")
      )
    )
  )
)


server <- function(input, output) {
  date_simulations <- "20240109_122850"
  dir_simulations <- file.path("results_simulations", date_simulations)
  
  df_population_parameters <- readRDS(file.path(dir_simulations, "df_population_parameters.RDS")) 
  
  path_results_experiment <- reactive({
    file.path(file.path(dir_simulations, input$num_experiment))
  })
    

  true_effect <- reactive({
    true_effectS <- readRDS(file.path(path_results_experiment(), "/average_outcome_df.RDS"))
    true_effectS[outcome_type == "conditional", AB, drop = TRUE]
  })
  combined_results <- reactive({
    list_files <- list.files(path_results_experiment(), pattern = "experiment_.*.RDS", full.names = TRUE)
    results_experiments <- lapply(list_files, function(file) readRDS(file))
    results_experiments |> 
      lapply(dplyr::bind_rows) |>
      dplyr::bind_rows(.id = "estimator_num") |> 
      # dplyr::group_by(estimator_num) |>
      dplyr::mutate(mean_estimate = mean(estimate),
                    lb = estimate - qnorm(0.975)*sqrt(variance),
                    ub = estimate + qnorm(0.975)*sqrt(variance),
                    includes_true_effect = (lb <= true_effect()) & (ub >= true_effect()),
                    includes_0 = (lb <= 0) & (ub >= 0),
                    correct_decision = (true_effect() == 0 & includes_0) | 
                      ((true_effect() != 0) & (sign(true_effect()) == sign(estimate)) & !includes_0)
      )
  })

  output$covariates_distribution_Output <- renderPlot({
    readRDS(file.path(path_results_experiment(), "covariates_distribution.RDS"))
  })
  
  # output$covariates_distribution_Output <- renderImage({
  #   # Display the png file as a plot
  #   list(src = file.path(path_results_experiment(), "covariates_distribution.png"),
  #        contentType = 'image/png',
  #        width = 1000,
  #        height = 500,
  #        alt = "This is alternate text") 
  #   # covariates_distribution()
  # }, deleteFile = TRUE)
  
  output$treatment_effect <- renderPlot({
    true_treatment_effect <- readRDS(file.path(dir_simulations, input$num_experiment, "average_outcome_df.RDS"))[outcome_type == "conditional", AB]
    
    df_population_parameters |>
      dplyr::filter(population_parameters_num == input$num_experiment) |>
      # dplyr::filter(population_parameters_num == 1) |> 
      dplyr::select(matches("bY.+X[12]")) |> 
      dplyr::mutate(conditional_AB_effect = true_treatment_effect) |> 
      dplyr::mutate(across(everything(), as.double)) |> 
      tidyr::pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") |> 
      ggplot() +
      geom_point(aes(y = parameter, x = value, color = parameter), size = 3) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
      geom_vline(xintercept = true_treatment_effect, linetype = "dashed", colour = "red") +
      guides(color = "none")
  })
  
  output$plots <- renderUI({
    combined_results_df <- combined_results()
    plot_output_list <- lapply(names(list_results_plots), function(name) {
      plotname <- paste0("Estimator ", name)
      plotOutput(plotname)
    })
    do.call(tagList, plot_output_list)
  })
  
  renderPlotForEstimatorNum <- function(chosen_estimator_num) {
    renderPlot({
      true_effect <- readRDS(file.path(dir_simulations, input$num_experiment, "average_outcome_df.RDS"))[outcome_type == "conditional", AB]
      result_plot <- combined_results() |> 
        dplyr::filter(estimator_num == chosen_estimator_num) |> 
        ggplot() +
        geom_rect(aes(xmin = -Inf, xmax = true_effect, ymin = true_effect, ymax = Inf), fill = "#A0D2AD", alpha = 0.02) +
        geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), fill = "grey", alpha = 0.01) +
        geom_point(aes(x = lb, y = ub, color = anchored), alpha = 0.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
        facet_wrap(~adjustment + model)
      
      return(result_plot)
    })
  }
  
  output$resultsPlot1 <- renderPlotForEstimatorNum("1")
  output$resultsPlot2 <- renderPlotForEstimatorNum("2")
  output$resultsPlot3 <- renderPlotForEstimatorNum("3")
    
}

shinyApp(ui = ui, server = server)
