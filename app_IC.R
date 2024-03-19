library(shiny)
library(ggplot2)
library(data.table)
library(dplyr)
library(bs4Dash)
library(DT)
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("flat")
if ("ggthemr" %in% dimnames(installed.packages())[[1]]) ggthemr::ggthemr("fresh")
theme_set(theme_bw() + theme(text = element_text(size = 16)))
source("env_variables.R")

ui <- bs4DashPage(
  header = bs4DashNavbar(title = "PAIC Results Exploration"),
  sidebar = bs4DashSidebar(
    side = "left",
    status = "primary",
    skin = "light",
    collapsed = TRUE, 
    bs4Dash::sidebarUserPanel(
      name = tags$span("Menu"),
    ),
    bs4Dash::sidebarMenu(
      bs4Dash::menuItem(
        text = "test",
        tabName = "results",
        icon = icon("cog")),
    flat = FALSE, compact = TRUE, childIndent = TRUE
    ),
  ),
  body = bs4DashBody(
    bs4TabItems(
      bs4TabItem(tabName = "results",
                 bs4Card(
                   uiOutput("dynamic_selectors")
                 ),
                 bs4Card(
                   DT::dataTableOutput("selected_parameters"),
                   width = 12
                 ),
                 bs4Card(width = 12, plotOutput("resultsPlot1")),
                 bs4Card(width = 12, plotOutput("resultsPlot2")),
                 bs4Card(width = 12, plotOutput("resultsPlot3")),
                 bs4Card(width = 12, plotOutput("treatment_effect")),
                 bs4Card(width = 12, plotOutput("covariates_distribution_Output"))
      )
    )
  ),
  controlbar = bs4DashControlbar(),
  footer = bs4DashFooter()
)




server <- function(input, output) {

  date_simulations <- DATE_EXPERIMENT
  dir_simulations <- file.path("results_simulations", date_simulations)
  
  df_population_parameters <- readRDS(file.path(dir_simulations, "df_population_parameters.RDS"))
  df_estimators_parameters <- readRDS(file.path(dir_simulations, "df_estimators_parameters.RDS"))

  output$dynamic_selectors <- renderUI({
    names_parameters <- names(df_population_parameters)[names(df_population_parameters) != "population_parameters_num"]
    ui_elements <- lapply(names_parameters, function(column) {
      if (length(unique(df_population_parameters[[column]])) > 1) {
        selectInput(
          inputId = paste0("select_", column),
          label = column,
          choices = as.character(unique(df_population_parameters[[column]])),
          selected = as.character(unique(df_population_parameters[[column]]))[1]
        )
      }
    })
    go_button <- actionButton("go", "Go")
    do.call(tagList, list(ui_elements, go_button))
  })
  
  # num_experiment <- reactive({
  #   unique(selected_row()[, population_parameters_num])
  # })
  # 
  # 
  
  # observe({
  #   input$go
  #   print("oh yeah")
  # })

  selected_row <- bindEvent(
    input$go,
    x = reactive({
      print("oh yeah")
      filtered_df <- copy(df_population_parameters)
      for (parameter_input in names(input)[startsWith(names(input), "select_")]) {
        colname <- sub("select_", "", parameter_input, fixed = TRUE)
        filtered_df <- filtered_df[filtered_df[[colname]] == input[[parameter_input]], ]
      }
      print(filtered_df)
      if(nrow(filtered_df) == 1) {
        print("oh yeah")
        filtered_df
      } else {
        NULL # In case no row matches or multiple rows match, though your setup should prevent the latter
      }
    })
  )

  output$selected_parameters <- DT::renderDataTable({
    selected_row()[, lapply(.SD, as.character)] |> DT::datatable(options = list(scrollX = TRUE))
  })
  num_experiment <- bindEvent(
    selected_row(), 
    x = reactive({
      selected_row()$population_parameters_num
    })
  )

  # output$num_experiment_picker <- renderUI({
  #   selectInput("num_experiment",
  #               "Experiment number",
  #               choices = 1:nrow(df_population_parameters),
  #               selected = 1)
  # })
  # bindEvent(input$num_experiment, x = reactive({print(input$num_experiment)}))

  # 
  path_results_experiment <- reactive({
    req(num_experiment())
    print(num_experiment())
    file.path(file.path(dir_simulations, num_experiment()))
  })
  # 
  # 
  true_effect <- bindEvent(
    path_results_experiment(),
    x = reactive({
      true_effectS <- readRDS(file.path(path_results_experiment(), "/average_outcome_df.RDS"))
      true_effectS[outcome_type == "conditional", AB, drop = TRUE]
    })
  )

  combined_results <- bindEvent(
    path_results_experiment(),
    x = reactive({
      list_files <- list.files(path_results_experiment(), pattern = "experiment_.*.RDS", full.names = TRUE)
      results_experiments <- lapply(list_files, function(file) readRDS(file))
      results_experiments |>
        lapply(dplyr::bind_rows) |>
        dplyr::bind_rows(.id = "estimator_num") |>
        # dplyr::group_by(estimator_num) |>
        dplyr::mutate(lb = estimate - qnorm(0.975)*sqrt(variance),
                      ub = estimate + qnorm(0.975)*sqrt(variance),
                      includes_true_effect = (lb <= true_effect()) & (ub >= true_effect()),
                      includes_0 = (lb <= 0) & (ub >= 0),
                      correct_decision = (true_effect() == 0 & includes_0) |
                        ((true_effect() != 0) & (sign(true_effect()) == sign(estimate)) & !includes_0)
        )
    })
  )
  # 
  output$covariates_distribution_Output <- bindEvent(
    path_results_experiment(),
    x = renderPlot({
      readRDS(file.path(path_results_experiment(), "covariates_distribution.RDS"))
    })
  )
  # output$covariates_distribution_Output <- renderImage({
  #   # Display the png file as a plot
  #   list(src = file.path(path_results_experiment(), "covariates_distribution.png"),
  #        contentType = 'image/png',
  #        width = 1000,
  #        height = 500,
  #        alt = "This is alternate text")
  #   # covariates_distribution()
  # }, deleteFile = TRUE)

  output$treatment_effect <- bindEvent(
    path_results_experiment(),
    x = renderPlot({
      true_treatment_effect <- readRDS(file.path(path_results_experiment(), "average_outcome_df.RDS"))[outcome_type == "conditional", AB]

      selected_row() |>
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
  )


  renderPlotForEstimatorNum <- function(chosen_estimator_num) {
    covariate_names <- df_estimators_parameters[chosen_estimator_num,]$covariate_names |> paste(collapse = " - ")
    outcome_regression_model <- df_estimators_parameters[chosen_estimator_num, ]$outcome_regression_model
    true_effect <- readRDS(file.path(path_results_experiment(), "average_outcome_df.RDS"))[outcome_type == "conditional", AB]
    combined_results_df <- combined_results()
    result_plot <- combined_results_df |>
      dplyr::filter(estimator_num == chosen_estimator_num) |>
      dplyr::filter(variance < 20) |> # filtering absurd variance estimates for graphical exploration
      ggplot() +
      geom_rect(aes(xmin = -Inf, xmax = true_effect, ymin = true_effect, ymax = Inf), fill = "#A0D2AD", alpha = 0.02) +
      geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), fill = "grey", alpha = 0.01) +
      geom_point(aes(x = lb, y = ub, color = anchored), alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
      facet_wrap(~adjustment + model) +
      labs(subtitle = paste("outcome regression model: ", outcome_regression_model, "\n",
                            "covariate names: ", covariate_names, sep = ""))

    return(result_plot)
  }



  output$resultsPlot1 <- bindEvent(path_results_experiment(), x = renderPlot(renderPlotForEstimatorNum(1)))
  output$resultsPlot2 <- bindEvent(path_results_experiment(), x = renderPlot(renderPlotForEstimatorNum(2)))
  output$resultsPlot3 <- bindEvent(path_results_experiment(), x = renderPlot(renderPlotForEstimatorNum(3)))
    
}

shinyApp(ui = ui, server = server)
