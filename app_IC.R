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
                 fluidRow(
                   bs4Card(
                     title = "Select Parameters",
                     uiOutput("dynamic_selectors"), 
                     width = 2
                   ),
                   bs4TabCard(width = 10, 
                              title = "Distributions",
                              tabPanel(
                                title = "Treatment effect",
                                plotOutput("treatment_effect")
                              ),
                              tabPanel(
                                title = "Outcome distribution",
                                plotOutput("outcome_distribution_Output")
                              ), 
                              tabPanel(
                                title = "Covariates distribution",
                                plotOutput("covariates_distribution_Output")
                              )
                   )
                 ),
                 bs4Card(
                   DT::dataTableOutput("selected_parameters"),
                   title = "Selected Parameters",
                   width = 12, 
                   collapsed = TRUE
                 ),
                 uiOutput("results_box")

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

  df_stats <- readRDS(file.path(dir_simulations, "processed_results", "df_stats.RDS"))
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
  output$outcome_distribution_Output <- bindEvent(
    path_results_experiment(),
    x = renderPlot({
      readRDS(file.path(path_results_experiment(), "outcomes_distribution.RDS"))
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
        geom_segment(aes(y = parameter, xend = value, yend = parameter, color = parameter), arrow = arrow(length = unit(0.2, "inches"), type = "closed", angle = 15), arrow.fill = "black", x = 0, linetype = "solid") +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
        geom_vline(xintercept = true_treatment_effect, linetype = "dashed", colour = "red") +
        guides(color = "none") +
        labs(title = "Outcome model: covariates coefficient values")
    })
  )


  renderCIPlot <- function(chosen_estimator_num) {
    covariate_names <- df_estimators_parameters[chosen_estimator_num,]$covariate_names |> paste(collapse = " - ")
    outcome_regression_model <- df_estimators_parameters[chosen_estimator_num, ]$outcome_regression_model
    true_effect <- readRDS(file.path(path_results_experiment(), "average_outcome_df.RDS"))[outcome_type == "conditional", AB]
    combined_results_df <- combined_results()
    result_plot <- combined_results_df |>
      dplyr::filter(estimator_num == chosen_estimator_num) |>
      dplyr::mutate(data_type = ifelse(model %in% c("stc", "maic", "unadjusted"), 
                                       "PAIC", "IPD") |> as.factor() |> relevel(ref = "PAIC"), 
                    adjustment = relevel(as.factor(adjustment), ref = "unadjusted")) |> 
      dplyr::filter(variance < 20) |> # filtering absurd variance estimates for graphical exploration
      ggplot() +
      geom_rect(aes(xmin = -Inf, xmax = true_effect, ymin = true_effect, ymax = Inf), fill = "#A0D2AD", alpha = 0.02) +
      geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf), fill = "grey", alpha = 0.01) +
      geom_point(aes(x = lb, y = ub, color = anchored), alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
      # facet_wrap(~adjustment + model) +
      facet_grid(data_type ~ adjustment) +
      labs(subtitle = paste("outcome regression model: ", outcome_regression_model, "\n",
                            "covariate names: ", covariate_names, sep = ""))

    return(result_plot)
  }

  renderPlotIndicators <- function(chosen_estimator_num) {
    df_stats |> 
      dplyr::filter(indicator != "correct_decision") |> 
      dplyr::filter(Population_parameters_num == 1 &
                      Estimator_num == 1) |>
      dplyr::mutate(objective = ifelse(indicator %in% c("bias", "rmse"), 0, 
                                       ifelse(indicator == "vr", 1, 
                                              ifelse(indicator == "cov_95", 0.95, NA))), 
                    Model = ifelse(Model %in% c("GLM", "ML"), "IPD", 
                                   ifelse(Model == "Unadjusted", "", Model)), 
                    Method = paste(Adjustment, Model, sep = ": ")) |> 
      ggplot() +
      facet_wrap(~indicator, scales = "free") +
      geom_point(aes(y = Method, x = values, color = Anchored), size = 3.5, alpha = 0.6) +
      geom_vline(aes(xintercept = objective), linetype = "dashed") +
      theme_bw() +
      theme(axis.text = element_text(size = 12), 
            strip.text = element_text(size = 12)
            )
  }
  


  output$CIPlot1 <- bindEvent(path_results_experiment(), x = renderPlot(renderCIPlot(1)))
  output$CIPlot2 <- bindEvent(path_results_experiment(), x = renderPlot(renderCIPlot(2)))
  output$CIPlot3 <- bindEvent(path_results_experiment(), x = renderPlot(renderCIPlot(3)))
  
  output$PlotIndicators1 <- bindEvent(path_results_experiment(),
                                    x = renderPlot(renderPlotIndicators(1)))
  output$PlotIndicators2 <- bindEvent(path_results_experiment(),
                                    x = renderPlot(renderPlotIndicators(2)))
  output$PlotIndicators3 <- bindEvent(path_results_experiment(),
                                    x = renderPlot(renderPlotIndicators(3)))
  
    
  output$results_box <- renderUI(
    tagList(
    bs4TabCard(width = 12,
               title = paste0("Estimator results: ", df_estimators_parameters[1,]$outcome_regression_model),
               tabPanel(
                 title = "CI", 
                 closable = TRUE,
                 plotOutput("CIPlot1")
               ),
               tabPanel(
                 title = "Indicators",
                 closable = TRUE,
                 plotOutput("PlotIndicators1")
               )
    ),
    bs4TabCard(width = 12,
               title = paste0("Estimator results: ", df_estimators_parameters[2,]$outcome_regression_model),
               tabPanel(
                 title = "CI", 
                 closable = TRUE,
                 plotOutput("CIPlot2")
               ),
               tabPanel(
                 title = "Indicators",
                 closable = TRUE,
                 plotOutput("PlotIndicators2")
               )
    ),
    bs4TabCard(width = 12,
               title = paste0("Estimator results: ", df_estimators_parameters[3,]$outcome_regression_model),
               tabPanel(
                 title = "CI", 
                 closable = TRUE,
                 plotOutput("CIPlot3")
               ),
               tabPanel(
                 title = "Indicators",
                 closable = TRUE,
                 plotOutput("PlotIndicators3")
               )
    )
    )
  )
  
}

shinyApp(ui = ui, server = server)
