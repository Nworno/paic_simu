# Load required libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(ggthemr)
library(patchwork)
library(data.table)
library(ggrepel)
library(plotly)
library(shinydashboard)
library(bs4Dash)
library(DT)

ggthemr::ggthemr("pale")

source("processing_results.R")
vline_data <- data.frame(
  value = c(0, 0, 1, 0.95, 1),
  indicator = factor(c("Bias", "RMSE", "VR", "95% coverage", "Correct decision"))
)

df_stats <- df_stats |>
  mutate(Model = ifelse(Model == "Unadjusted", "Unadjusted", toupper(Model)),
         indicator =
           case_match(indicator, "bias" ~ "Bias",
                      "rmse" ~ "RMSE",
                      "vr" ~ "VR",
                      "cov_95" ~ "95% coverage",
                      "correct_decision" ~ "Correct decision") |>
           factor(levels = c("Bias", "RMSE", "VR", "95% coverage", "Correct decision"))
  ) |>
  mutate(combined = factor(paste(Data, Anchored)),
         annotation = ifelse(Model == "MAIC", "MAIC", ifelse(Model == "STC", "STC", NA_character_)))

classification_scenario <- classification_scenario |>
  dplyr::select(f_X1, f_X2, bY_A_X1, bY_A_X2, covariate_names, outcome_regression_model, model_marginal, model_conditional, any_nonnormal, everything(), -key) |>
  dplyr::select(!contains(c("X3", "X4")))

ui <- tagList(
  tags$head(
    tags$style(HTML("
                      .sidebar-open .sidebar {
                        display: none;
                      }
                      .sidebar-open:hover .sidebar {
                        display: block;
                      }
                      "))
  ),
  dashboardPage(
  dashboardHeader(title = "Explore Results"),
  dashboardSidebar(
    collapsed = TRUE, 
    sidebarMenu(
      menuItem("Explore Results", tabName = "explore_results", icon = icon("chart-bar")), 
      menuItem("Explore IC", tabName = "explore_IC", icon = icon("chart-bar"))
      # selectInput("population", "Select Population",
      #             choices = unique(df_stats$Population_parameters_num)),
      # selectInput("estimator", "Select Estimator",
      #             choices = unique(df_stats$Estimator_num))
      
      # selectInput("indicator", "Select Indicator",
      #             choices = unique(df_stats$indicator))
    ),
    dashboardBody(
      DTOutput("scenarios_table"),
      uiOutput("list_plots")
    )
  )
  )
)

server <- function(input, output) {

  # filtered_data <- reactive({
  #   df_filtered <- df_stats %>%
  #     dplyr::filter(Population_parameters_num == input$population,
  #                   Estimator_num == input$estimator
  #                   # indicator == input$indicator
  #                   )
  #   return(df_filtered)
  # })

  # output$resultPlot <- renderPlotly({
  #
  #   plot_ly(data = filtered_data(),
  #           x = ~values,
  #           y = ~Adjustment,
  #           color = ~Data,
  #           symbol = ~Anchored,
  #           type = "scatter",
  #           text = ~paste(values),
  #           mode = "markers") |>
  #     layout(title = "Exploration of Results",
  #            xaxis = list(title = "Values"),
  #            yaxis = list(title = "Adjustment")
  #     )
  # })
  output$scenarios_table <- DT::renderDataTable(
    DT::datatable(
      classification_scenario,
      selection = list(mode = 'multiple', selected = c(1), target = 'row'),
      options = list(pageLength = 5, scrollX = TRUE),
      filter = "top",
      )
  )

  selected_row_data <- reactive({
    selected_row <- input$scenarios_table_rows_selected
    if (length(selected_row) == 0) return(NULL)
    dplyr::semi_join(
      df_stats,
      classification_scenario[selected_row, ],
      by = c("Population_parameters_num" = "population_parameters_num",
             "Estimator_num" = "estimator_num"))
  })

  plot_function <- function(df) {
    ggplot_result <- df |>
      ggplot(aes(x = values, y = Adjustment, shape = combined, fill = Data, color = Data)) +
      geom_jitter(height = 0.1, size = 3, stroke = 1.5) +
      # facet_wrap("indicator", nrow = 4, scales = "free") +
      scale_fill_manual(values = alpha(c("#de6757", "#eb9050", "#3262ab"), alpha = 0.5)) +
      scale_color_manual(values = alpha(c("#de6757", "#eb9050", "#3262ab"))) +
      scale_shape_manual(values = c(2, 1, 17, 16, 24, 21), guide = "none") +
      # geom_text_repel(aes(label = annotation),
      #                 point.padding = 0.2,
      #                 box.padding = 0.3,
      #                 segment.curvature = -1e-20,
      #                 arrow = arrow(length = unit(0.015, "npc")), min.segment.length = 0, size = 3) +
      labs(x = NULL, y = NULL) +
      theme(legend.position = "bottom")
  }
  output$list_plots <- renderUI({
    req(selected_row_data())

    lapply(c("Bias", "RMSE", "VR", "95% coverage", "Correct decision"),
           function(indic) {
             renderPlotly({
               ggplot_result <- selected_row_data() |>
                 dplyr::filter(indicator == indic) |>
                 plot_function() +
                 labs(title = indic) +
                 geom_vline(aes(xintercept = value), linetype = "dashed",
                            data = vline_data[vline_data$indicator == indic, ])
               ggplotly(ggplot_result,
                        text = ~paste(values))
             })
           }
    )
  }
  )
}

shinyApp(ui, server)
