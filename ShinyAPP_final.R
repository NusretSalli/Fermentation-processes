

rm(list=ls())

require(shiny)

require(deSolve)

require(ggplot2)

require(epiR)

source("Models.R")

source("Sensitivity_analysis.R")



# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Final model demonstration"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      sliderInput(inputId = "N0",
                  label = "Intial value for N",
                  min = 0, max = 500, value = 5),
      
      sliderInput(inputId = "G0",
                  label = "Intial value for G",
                  min = 10, max = 500, value = 400),
      
      sliderInput(inputId = "L0",
                  label = "Intial value for L",
                  min = 0, max = 500, value = 0),
      
      sliderInput(inputId = "rate", 
                  label = "Rate value",
                  min = 0, max = 0.8, value = 0.1),
      
      sliderInput(inputId = "flow", 
                  label = "Flow value:",
                  min = 0, max = 1, value = 0.75),
      
      sliderInput(inputId = "G_medium", 
                  label = "G_medium value:",
                  min = 0, max = 700, value = 400),
      
      sliderInput(inputId = "G50", 
                  label = "G50-value",
                  min = 0, max = 500, value = 50),
      
      sliderInput(inputId = "N_rate_inhib_growth", 
                  label = "N inhibition growth rate",
                  min = 0, max = 4, value = 0.5, step = 0.1),
      
      sliderInput(inputId = "lac_con_growth", 
                  label = "Lactate consumption growth rate",
                  min = 0, max = 4, value = 0.5, step = 0.1),
      
      sliderInput(inputId = "lac_prod_growth", 
                  label = "lactate production growth rate",
                  min = 0, max = 4, value = 0.5, step = 0.1),
      
      sliderInput(inputId = "N_rate_inhib_mid", 
                  label = "N inhibition 50% point",
                  min = 10, max = 120, value = 20, step = 1),
      
      sliderInput(inputId = "lac_con_mid", 
                  label = "Lactate consumption 50% point",
                  min = 10, max = 120, value = 20, step = 1),
      
      sliderInput(inputId = "lac_prod_mid", 
                  label = "Lactate production 50% point",
                  min = 10, max = 120, value = 20, step = 1),
      
      sliderInput(inputId = "time_end", 
                  label = "time end:",
                  min = 0, max = 100, value = 30),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: plots ----
      plotOutput(outputId = "N_plot"),
      
      plotOutput(outputId = "G_plot"),
      
      plotOutput(outputId = "L_plot"),
      
      plotOutput(outputId = "PRCC_plot"),
      
    )
  )
)



# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  
  
  # output$distPlot <- renderPlot({
  #   
  #   x    <- faithful$waiting
  #   bins <- seq(min(x), max(x), length.out = input$bins + 1)
  #   
  #   hist(x, breaks = bins, col = "#007bc2", border = "white",
  #        xlab = "Waiting time to next eruption (in mins)",
  #        main = "Histogram of waiting times")
  #   
  # })
  
  diff_eq_values <- reactive({
    
    list("rate" = input$rate,
         "flow" = input$flow,
         "G_medium" = input$G_medium,
         "G50" = input$G50,
         "N_rate_inhib_growth" = input$N_rate_inhib_growth,
         "lac_con_growth" = input$lac_con_growth,
         "lac_prod_growth" = input$lac_prod_growth,
         "N_rate_inhib_mid" = input$N_rate_inhib_mid,
         "lac_con_mid" = input$lac_con_mid,
         "lac_prod_mid" = input$lac_prod_mid)
    
  })
  
  initial_states <- reactive({
    
    c(N = input$N0, G = input$G0, L = input$L0)
    
  })
  
  time_end <- reactive({
    
    seq(0,input$time_end,0.1)
    
  })
  
  
  output$N_plot <- renderPlot({
    
    sol <- ode(initial_states(),time_end(),final_model,diff_eq_values())
    
    output <- data.frame(sol)
    
    ggplot(data = output, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
      labs(title = "Number of cells", x = "time", y = "number of cells")
    
    
  })
  
  output$G_plot <- renderPlot({
    
    sol <- ode(initial_states(),time_end(),final_model,diff_eq_values())
    
    output <- data.frame(sol)
    
    ggplot(data = output, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
      labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
      ylim(0, max(output$G)+5)
    
    
  })
  
  output$L_plot <- renderPlot({
    
    sol <- ode(initial_states(),time_end(),final_model,diff_eq_values())
    
    output <- data.frame(sol)
    
    ggplot(data = output, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
      labs(title = "Lactate levels", x = "time", y = "Lactate") + 
      ylim(0, max(output$L)+5)
    
  })
  
  
  output$PRCC_plot <- renderPlot({
    
    
    
    state_name <- c("N", "G", "L")
    
    param_name <- c("rate",
                    "flow",
                    "G_medium",
                    "G50",
                    "N_rate_inhib_growth",
                    "lac_con_growth",
                    "lac_prod_growth",
                    "N_rate_inhib_mid",
                    "lac_con_mid",
                    "lac_prod_mid")
    
    n_iterations <- 100
    
    rate_list <- runif(n_iterations, min = 0.001, max = 0.2)
    
    flow_list <- runif(n_iterations, min = 0.2, max = 0.95)
    
    G_medium_list <- runif(n_iterations, min = 50, max = 700)
    
    G50_list <- runif(n_iterations, min = 0.01, max = 400)
    
    N_rate_inhib_growth_list <- runif(n_iterations, min = 0.1, max = 2)
    
    lac_con_growth_list <- runif(n_iterations, min = 0.1, max = 2)
    
    lac_prod_growth_list <- runif(n_iterations, min = 0.1, max = 2)
    
    N_rate_inhib_mid_list <- runif(n_iterations, min = 10, max = 120)
    
    lac_con_mid_list <- runif(n_iterations, min = 10, max = 120)
    
    lac_prod_mid_list <- runif(n_iterations, min = 10, max = 120)
    
    param_data_frame <- cbind(rate_list,
                              flow_list,
                              G_medium_list,
                              G50_list,
                              N_rate_inhib_growth_list,
                              lac_con_growth_list,
                              lac_prod_growth_list,
                              N_rate_inhib_mid_list,
                              lac_con_mid_list,
                              lac_prod_mid_list)
    
    results_PRCC <- PRCC_calc(final_model,
                              initial_states(),
                              state_name,
                              param_name,
                              time_end(),
                              param_data_frame,
                              n_iterations) 
    
    
    colnames(results_PRCC)[1:length(param_name)] <- param_name 
    
    
    PRCC_data <- PRCC_data_maker(results_PRCC, state_name, param_name)
    
    status <- rep(state_name, times = 1, each = length(param_name))
    
    PRCC_plot <- PRCC_plot(PRCC_data,status)
    
    PRCC_plot
    
    
  })
  
  
}

shinyApp(ui = ui, server = server)


