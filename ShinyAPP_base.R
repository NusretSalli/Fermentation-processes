
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
  titlePanel("Base model demonstration"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      sliderInput(inputId = "N0",
                  label = "Intial value for N",
                  min = 0, max = 500, value = 5),
      
      sliderInput(inputId = "G0",
                  label = "Intial value for G",
                  min = 0, max = 500, value = 50),
      
      sliderInput(inputId = "rate", 
                  label = "Rate value",
                  min = 0, max = 0.2, value = 0.1),
      
      sliderInput(inputId = "flow", 
                  label = "Flow value:",
                  min = 0, max = 1, value = 0.75),
      
      sliderInput(inputId = "G_medium", 
                  label = "G_medium value:",
                  min = 0, max = 700, value = 400),
      
      sliderInput(inputId = "time_end", 
                  label = "time end:",
                  min = 0, max = 100, value = 30),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: plots ----
      plotOutput(outputId = "N_plot"),
      
      plotOutput(outputId = "G_plot"),
      
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
         "G_medium" = input$G_medium)
    
  })
  
  initial_states <- reactive({
    
    c(N = input$N0, G = input$G0)
    
  })
  
  time_end <- reactive({
    
    seq(0,input$time_end,1)
    
  })
  
  
  output$N_plot <- renderPlot({
    
    sol <- ode(initial_states(),time_end(),base_model,diff_eq_values())
    
    output <- data.frame(sol)
    
    ggplot(data = output, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
      labs(title = "Number of cells", x = "time", y = "number of cells")
    
    
  })
  
  output$G_plot <- renderPlot({
    
    sol <- ode(initial_states(),time_end(),base_model,diff_eq_values())
    
    output <- data.frame(sol)
    
    ggplot(data = output, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
      labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
      ylim(0, max(output$G)+5)
    
    
  })
  
  
  output$PRCC_plot <- renderPlot({
    
    
    state_name <- c("N", "G")
    
    param_name <- c("rate", "flow", "G_medium")
    
    n_iterations <- 100
    
    rate_list <- runif(n_iterations, min = 0.001, max = 0.2)
    
    flow_list <- runif(n_iterations, min = 0.2, max = 0.95)
    
    G_medium_list <- runif(n_iterations, min = 50, max = 700)
    
    param_data_frame <- cbind(rate_list,
                              flow_list,
                              G_medium_list)
    
    results_PRCC <- PRCC_calc(base_model,
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
