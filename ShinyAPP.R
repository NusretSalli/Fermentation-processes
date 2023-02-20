
rm(list=ls())

require(shiny)

require(deSolve)

require(ggplot2)

source("Models.R")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Short demonstration"),
  
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
      
      sliderInput(inputId = "L0",
                  label = "Intial value for L",
                  min = 0, max = 500, value = 0),
      
      sliderInput(inputId = "rate", 
                  label = "Rate value",
                  min = 0, max = 0.2, value = 0.1),
      
      sliderInput(inputId = "flow", 
                  label = "Flow value:",
                  min = 0, max = 1, value = 0.75),
      
      sliderInput(inputId = "G_medium", 
                  label = "G_medium value:",
                  min = 0, max = 700, value = 400)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: plots ----
      plotOutput(outputId = "N_plot"),
      
      plotOutput(outputId = "G_plot"),
      
      plotOutput(outputId = "L_plot")
      
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
  
  
  output$N_plot <- renderPlot({
    
    time <- seq(0,30,1)
    
    sol <- ode(initial_states(),time,base_model,diff_eq_values())
    
    output <- data.frame(sol)
    
    ggplot(data = output, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
      labs(title = "Number of cells", x = "time", y = "number of cells")
    
    
  })
  
  
}

shinyApp(ui = ui, server = server)
