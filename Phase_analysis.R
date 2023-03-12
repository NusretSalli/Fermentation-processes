

rm(list=ls())

require(deSolve)

require(ggplot2)

require(MASS)

source("Phase_analysis_func.R")

source("Models.R")

require(plotly)

## Initial state
N0 <- 5
G0 <- 395
L0 <- 0

x0 <- c(N = N0, G = G0, L = L0)

p <- list()

p$rate <- 0.05

p$flow <- 0.75 # this shouldn't change from 0.75

p$G_medium <- 400

p$G50 <- 50

p$N_rate_inhib_growth <- 0.2

p$lac_con_growth <- 0.2

p$lac_prod_growth <- 0.2

p$N_rate_inhib_mid <- 108

p$lac_con_mid <- 20

p$lac_prod_mid <- 130

p$n_rate_inhib_max <- 0.8

p$lac_con_max <- 0.9

p$lac_prod_max <- 0.9

time <- seq(0,30,0.1)

sol <- ode(x0,time,final_model,p)

output <- data.frame(sol)


state_output <- output[,c(2,3,4)]

m <- ggplot(state_output, aes(x = N, y = G)) +
  geom_point() +
  xlim(0, 350) +
  ylim(0, 400)

# contour lines
m + geom_density_2d_filled(alpha = 0.5) +
  geom_density_2d(linewidth = 0.25, colour = "black")


N0_list <- runif(1000, 5, 50)

G0_list <- runif(1000, 100, 350)

L0_list <- p$G_medium - N0_list - G0_list


x0_state <- c(N = N0_list[1], G = G0_list[1], L = L0_list[1])
  
sol_iter <- ode(x0_state,time,final_model,p)

output <- data.frame(sol_iter)

ggplot(data = output, aes(x = N, y = G)) + geom_point(size = 3, color = "blue") + 
  labs(title = "Phase N / G", x = "N", y = "G")

ggplot(data = output, aes(x = N, y = L)) + geom_point(size = 3, color = "blue") + 
  labs(title = "Phase N / L", x = "N", y = "L")

ggplot(data = output, aes(x = G, y = L)) + geom_point(size = 3, color = "blue") + 
  labs(title = "Phase G / L", x = "G", y = "L")  

fig <- plot_ly(output, x = ~N, y = ~G, z = ~L, type = 'scatter3d', mode = 'lines',
               opacity = 1, line = list(width = 6))

fig  






