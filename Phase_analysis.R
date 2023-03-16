

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

output_matrix <- as.matrix(output)

output_int <- data.frame(output)

output_int[,c(1,2,3,4)] <- lapply(output[,c(1,2,3,4)], function(x) ceiling(x))

x_grid <- seq(min(output_int$G),max(output_int$N),length.out = length(output$L))

y_grid <- seq(min(output_int$G),max(output_int$G),length.out = length(output$L))



## WHAT UFFE WANTS ##

ggplot(data = output, aes(x = N, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")


ggplot(data = output, aes(x = N, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output, aes(x = G, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")


# when doing with different initial values #

n_iterations <- 10

N0_list <- runif(n_iterations, min = 1, max = 60)

G0_list <- runif(n_iterations, min = 200, max = 340)

L0_list <- p$G_medium - N0_list - G0_list


for (i in 1:n_iterations){
  
  x0_list <- c(N = N0_list[i], G = G0_list[i], L = L0_list[i])
  
  sol_list <- ode(x0_list,time,final_model,p)
  
  if (i == 1){
    
    output_final_list <- data.frame(sol_list)
    
  }
  
  output_list <- data.frame(sol_list)
  
  output_final_list <- rbind(output_final_list,output_list)
  
}

run <- rep(c(rep(1, 301), rep(2, 301), rep(3, 301), rep(4, 301), rep(5,301), rep(6,301),rep(7,301),rep(8,301),rep(9,301),rep(10,301),rep(11,301)))

output_final_final_list <- cbind(output_final_list, run)


output_compare <- phase_plane_plot_all(n_iterations)

ggplot(data = output_compare, aes(x = N, y = G)) + geom_point(size = 3, color =) +
  labs(title = "NG phaseplot", x = "N", y = "G")


sol_list <- ode(x0_list,time,final_model,p)

output_list <- data.frame(sol_list)

ggplot(data = output_list, aes(x = N, y = G)) + geom_point(size = 3, color = "blue") +
  labs(title = "NG phaseplot", x = "N", y = "G")


ggplot(data = output_list, aes(x = N, y = L)) + geom_point(size = 3, color = "blue") + 
  labs(title = "NL phaseplot", x = "N", y = "L")

ggplot(data = output_list, aes(x = G, y = L)) + geom_point(size = 3, color = "blue") + 
  labs(title = "GL phaseplot", x = "G", y = "L")



## Trying to make contour plots need to get L-data for all N / G combination

full_grid <- expand.grid(x_grid,y_grid)


full_grid$L <- output_int$L


# output[,c(2,3,4)] <- lapply(output[,c(2,3,4)], function(x) ceiling(x))

plot1 <- ggplot(full_grid, aes(x = Var1, y = Var2, z = L)) +
  geom_contour_filled()


plot1





# 
# plot_N_G <- ggplot(output, aes(x = N, y = G)) +
#   geom_point() +
#   xlim(-60, 370) +
#   ylim(-50, 470) +
#   stat_density_2d_filled(alpha = 0.5) +
#   stat_density_2d(linewidth = 0.25, colour = "black")
#   
# 
# plot_N_G
# 
# plot_N_L <- ggplot(output, aes(x = N, y = L)) +
#   geom_point() +
#   xlim(-30, 370) +
#   ylim(-30, 120) +
#   stat_density_2d_filled(alpha = 0.5) +
#   stat_density_2d(linewidth = 0.25, colour = "black")
# 
# plot_N_L
# 
# plot_G_L <- ggplot(output, aes(x = L, y = G)) +
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)+
#   stat_density_2d(linewidth = 0.25, colour = "black")
# 
# plot_G_L


# N0_list <- runif(1000, 5, 50)
# 
# G0_list <- runif(1000, 100, 350)
# 
# L0_list <- p$G_medium - N0_list - G0_list
# 
# 
# x0_state <- c(N = N0_list[1], G = G0_list[1], L = L0_list[1])
#   
# sol_iter <- ode(x0_state,time,final_model,p)
# 
# output <- data.frame(sol_iter)
# 
# ggplot(data = output, aes(x = N, y = G)) + geom_point(size = 3, color = "blue") + 
#   labs(title = "Phase N / G", x = "N", y = "G")
# 
# ggplot(data = output, aes(x = N, y = L)) + geom_point(size = 3, color = "blue") + 
#   labs(title = "Phase N / L", x = "N", y = "L")
# 
# ggplot(data = output, aes(x = G, y = L)) + geom_point(size = 3, color = "blue") + 
#   labs(title = "Phase G / L", x = "G", y = "L")  
# 
# fig <- plot_ly(output, x = ~N, y = ~G, z = ~L, type = 'scatter3d', mode = 'lines',
#                opacity = 1, line = list(width = 6))
# 
# fig  


# ##
# 
# fig <- plot_ly(
#   x = seq(0,max(plotly_state_output[,1]),50), 
#   y = seq(0,max(plotly_state_output[,2]),50),
#   z = plotly_state_output, 
#   type = "contour" 
# )
# 
# fig


# plotly_state_output <- as.matrix(output[,c(2,3,4)])
# 
# fig <- plot_ly(z = ~plotly_state_output, type = "contour")
# 
# fig

# fig <- plot_ly(
#   x = seq(0,max(output_matrix[,4]),10), 
#   y = seq(0,max(output_matrix[,3]),10),
#   z = output_matrix[,c(3,4)], 
#   type = "contour" 
# )
# 
# fig

# gg = ggplot(output_int, aes(x = N, y = G)) + 
#   geom_raster(aes(fill = L)) + 
#   geom_contour(aes(z = L), colour = "white", size = 0.2, alpha = 0.5) + 
#   geom_text_contour(aes(z = L),  colour = "white" ) +
#   labs(x = "Filtered light intensity (umol photons/m2/s)", 
#        y = "Dry Biomass Concentration (g/L)", 
#        fill = "Productivity (umol O2/g biomass/s)") + 
#   theme(legend.title = element_text(size = 10, face = "bold"), 
#         legend.position = "top", panel.background = element_blank(), 
#         axis.text = element_text(colour = "black", size = 10, face = "bold"), 
#         axis.title = element_text(size = 12, face = "bold"), 
#         legend.text = element_text(size = 11), legend.key = element_blank()) + 
#   scale_fill_continuous(low = "#BFE1B0", high = "#137177") + 
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_continuous(expand = c(0,0)) 
# 
# gg



