


rm(list=ls())

require(deSolve)

require(ODEsensitivity)

require(epiR)

require(ggplot2)

require(ggpubr)

source("Models.R")

source("Phase_analysis_func.R")

source("Sensitivity_analysis_func.R")



bound_var <- c("rate",
               "flow",
               "G_medium",
               "G50",
               "N_rate_inhib_growth",
               "lac_con_growth",
               "lac_prod_growth",
               "N_rate_inhib_mid",
               "lac_con_mid",
               "lac_prod_mid",
               "n_rate_inhib_max",
               "lac_con_max",
               "lac_prod_max")

bound_min_var <- c(0.001,
                   0.20,
                   50,
                   0.01,
                   0.1,
                   0.1,
                   0.1,
                   30,
                   30,
                   30,
                   0.1,
                   0.1,
                   0.1)

bound_max_var <- c(0.7,
                   0.95,
                   700,
                   400,
                   2,
                   2,
                   2,
                   200,
                   200,
                   200,
                   1,
                   1,
                   1)

time_val <- c(0.01, seq(0.1,30, by = 0.1))

N0 <- 5
G0 <- 400
L0 <- 0

x0 <- c(N = N0, G = G0, L = L0)


sensitive_sobol_final <- sobol_sensitivity(final_model_analysis,
                                           bound_var,
                                           x0,
                                           bound_min_var,
                                           bound_max_var,
                                           time_val)


plot(sensitive_sobol_final, pars_plot = bound_var, state_plot = "N", main_title = "N sensitivity - SOBOL", type = "l", lwd = 2)

plot(sensitive_sobol_final, pars_plot = bound_var, state_plot = "G", main_title = "G sensitivity - SOBOL", type = "l", lwd = 2)

plot(sensitive_sobol_final, pars_plot = bound_var, state_plot = "L", main_title = "L sensitivity - SOBOL", type = "l", lwd = 2)


# morris # 


sensitive_morris_final <- morris_sensitivity(final_model_analysis,
                                             bound_var,
                                             x0,
                                             bound_min_var,
                                             bound_max_var,
                                             time_val)


plot(sensitive_morris_final, pars_plot = bound_var, state_plot = "N", kind = "sep", main_title = "N sensitivity - Morris", type = "l")

plot(sensitive_morris_final, pars_plot = bound_var, state_plot = "G", kind = "sep", main_title = "G sensitivity - Morris", type = "l")

plot(sensitive_morris_final, pars_plot = bound_var, state_plot = "L", kind = "sep", main_title = "L sensitivity - Morris", type = "l")

# PRCC #

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
                "lac_prod_mid",
                "n_rate_inhib_max",
                "lac_con_max",
                "lac_prod_max")

n_iterations <- 1000

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

n_rate_inhib_max_list <- runif(n_iterations, min = 0.1, max = 1)

lac_con_max_list <- runif(n_iterations, min = 0.1, max = 1)

lac_prod_max_list <- runif(n_iterations, min = 0.1, max = 1)


param_data_frame <- cbind(rate_list,
                          flow_list,
                          G_medium_list,
                          G50_list,
                          N_rate_inhib_growth_list,
                          lac_con_growth_list,
                          lac_prod_growth_list,
                          N_rate_inhib_mid_list,
                          lac_con_mid_list,
                          lac_prod_mid_list,
                          n_rate_inhib_max_list,
                          lac_con_max_list,
                          lac_prod_max_list)

results_PRCC <- PRCC_calc(final_model,
                          x0,
                          state_name,
                          param_name,
                          time_val,
                          param_data_frame,
                          n_iterations) 


colnames(results_PRCC)[1:length(param_name)] <- param_name 


# pair plot scatterplot # 

pairs_plot_sim(results_PRCC, state_name, param_name)


# PRCC plot #

PRCC_data <- PRCC_data_maker(results_PRCC, state_name, param_name)


# PRCC plot

PRCC_plot <- PRCC_plot(PRCC_data,status)

PRCC_plot



### Sensitivity of protein model

## Initial state
N0 <- 5
G0 <- 400
L0 <- 0
P0 <- 0

x0 <- c(N = N0, G = G0, L = L0, P = P0)

p <- list()

p$rate <- 0.04

p$flow <- 0.75 # this shouldn't change from 0.75

p$G_medium <- 400

p$G50 <- 50

p$N_rate_inhib_growth <- 0.5

p$lac_con_growth <- 0.5

p$lac_prod_growth <- 0.5

p$N_rate_inhib_mid <- 120

p$lac_con_mid <- 20

p$lac_prod_mid <- 120

p$n_rate_inhib_max <- 1

p$lac_con_max <- 0.9

p$lac_prod_max <- 0.75

p$epsilon <- 0.99

time <- seq(0,30,0.1)

sol <- ode(x0,time,final_protein_model,p)

output <- data.frame(sol)

#### plotting section ####

ggplot(data = output, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
  ylim(0, max(output$G)+5)

ggplot(data = output, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Lactate levels", x = "time", y = "Lactate") + 
  ylim(0, max(output$L)+5)

ggplot(data = output, aes(x = time, y = P)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Protein levels", x = "time", y = "Product") + 
  ylim(0, max(output$L)+5)


##########

param_lower <- 0.4

param_upper <- 0.99

param_name <- c("rate",
                "flow",
                "G_medium",
                "G50",
                "N_rate_inhib_growth",
                "lac_con_growth",
                "lac_prod_growth",
                "N_rate_inhib_mid",
                "lac_con_mid",
                "lac_prod_mid",
                "n_rate_inhib_max",
                "lac_con_max",
                "lac_prod_max",
                "epsilon")


param_optim_output <- param_simulator(param_lower, param_upper, p, param_name, 2, 20)


ggplot(data = param_optim_output, aes(x = time, y = P)) + geom_point(aes(color = flow),size = 1) +
  labs(title = "Process optimization", x = "time", y = "P")+
  scale_color_gradientn(colours = rainbow(10))

ggplot(data = param_optim_output, aes(x = time, y = N)) + geom_point(aes(color = flow),size = 1) +
  labs(title = "Process optimization", x = "time", y = "N")+
  scale_color_gradientn(colours = rainbow(10))

ggplot(data = param_optim_output, aes(x = time, y = G)) + geom_point(aes(color = flow),size = 1) +
  labs(title = "Process optimization", x = "time", y = "G")+
  scale_color_gradientn(colours = rainbow(10))

ggplot(data = param_optim_output, aes(x = time, y = L)) + geom_point(aes(color = flow),size = 1) +
  labs(title = "Process optimization", x = "time", y = "N")+
  scale_color_gradientn(colours = rainbow(10))


##### 

# simulating all of them! 

param_lower_list <- c(0.01, # rate
                      0.4, # flow
                      200, # G_medium
                      10, # G50
                      0.1, # N_rate_inhib_growth
                      0.1, # lac_con_growth
                      0.1, # lac_prod_growth
                      40, # N_rate_inhib_mid
                      40, # lac_con_mid
                      40, # lac_prod_mid
                      0.5, # n_rate_inhib_max
                      0.5, # lac_con_max
                      0.5, # lac_prod_max
                      0.4) # epsilon

param_upper_list <- c(0.4, # rate
                      0.90, # flow
                      600, # G_medium
                      100, # G50
                      1, # N_rate_inhib_growth
                      1, # lac_con_growth
                      1, # lac_prod_growth
                      220, # N_rate_inhib_mid
                      220, # lac_con_mid
                      220, # lac_prod_mid
                      1, # n_rate_inhib_max
                      1, # lac_con_max
                      1, # lac_prod_max
                      0.99) # epsilon


result <- param_simulator_plotter(param_lower_list, param_upper_list, p, param_name, 30)

result[["epsilon"]]


# show the plot by typing result[["name of the variable"]]


## Contour plot optimization ##




















