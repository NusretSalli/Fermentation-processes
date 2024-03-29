
rm(list=ls())

require(deSolve)

require(ggplot2)

require(tidyr)

require(phaseR)

require(ODEsensitivity)

require(epiR)

require(FME)

require(dplyr)

require(plotly)

source("Models.R")

source("Phase_analysis_func.R")

source("Sensitivity_analysis_func.R")

#source("helper functions/3dplot.R")

#source("helper functions/2dplot.R")

source("Parameter_estimation_func.R")


## Initial state
N0 <- 5
G0 <- 400
L0 <- 0

x0 <- c(N = N0, G = G0, L = L0)

p <- list()

p$rate <- 0.1

p$flow <- 0.75 # this shouldn't change from 0.75

p$G_medium <- 400

p$G50 <- 50

p$N_rate_inhib_growth <- 0.2

p$lac_con_growth <- 0.2

p$lac_prod_growth <- 0.2

p$N_rate_inhib_mid <- 50

p$lac_con_mid <- 60

p$lac_prod_mid <- 30

p$n_rate_inhib_max <- 0.8

p$lac_con_max <- 0.9

p$lac_prod_max <- 0.75

time <- seq(0,30,0.1)

sol <- ode(x0,time,final_model,p)

output <- data.frame(sol)


p_new <- list()

p_new$rate <- 0.1

p_new$flow <- 0.75 # this shouldn't change from 0.75

p_new$G_medium <- 800

p_new$G50 <- 30

p_new$N_rate_inhib_growth <- 0.2

p_new$lac_con_growth <- 0.2

p_new$lac_prod_growth <- 0.2

p_new$N_rate_inhib_mid <- 50

p_new$lac_con_mid <- 60

p_new$lac_prod_mid <- 30

#### plotting section ####


ggplot(data = output, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
  ylim(0, max(output$G)+5)

ggplot(data = output, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Lactate levels", x = "time", y = "Lactate") + 
  ylim(0, max(output$L)+5)

#### PHASE PLANE ANALYSIS ####

#plotting_phase <- flowField(plane_plot_base, xlim = c(-3,3), ylim = c(-3,3), parameters = 1, points = 21, add= FALSE)

#phaseplan_analysis <- phasePlaneAnalysis(plane_plot_base, xlim = c(-3,3), ylim = c(-3,3), tend = 100, parameters = -2, add= FALSE)

## PLOTLY ##

fig <- plot_ly(output, x = ~N, y = ~G, z = ~L, type = 'scatter3d', mode = 'lines',
               opacity = 1, line = list(width = 6))

fig


# phase plan analysis in 2D #

xlim <- c(0,800)

ylim <- c(0,100)

param_list <- c(p$rate,
                p$flow,
                p$G_medium)

state_names <- c("N", "G")

phaseplan <- phase_plan_analysis(phase_plane_plot_base,
                                 xlim,
                                 ylim,
                                 300,
                                 param_list,
                                 state_names)


## SENSITIVITY ANALYSIS ##

# final model sobol #


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
                          time,
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



# parameter estimation

new_sol <- ode(x0,time,final_model,p_new)

output_new <- data.frame(new_sol)

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

n_iterations = 10000

step_size <- c(0.1, 0.1, 10, 5, 0.1, 0.1, 0.1, 5, 5, 5)

resulting_param <- stochastic_estimation(final_model,x0,p,param_name,time,step_size,output_new,n_iterations)


## PARAMETER ESTIMATION - WITH FDE PACKAGE and its FIT() function ##


parameters <- c(rate = 0.199,
                flow = 0.75,
                G_medium = 800,
                G50 = 50,
                N_rate_inhib_growth = 0.2,
                lac_con_growth = 0.2,
                lac_prod_growth = 0.2,
                N_rate_inhib_mid = 50,
                lac_con_mid = 60,
                lac_prod_mid = 50)

N0 <- 5
G0 <- 400
L0 <- 0

init <- c(N = N0, G = G0, L = L0)

time <- seq(0,30,0.1)

output_est <- ode(init,time,final_model_estimation,parameters)


new_param <- c(rate = 0.6,
               flow = 0.75,
               G_medium = 600,
               G50 = 30,
               N_rate_inhib_growth = 0.2,
               lac_con_growth = 0.2,
               lac_prod_growth = 0.2,
               N_rate_inhib_mid = 50,
               lac_con_mid = 60,
               lac_prod_mid = 50)

sol_real <- ode(init,time,final_model_estimation,new_param)

output_real <- data.frame(sol_real)

modcost_test <- modCost(output_est,sol_real, x = "time")


error_functional_test <- function(param){
  
  est <- ode(init,time,final_model_estimation,param)
  
  return(modCost(est,sol_real, x = "time"))
  
  
}

bound_min_var <- c(0.001,
                   0.20,
                   50,
                   10,
                   0.1,
                   0.1,
                   0.1,
                   10,
                   10,
                   10)

bound_max_var <- c(0.7,
                   0.95,
                   1000,
                   400,
                   2,
                   2,
                   2,
                   120,
                   120,
                   120)


fit <- modFit(error_functional_test,
              p = parameters,
              lower = bound_min_var,
              upper = bound_max_var,
              method = "Pseudo")

sol_estimated <- ode(init, time, final_model_estimation,fit$par)

output_estimated <- data.frame(sol_estimated)

fit

ggplot(data = output_real, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output_real, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
  ylim(0, max(output_real$G)+5)

ggplot(data = output_real, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Lactate levels", x = "time", y = "Lactate") + 
  ylim(0, max(output$L)+5)

ggplot(data = output_estimated, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output_estimated, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
  ylim(0, max(output$G)+5)

ggplot(data = output_estimated, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Lactate levels", x = "time", y = "Lactate") + 
  ylim(0, max(output$L)+5)


### demonstration of lactate switch

# test cases

lac_prod_max <- 0.9

lac_prod_mid <- 120

lac_prod_mid2 <- 115

lac_prod_growth <- 0.5

lac_prod_growth2 <- 4.5

G_vector <- seq(50,200,0.5)

lac_prod_example <- (lac_prod_max / (1 + exp(-lac_prod_growth*(G_vector-lac_prod_mid))))

lac_prod_example2 <- (lac_prod_max / (1 + exp(-lac_prod_growth2*(G_vector-lac_prod_mid2))))

plot(G_vector, lac_prod_example, main = "Total",
     col = "blue",
     cex = 1,
     pch = 18,
     xlab = "Glucose levels",
     ylab = "Inhibition factor")
lines(G_vector, lac_prod_example, col = "red",
      lwd = 2)
points(G_vector,lac_prod_example2, col = "black",
       cex = 1,
       pch = 18)
lines(G_vector, lac_prod_example2, col = "green",
      lwd = 2)
axis(side = 1, at = seq(50,200,10))

