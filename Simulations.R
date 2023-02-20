
rm(list=ls())

require(deSolve)

require(ggplot2)

require(tidyr)

require(phaseR)

require(ODEsensitivity)

require(epiR)

require(dplyr)

require(plotly)

source("Models.R")

source("Phase_analysis.R")

source("Sensitivity_analysis.R")



## Initial state
N0 <- 5
G0 <- 50
L0 <- 0

x0 <- c(N = N0, G = G0, L = L0)

x0_base <- c(N = N0, G = G0)

p <- list()

p$rate <- 0.06

p$flow <- 0.75 # this shouldn't change from 0.75

p$G_medium <- 400

p$l_rate <- 0.01

p$L_log_growth <- 0.5

p$N_log_growth <- 0.5

p$L_log_mid <- 45

p$N_log_mid <- 20


time <- seq(0,30,0.1)

sol <- ode(x0,time,lactate_model,p)

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

ggplot(data = output, aes(x = time, y = logisticN.L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "LOGN levels", x = "time", y = "factor") + 
  ylim(0, max(output$logisticN.L))

ggplot(data = output, aes(x = time, y = logisticL.G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "LOGL levels", x = "time", y = "factor") + 
  ylim(0, max(output$logisticL.G))

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

param_list <- c(p$rate, p$flow, p$G_medium)

state_names <- c("N", "G")

phaseplan <- phase_plan_analysis(phase_plane_plot_base,
                                 xlim,
                                 ylim,
                                 300,
                                 param_list,
                                 state_names)


## SENSITIVITY ANALYSIS ##

# base model - sobol # 

bound_pars <- c("rate", "flow", "G_medium")

bound_min <- c(0.001, 0.20, 50)

bound_max <- c(0.15, 0.95, 700)

time_val <- c(0.01, seq(1, 30, by = 1))

sensitive_sobol_base <- sobol_sensitivity(base_model_sensitivity, bound_pars, x0_base, bound_min, bound_max, time_val)

plot(sensitive_sobol_base, pars_plot = bound_pars, state_plot = "N", main_title = "N sensitivity - SOBOL", type = "l", lwd = 2)

plot(sensitive_sobol_base, pars_plot = bound_pars, state_plot = "G", main_title = "G sensitivity - SOBOL", type = "l", lwd = 2)


# lactate model - sobol # 

# change the lactate model such that the logN and LogL's min / max value are changing from 0 to 1 - maybe that will fix it

bound_pars_lac <- c("rate", "flow", "G_medium", "l_rate", "L_log_growth", "N_log_growth", "L_log_mid", "N_log_mid")

bound_min_lac <- c(0.001, 0.20, 50, 0.001, 0.2, 0.2, 20, 10)

bound_max_lac <- c(0.15, 0.95, 700, 0.15, 1.5, 2.5, 70, 60)

time_val <- c(0.01, seq(1,30, by = 1))


sensitive_sobol_lac <- sobol_sensitivity(lactate_model_sensitivity,
                                         bound_pars_lac,
                                         x0,
                                         bound_min_lac,
                                         bound_max_lac,
                                         time_val)


plot(sensitive_sobol_lac, pars_plot = bound_pars_lac, state_plot = "N", main_title = "N sensitivity - SOBOL", type = "l", lwd = 2)

plot(sensitive_sobol_lac, pars_plot = bound_pars_lac, state_plot = "G", main_title = "G sensitivity - SOBOL", type = "l", lwd = 2)

plot(sensitive_sobol_lac, pars_plot = bound_pars_lac, state_plot = "L", main_title = "L sensitivity - SOBOL", type = "l", lwd = 2)


# morris # 


sensitive_morris_lac <- morris_sensitivity(lactate_model_sensitivity,
                                       bound_pars_lac, x0,
                                       bound_min_lac,
                                       bound_max_lac,
                                       time_val)

plot(sensitive_morris_lac, pars_plot = bound_pars_lac, state_plot = "N", kind = "sep", main_title = "N sensitivity - Morris", type = "l")

plot(sensitive_morris_lac, pars_plot = bound_pars_lac, state_plot = "G", kind = "sep", main_title = "G sensitivity - Morris", type = "l")

plot(sensitive_morris_lac, pars_plot = bound_pars_lac, state_plot = "L", kind = "sep", main_title = "L sensitivity - Morris", type = "l")

# PRCC #

state_name <- c("N", "G", "L")

param_name <- c("rate", "flow", "G_medium", "l_rate", "L_log_growth", "N_log_growth", "L_log_mid", "N_log_mid")

n_iterations <- 1000

rate_list <- runif(n_iterations, min = 0.001, max = 0.2)

flow_list <- runif(n_iterations, min = 0.2, max = 0.95)

G_medium_list <- runif(n_iterations, min = 50, max = 700)

l_rate_list <- runif(n_iterations, min = 0.001, max = 0.15)

L_log_growth_list <- runif(n_iterations, min = 0.2, max = 1.5)

N_log_growth_list <- runif(n_iterations, min = 0.2, max = 2.5)

L_log_mid_list <- runif(n_iterations, min = 20, max = 70)

N_log_mid_list <- runif(n_iterations, min = 10, max = 60)

param_data_frame <- cbind(rate_list,
                          flow_list,
                          G_medium_list,
                          l_rate_list,
                          L_log_growth_list,
                          N_log_growth_list,
                          L_log_mid_list,
                          N_log_mid_list)

results_PRCC <- PRCC_calc(lactate_model,
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

# N_PRCC <- epi.prcc(results_PRCC[,c(1:length(param_name),length(param_name)+1)])
# 
# G_PRCC <- epi.prcc(results_PRCC[,c(1:length(param_name),length(param_name)+2)])
# 
# L_PRCC <- epi.prcc(results_PRCC[,c(1:length(param_name),length(param_name)+3)])
# 
# states <- c("N","G","L")
# 
# status <- rep(states, times = 1, each = length(param_name))
# 
# NGL_PRCC <- rbind(N_PRCC[1:4],
#                   G_PRCC[1:4],
#                   L_PRCC[1:4])
# 
# data <- cbind(NGL_PRCC, status)
# 
# colnames(data)[1] <- "variables" 

PRCC_data <- PRCC_data_maker(results_PRCC, state_name, param_name)


# PRCC plot

PRCC_plot <- PRCC_plot(PRCC_data,status)

PRCC_plot







