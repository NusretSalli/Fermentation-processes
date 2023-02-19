
rm(list=ls())

require(deSolve)

require(ggplot2)

require(tidyr)

require(phaseR)

require(ODEsensitivity)

require(epiR)

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


time <- seq(0,30,1)

sol <- ode(x0_base,time,base_model,p)

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

# base model # 

bound_pars <- c("rate", "flow", "G_medium")

bound_min <- c(0.001, 0.20, 50)

bound_max <- c(0.15, 0.95, 700)

time_val <- c(0.01, seq(1, 30, by = 1))

# sobol #

sensitive_sobol <- sobol_sensitivity(base_model_sensitivity, bound_pars, x0_base, bound_min, bound_max, time_val)


plot(sensitive_sobol, pars_plot = bound_pars, state_plot = "N", main_title = "N sensitivity - SOBOL", type = "l", lwd = 2)

plot(sensitive_sobol, pars_plot = bound_pars, state_plot = "G", main_title = "G sensitivity - SOBOL", type = "l", lwd = 2)

# morris # 


sensitive_morris <- morris_sensitivity(base_model_sensitivity,
                                       bound_pars, x0_base,
                                       bound_min,
                                       bound_max,
                                       time_val)

plot(sensitive_morris, pars_plot = bound_pars, state_plot = "N", kind = "sep", main_title = "N sensitivity - Morris", type = "l")

plot(sensitive_morris, pars_plot = bound_pars, state_plot = "G", kind = "sep", main_title = "G sensitivity - Morris", type = "l")

# PRCC #

state_name <- c("N", "G")

param_name <- c("rate", "flow", "G_medium")

n_iterations <- 500

rate_list <- runif(n_iterations, min = 0.001, max = 0.2)

flow_list <- runif(n_iterations, min = 0.2, max = 0.95)

G_medium_list <- runif(n_iterations, min = 50, max = 700)

param_data_frame <- cbind(rate_list, flow_list, G_medium_list)

results_PRCC <- PRCC_calc(base_model, x0_base, state_name, param_name, time, param_data_frame, 500) 



# pair plot scatterplot # 

pairs(sim_result)


N_PRCC <- epi.prcc(sim_result[,c(1,2,3,4)])

G_PRCC <- epi.prcc(sim_result[,c(1,2,3,5)])

state <- c("N", "N", "N", "G", "G", "G")

NG_PRCC <- rbind(N_PRCC[1:4], G_PRCC[1:4])

NG_PRCC$est <- round(NG_PRCC$est, digits = 2)

total_PRCC <- cbind(NG_PRCC, state)



# plotting #

ggplot(data = total_PRCC, aes(x = var,
                              weight = est,
                              ymin = lower,
                              ymax = upper,
                              fill = state)) +
  geom_bar(position=position_dodge(),
           aes(y=est),
           stat="identity") +
  
  geom_errorbar(position=position_dodge(width=0.9),
                colour="dark green",
                size = 1.1)+
  labs(y = "PRCC coefficient", title = "PRCC bar plot with 95% conf interval")





