

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

source("Parameter_estimation_func.R")


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
              method = "Marq")

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



### TRYING FOR IT ALL TO WORK ###

parameters <- c(rate = 0.199,
                flow = 0.75,
                G_medium = 800,
                G50 = 30,
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


solve_model <- function(pars) {
  
  # initial values
  
  state <- c(N = N0, G = G0, L = L0)
  
  time <- seq(0,30,0.1)
  
  output <- ode(y = state, time, final_model_estimation, parms = pars)
  
  return(output)
}


objective_func <- function(x, parset = names(x)) {
  
  parameters[parset] <- x
  
  output <- solve_model(parameters)
  
  modCost(output, output_real, x = "time")
  
}


param_to_fit <- c(rate = 0.01,
                  flow = 0.60,
                  G_medium = 300)

bound_min_var <- c(0.001,
                   0.20,
                   50)

bound_max_var <- c(0.7,
                   0.95,
                   1000)


fit <- modFit(objective_func,
              p = param_to_fit,
              lower = bound_min_var,
              upper = bound_max_var,
              method = "Marq",
              control = list(
                maxiter = 400,
                ftol = 1e-06,
                ptol = 1e-06,
                gtol = 1e-06))

fit

