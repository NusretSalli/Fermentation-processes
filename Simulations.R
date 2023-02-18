
rm(list=ls())

require(deSolve)

require(ggplot2)

require(tidyr)

require(phaseR)

require(ODEsensitivity)

require(epiR)

base_model <- function(t,x,p)
{
  ## Unpack state by hand 
  n <- length(x) / 2
  N <- x[1:n]
  G <- x[n + (1:n)]
  
  dN <- numeric(n)
  dG <- numeric(n)
  
  with(p,
       {
         dN <- rate * N * G - flow * N
         dG <- flow * (G_medium - G) - rate * N * G
         
         return(list(c(dN, dG)))
       }
  )
  
}


## Initial state
N0 <- 5
G0 <- 50

x0 <- c(N = N0, G = G0)

p <- list()

p$rate <- 0.06

p$flow <- 0.75 # this shouldn't change from 0.75

p$G_medium <- 400


time <- seq(0,30,1)

sol <- ode(x0,time,base_model,p)

output <- data.frame(sol)



#### plotting section ####


ggplot(data = output, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
  ylim(0, max(output$G)+5)


#### PHASE PLANE ANALYSIS ####


plane_plot <- function(t, y, parameters) {
  I_0 <- parameters
  dy <- numeric(2)
  dy[1] <- 2 *  (y[2] + y[1] - 1/3 * y[1]^3) + I_0
  dy[2] <- (1 - y[1] - y[2])/2   
  return(list(dy))
}

phase_plane_plot <- function(t, func_var, parameters){
  
  rate <- parameters[1]
  
  flow <- parameters[2]
  
  G_medium <- parameters[3]
  
  
  diff_eq <- numeric(2)
  
  
  diff_eq[1] <- rate * func_var[1] * func_var[2] - flow * func_var[1]
  
  diff_eq[2] <- flow * (G_medium - func_var[2]) - rate * func_var[1] * func_var[2]
  
  return(list(diff_eq))
  
}

#plotting_phase <- flowField(plane_plot, xlim = c(-3,3), ylim = c(-3,3), parameters = 1, points = 21, add= FALSE)

#phaseplan_analysis <- phasePlaneAnalysis(plane_plot, xlim = c(-3,3), ylim = c(-3,3), tend = 100, parameters = -2, add= FALSE)

phasplane_analysis <- phasePlaneAnalysis(phase_plane_plot,
                                         xlim = c(0,400),
                                         ylim = c(0,100),
                                         tend = 300,
                                         parameters = c(p$rate, p$flow, p$G_medium),
                                         add= FALSE,
                                         state.names = c("N", "G"))


## SENSITIVITY ANALYSIS ##

# pre-analysis - setting everything up #

base_model_sensitivity <- function(t,x,p){
  with(as.list(c(x,p)), {
    
    dN <- rate * N * G - flow * N
    
    dG <- flow * (G_medium - G) - rate * N * G
    
    return(list(c(dN, dG)))
    
  })
}


bound_pars <- c("rate", "flow", "G_medium")

bound_min <- c(0.001, 0.20, 50)

bound_max <- c(0.15, 0.95, 700)

time_val <- c(0.01, seq(1, 30, by = 1))

# sobol #

sensitive_sobol <- ODEsobol(mod = base_model_sensitivity,
                            pars = bound_pars,
                            state_init = x0,
                            times = time_val,
                            n = 2000,
                            rfuncs = "runif",
                            rargs = paste0("min = ", bound_min,
                                           ", max = ", bound_max),
                            sobol_method = "Martinez",
                            ode_method = "lsoda",
                            parallel_eval = TRUE,
                            parallel_eval_ncores = 2)


plot(sensitive_sobol, pars_plot = bound_pars, state_plot = "N", main_title = "N sensitivity - SOBOL", type = "l", lwd = 2)

plot(sensitive_sobol, pars_plot = bound_pars, state_plot = "G", main_title = "G sensitivity - SOBOL", type = "l", lwd = 2)

# morris # 

sensitive_morris <- ODEmorris(base_model_sensitivity,
          pars = bound_pars,
          state_init = x0,
          times = time_val,
          binf = bound_min,
          bsup = bound_max,
          r = 500,
          design = list(type = "oat",
                        levels = 10,
                        grid.jump =1),
          scale = TRUE,
          ode_method = "lsoda",
          parallel_eval = TRUE,
          parallel_eval_ncores = 2)

sensitive_morris$N

plot(sensitive_morris, pars_plot = bound_pars, state_plot = "N", kind = "sep", main_title = "N sensitivity - Morris", type = "l")

plot(sensitive_morris, pars_plot = bound_pars, state_plot = "G", kind = "sep", main_title = "G sensitivity - Morris", type = "l")

# PRCC #

n_iterations <- 1000

rate_list <- runif(n_iterations, min = 0.001, max = 0.2)

flow_list <- runif(n_iterations, min = 0.2, max = 0.95)

G_medium_list <- runif(n_iterations, min = 50, max = 700)

N_end <- numeric(n_iterations)

G_end <- numeric(n_iterations)

for ( i in 1:n_iterations){
  
  # simulate the differential equations and solve them
  
  parameter_list <- list()
  
  parameter_list$rate <- rate_list[i]
  
  parameter_list$flow <- flow_list[i]
  
  parameter_list$G_medium <- G_medium_list[i]
  
  output <- ode(x0,time,base_model,parameter_list)
  
  N_end[i] <- mean(tail(output[,"N"],n=3))
  
  G_end[i] <- mean(tail(output[,"G"],n=3))
  
  
  
}


sim_result <- data.frame(rate = rate_list,
                         flow = flow_list,
                         G_medium = G_medium_list,
                         N_end,
                         G_end)

# pair plot scatterplot # 

pairs(sim_result)


N_PRCC <- epi.prcc(sim_result[,c(1,2,3,4)]) #

G_PRCC <- epi.prcc(sim_result[,c(1,2,3,5)]) #


total_PRCC <- 



# plotting #

ggplot(N_PRCC) +
  geom_bar(aes(x = N_PRCC$var, y = N_PRCC$est),
           stat = "identity",
           fill="skyblue",
           alpha=0.7)



