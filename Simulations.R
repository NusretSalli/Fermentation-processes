
rm(list=ls())

require(deSolve)

require(ggplot2)

require(tidyr)

require(phaseR)

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


ggplot(data = output, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", size = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", size = 1.5) + 
  labs(title = "Glucose levels", x = "time", y = "glucose levels") + 
  ylim(0, max(output$G)+5)


plot(sol[,2], sol[,3])
lines(sol[,2], sol[,3])


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









