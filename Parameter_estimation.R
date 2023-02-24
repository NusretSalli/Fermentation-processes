
rm(list=ls())

require(deSolve)


stochastic_estimation <- function(model, initial, parameters, time_interval,n_iterations){
  
  
  for (i in 1:n_iterations){
    
    output <- ode(initial,time_interval,model,parameters)
    
    
  }
  
  
  
  
}



