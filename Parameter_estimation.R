
require(deSolve)



error_function <- function(result, data){
  
  
  error <- sum(abs(result - data)^2)
  
  
  return(error)
  
}




stochastic_estimation <- function(model, initial, parameters, time_interval,step_size,data,n_iterations){
  
  for (i in 1:n_iterations){
    
    output <- ode(initial,time_interval,model,parameters)
    
    result_output <- output[,2:c(1+length(initial))]
    
    
    param_changing <- unlist(parameters)
    
    
    # calculating the error
    
    error <- error_function(result_output, data)
    
    # make a random sign generator that gives either 1 or -1
    
    rand_sign <- 2*sample(c(0,1), replace=TRUE, size=1)-1
    
    # creating a random number from 1 to the number of 
    
    rand_index <- sample(1:length(parameters),1)
    
    
    param_changing[rand_index] <- param_changing[rand_index] + step_size * rand_sign
    
    # error here - we need to have it as a list in param_changing
    
    new_output <- ode(initial, time_interval, model, param_changing)
    
    new_result_output <- new_output[,2:c(1+length(initial))]
    
    new_error <- error_function(new_result_output, data)
    
    if (new_error < error){
      
      parameters[rand_index] <- param_changing[rand_index]
      
    }
    
  }
  
  
  return(c(parameters,error))
  
}



