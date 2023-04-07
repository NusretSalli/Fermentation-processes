
require(deSolve)


final_model_estimation <- function(t,x,p){
  
  with(as.list(c(x,p)), {
    
    N_rate_inhib <- (1 / (1 + exp(N_rate_inhib_growth*(L-N_rate_inhib_mid))))
    
    lac_con <- (1 / (1 + exp(lac_con_growth*(G-lac_con_mid))))
    
    lac_prod <- (1 / (1 + exp(-lac_prod_growth*(G-lac_prod_mid))))
    
    
    dN <- rate * N * G / (1 + G / G50) * N_rate_inhib + N * L * lac_con - flow * N
    
    dG <- - rate * N * G / (1 + G / G50) * N_rate_inhib + flow * (G_medium - G) - N * lac_prod
    
    dL <- N * lac_prod - flow * L - N * L * lac_con
    
    
    return(list(c(dN, dG, dL)))
    
    
  })
  
}


add_noise <- function(real_data, mean_val, sd_val){
  
  nrow <- dim(real_data)[1]
  
  ncol <- dim(real_data)[2]
  
  noise <- matrix(rnorm(nrow * (ncol-1), mean = mean_val, sd = sd_val), ncol = ncol-1)
  
  noised_data <- real_data[,2:ncol] + noise
  
  real_data[,2:ncol] <- noised_data
  
  return(real_data)
  
}





stochastic_estimation <- function(model, initial, parameters,param_name, time_interval,step_size,data,n_iterations){
  
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
    
    
    param_changing[rand_index] <- param_changing[rand_index] + step_size[rand_index] * rand_sign
    
    # error here - we need to have it as a list in param_changing
    
    param_changed <- list()
    
    for (l in 1:length(param_name)){
      
      param_changed[l] <- param_changing[l]
      
    }
    
    names(param_changed) <- param_name
    
    
    new_output <- ode(initial, time_interval, model, param_changed)
    
    new_result_output <- new_output[,2:c(1+length(initial))]
    
    new_error <- error_function(new_result_output, data)
    
    if (new_error < error){
      
      parameters[rand_index] <- param_changed[rand_index]
      
    }
    
  }
  
  
  return(c(parameters,error))
  
}


fit_simulation <- function(objective_func, real_data, param_to_fit, bound_min_var, bound_max_var, n_simulations){
  
  
  for (i in 1:n_simulations){
    
    noised_data <- add_noise(real_data, mean = 12, sd = 6)
    
    fit <- modFit(objective_func,
                  p = param_to_fit,
                  lower = bound_min_var,
                  upper = bound_max_var,
                  method = "L-BFGS-B") 
    
  }
  
  
  
}
  
  
  
  





