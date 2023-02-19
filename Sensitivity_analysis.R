

base_model_sensitivity <- function(t,x,p){
  with(as.list(c(x,p)), {
    
    dN <- rate * N * G - flow * N
    
    dG <- flow * (G_medium - G) - rate * N * G
    
    return(list(c(dN, dG)))
    
  })
}

lactate_model_sensitivity <- function(t,x,p){
  
  with(as.list(c(x,p)), {
    
    logL <- (1 / (1 + exp(L_log_growth*(G-L_log_mid))))
    
    logN <- (0.7 / (1 + exp(N_log_growth*(L-N_log_mid)))) + 0.3
    
    dN <- rate * N * G * logN - flow * N
    
    dG <- flow * (G_medium - G) - rate * N * G + L * logL
    
    dL <- G * N * l_rate - flow * L - L * logL
    
    return(list(c(dN, dG, dL)))
    
    
  })
  
}


sobol_sensitivity <- function(model, var_pars, x0_init, var_min, var_max, time_val){
  
  sobol_result <- ODEsobol(mod = model,
                           pars = var_pars,
                           state_init = x0_init,
                           times = time_val,
                           n = 1000,
                           rfuncs = "runif",
                           rargs = paste0("min = ", bound_min,
                                          ", max = ", bound_max),
                           sobol_method = "Martinez",
                           ode_method = "lsoda",
                           parallel_eval = TRUE,
                           parallel_eval_ncores = 2)
  
  return(sobol_result)
  
}

morris_sensitivity <- function(model, var_pars, x0_init, var_min, var_max, time_val){
  
  ODEmorris(model,
            pars = var_pars,
            state_init = x0_init,
            times = time_val,
            binf = var_min,
            bsup = var_max,
            r = 1000,
            design = list(type = "oat",
                          levels = 10,
                          grid.jump =1),
            scale = TRUE,
            ode_method = "lsoda",
            parallel_eval = TRUE,
            parallel_eval_ncores = 2)
}


PRCC_calc <- function(model, state_init, state_name, param_name, time_val, param_data_frame, n_iterations){
  
  state_value <- matrix(data = 0, nrow = n_iterations, ncol = length(state_name))
  
  colnames(state_value) <- state_name
  
  for ( i in 1:n_iterations){
    
    
    parameter_list <- list()
    
    # simulate the differential equations and solve them
    
    for (l in 1:length(param_name)){
      
      parameter_list[l] <- param_data_frame[i,l]
      
      
    }
    
    names(parameter_list) <- param_name
    
    output <- ode(state_init,time_val,model,parameter_list)
    
    for (k in 1:length(state_name)){
      
      state_value[i,k] <- mean(tail(output[,k+1],n=3))
      
      
    }
    
  }
  
  return(state_value)
  
}





# N_end <- numeric(n_iterations)
# 
# G_end <- numeric(n_iterations)
# 
# for ( i in 1:n_iterations){
#   
#   # simulate the differential equations and solve them
#   
#   parameter_list <- list()
#   
#   parameter_list$rate <- rate[i]
#   
#   parameter_list$flow <- flow[i]
#   
#   parameter_list$G_medium <- G_medium[i]
#   
#   output <- ode(x0_base,time,base_model,parameter_list)
#   
#   N_end[i] <- mean(tail(output[,"N"],n=3))
#   
#   G_end[i] <- mean(tail(output[,"G"],n=3))
#   
#   
#   
# }
# 
# 
# sim_result <- data.frame(rate = rate,
#                          flow = flow,
#                          G_medium = G_medium,
#                          N_end,
#                          G_end)


