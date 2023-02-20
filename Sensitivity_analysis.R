

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
                           n = 2000,
                           rfuncs = "runif",
                           rargs = paste0("min = ", var_min,
                                          ", max = ", var_max),
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
            r = 2000,
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
  
  sim_result <- data.frame(param_data_frame,
                           state_value)
  
  colnames(sim_result)[1:length(param_name)] <- param_name
  
  return(sim_result)
  
}


pairs_plot_sim <- function(sim_dataframe, state_name, param_name){
  
  for (i in 1:length(param_name)){
    
    reversed_results <- rev(sim_dataframe)
    
    selection_vector <- c(1:length(state_name), length(state_name)+i)
    
    pairs_plotting_data <- reversed_results[,selection_vector]
    
    plot <- pairs(pairs_plotting_data)
    
  }
}


PRCC_data_maker <- function(sim_result, state_name, param_name){
  
  total_PRCC <- c()
  
  for (i in 1:length(state_name)){
    
    state_PRCC <- epi.prcc(sim_result[,c(1:length(param_name),length(param_name)+i)])
    
    total_PRCC <- rbind(total_PRCC, state_PRCC[1:4])
    
  }
  
  status <- rep(state_name, times = 1, each = length(param_name))
  
  data <- cbind(total_PRCC, status)
  
  colnames(data)[1] <- "variables" 
  
  return(data)
  
}



PRCC_plot <- function(PRCC_data, status){
  
  plot <- ggplot(data = PRCC_data, aes(x = variables,
                                weight = est,
                                ymin = lower,
                                ymax = upper,
                                fill = status)) +
    geom_bar(width = 0.6, position=position_dodge(),
             aes(y=est),
             stat="identity") +
    
    geom_errorbar(position=position_dodge(width=0.6),
                  colour="black",
                  size = 1.1, width = 0.4) +
    labs(y = "PRCC coefficient", title = "PRCC bar plot with 95% conf interval") +
    theme(axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 15, face = "bold")) +
    scale_y_continuous(breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
  
  return(plot)
  
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


