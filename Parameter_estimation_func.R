
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

final_model_estimation_lactate_switch <- function(t,x,p){
  
  with(as.list(c(x,p)), {
    
    N_rate_inhib <- (N_rate_inhib_max/ (1 + exp(N_rate_inhib_growth*(L-N_rate_inhib_mid))))
    
    lac_con <- (lac_con_max / (1 + exp(lac_con_growth*(G-lac_con_mid))))
    
    lac_prod <- (lac_prod_max / (1 + exp(-lac_prod_growth*(G-lac_prod_mid))))
    
    
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


## PARAMTER SIMULATION FUNCTIONS ##


solve_model <- function(pars) {
  
  # function to use to solve the model based on initial states and parameters provided
  
  
  # initial values
  
  state <- c(N = N0, G = G0, L = L0)
  
  # time vector
  
  time <- seq(0,30,0.1)
  
  # store the output and returning it
  
  output <- ode(y = state, time, final_model_estimation, parms = pars)
  
  return(output)
}


solve_model_lactate <- function(pars) {
  
  # function to use to solve the model based on initial states and parameters provided
  
  
  # initial values
  
  state <- c(N = N0, G = G0, L = L0)
  
  # time vector
  
  time <- seq(0,30,0.1)
  
  # store the output and returning it
  
  output <- ode(y = state, time, final_model_estimation_lactate_switch, parms = pars)
  
  return(output)
}



fit_simulation <- function(real_data,
                          param_to_fit,
                          bound_min_var,
                          bound_max_var,
                          method,
                          n_simulations){
  
  
  # function that simulates multiple fitting / parameter estimation with addition to noise to the original
  # data set
  
  # first defining the error vector
  
  ssr_vec <- numeric(n_simulations)
  
  # defining the progress-bar function
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_simulations, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  # preallocating the parameters dataframe
  
  parameters_data_frame <- data.frame(matrix(NA, nrow = n_simulations, ncol = length(param_to_fit)))
  
  # changing the column name to the parameters' name.
  
  colnames(parameters_data_frame) <- names(param_to_fit)
  
  # due to noised data needing to be defined constantly we define it here (change it so that we have it outside)
  
  objective_func_sim <- function(x, noised_data, parset = names(x)) {
    
    # we substitute the parameter values with x ()
    
    parameters[parset] <- x
    
    output <- solve_model(parameters) # solving the model
    
    modCost(output, noised_data, x = "time") # calculating the loss
    
  }
  
  # we loop from 1 to the number of simulaionts we want
  
  for (i in 1:n_simulations){
    
    # we create our data by adding noise to the original data
    
    noised_data <- add_noise(real_data, mean = 0, sd = 1)
    
    # we now fit our model to estimate the parameters
    
    fit <- modFit(objective_func_sim,
                  p = param_to_fit,
                  lower = bound_min_var,
                  upper = bound_max_var,
                  method = method,
                  noised_data = noised_data)
    
    # change this later - this is not needed
    
    if (i == 1){
      
      parameters_data_frame[i,] <- fit$par
      
      ssr_vec[i] <- fit$ssr
      
    } else {
      
      parameters_data_frame[i,] <- fit$par
      
      ssr_vec[i] <- fit$ssr
      
    }
    
    # progress bar defined
    
    setTxtProgressBar(pb, i)
    
  }
  
  # returning the parameters and error vector in a list.
  
  return(list(parameters_data_frame, ssr_vec))
}


fit_simulation_lactate <- function(real_data,
                          param_to_fit,
                          bound_min_var,
                          bound_max_var,
                          method,
                          n_simulations){
  
  
  # function that simulates multiple fitting / parameter estimation with addition to noise to the original
  # data set
  
  # first defining the error vector
  
  ssr_vec <- numeric(n_simulations)
  
  # defining the progress-bar function
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_simulations, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  # preallocating the parameters dataframe
  
  parameters_data_frame <- data.frame(matrix(NA, nrow = n_simulations, ncol = length(param_to_fit)))
  
  # changing the column name to the parameters' name.
  
  colnames(parameters_data_frame) <- names(param_to_fit)
  
  # due to noised data needing to be defined constantly we define it here (change it so that we have it outside)
  
  objective_func_sim <- function(x, noised_data, parset = names(x)) {
    
    # we substitute the parameter values with x ()
    
    parameters[parset] <- x
    
    output <- solve_model_lactate(parameters) # solving the model
    
    modCost(output, noised_data, x = "time") # calculating the loss
    
  }
  
  # we loop from 1 to the number of simulaionts we want
  
  for (i in 1:n_simulations){
    
    # we create our data by adding noise to the original data
    
    noised_data <- add_noise(real_data, mean = 0, sd = 1)
    
    # we now fit our model to estimate the parameters
    
    fit <- modFit(objective_func_sim,
                  p = param_to_fit,
                  lower = bound_min_var,
                  upper = bound_max_var,
                  method = method,
                  noised_data = noised_data)
    
    # change this later - this is not needed
    
    if (i == 1){
      
      parameters_data_frame[i,] <- fit$par
      
      ssr_vec[i] <- fit$ssr
      
    } else {
      
      parameters_data_frame[i,] <- fit$par
      
      ssr_vec[i] <- fit$ssr
      
    }
    
    # progress bar defined
    
    setTxtProgressBar(pb, i)
    
  }
  
  # returning the parameters and error vector in a list.
  
  return(list(parameters_data_frame, ssr_vec))
}


histogram_sim_maker <- function(parameter_list, real_parameter_val){
  
  # function that constructs a histogram of the simulated parameters and displays other useful values
  
  # simulating for each parameter that was simulated
  
  for (i in 1:dim(parameter_list)[2]){
    
    # creating the plot instance
    
    plot <- ggplot(parameter_list, aes(x = parameter_list[,i])) + 
      
      # geom histogram where we plot the density / frequency
      geom_histogram(aes(y = after_stat(count / sum(count))), bins = 40, color = "black", fill = "red", alpha = 0.5) +
      
      # we scale the y-axis such that we have it in percent
      scale_y_continuous(labels = scales::percent)+
      
      # we create 3 vertical lines that display the real parameter value, the mean and the median
      geom_vline(aes(xintercept=real_parameter_val[i], color="Real"), linetype = "dashed", linewidth = 1.5) +
      geom_vline(aes(xintercept = mean(parameter_list[,i]),color = "Mean"), linetype = "dashed", linewidth = 1.5) +
      geom_vline(aes(xintercept = median(parameter_list[,i]),color = "Median"), linetype = "dashed", linewidth = 1.5) +
      
      # we create a legend that 
      scale_color_manual(name = "Info", values = c(Real = "blue", Mean = "Red", Median = "black")) +
      
      # displaying the xlabel and ylabel
      xlab(colnames(parameter_list)[i]) +
      ylab("Percent (%)")+
      ggtitle(paste0("results for parameter ",colnames(parameter_list)[i]))
    
    # printing the plots one by one
    print(plot)
    
  }
  
}

# specifically for the initial_state_checker - might change this later


param_initial_state_tester <- function(real_data,
                                       param_names,
                                       bound_min_var,
                                       bound_max_var,
                                       method,
                                       n_simulations){
  
  
  # function that simulates multiple fitting / parameter estimation with addition to noise to the original
  # data set
  
  
  # defining the progress-bar function
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_simulations, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  # preallocating the parameters dataframe
  
  parameters_end_data_frame <- data.frame(matrix(NA, nrow = n_simulations, ncol = length(param_names)+1))
  
  parameters_end_data_frame[,ncol(parameters_end_data_frame)] <- rep("end",n_simulations)
  
  # changing the column name to the parameters' name.
  
  colnames(parameters_end_data_frame) <- c(param_names,"state")
  
  # due to noised data needing to be defined constantly we define it here (change it so that we have it outside)
  
  objective_func_sim_test <- function(x, real_data, parset = names(x)) {
    
    # we substitute the parameter values with x ()
    
    parameters[parset] <- x
    
    output <- solve_model(parameters) # solving the model
    
    modCost(output, real_data, x = "time") # calculating the loss
    
  }
  
  parameters_initial_data_frame <- data.frame(matrix(NA, nrow = n_simulations, ncol = length(param_names)+1))
  
  parameters_initial_data_frame[,ncol(parameters_initial_data_frame)] <- rep("initial",n_simulations)
  
  colnames(parameters_initial_data_frame) <- c(param_names,"state")
  
  for (j in 1:length(param_names)){
    
    ran_param_initial <- runif(n_simulations,bound_min_var[j],bound_max_var[j])
    
    parameters_initial_data_frame[,j] <- ran_param_initial
    
  }
  
  
  # we loop from 1 to the number of simulaionts we want
  
  for (i in 1:n_simulations){
    
    param_to_fit <- parameters_initial_data_frame[i,1:length(param_names)]
    
    names(param_to_fit) <- param_names
    
    
    # we now fit our model to estimate the parameters
    
    fit <- modFit(objective_func_sim_test,
                  p = param_to_fit,
                  lower = bound_min_var,
                  upper = bound_max_var,
                  method = method,
                  real_data = real_data)
    
    # get the parameters
    
    parameters_end_data_frame[i,1:length(param_names)] <- fit$par
    
    
    
    # progress bar defined
    
    setTxtProgressBar(pb, i)
    
  }
  
  # returning the parameters and error vector in a list.
  
  
  total_data_frame <- rbind(parameters_initial_data_frame,parameters_end_data_frame)
  
  return(total_data_frame)
}


# construct error function

error_function <- function(parameters,param_to_plot, parset = names(param_to_plot), real_data){
  
  parameters[parset] <- param_to_plot
  
  output <- solve_model(parameters)
  
  total_err <- sum((real_data - output)^2)
  
  return(total_err)
  
}


param_plot_contour <- function(parameters, param_names, bound_min_var, bound_max_var, dimension, true_data){
  
  
  # construct the sequence for the 2 parameters
  
  param_1 <- seq(bound_min_var[1], bound_max_var[1], length = dimension[1])
  
  names(param_1) <- rep(param_names[1],length(param_1))
  
  param_2 <- seq(bound_min_var[2], bound_max_var[2], length = dimension[2])
  
  names(param_2) <- rep(param_names[2],length(param_2))
  
  results_data_frame <- matrix(NA, nrow = length(param_2), ncol = length(param_1))
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(param_1), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  for (i in 1:length(param_1)){
    
    for (j in 1:length(param_2)){
      
      param_to_plot <- c(param_1[i], param_2[j])
      
      results_data_frame[j,i] <- round(error_function(parameters, param_to_plot, real_data = true_data),3)
      
      
    }
    
    setTxtProgressBar(pb, i)
    
  }
  
  return(list(param_1,param_2,results_data_frame))
}



# error_function_grid <- function(param1,param2){
#   
#   parset = names(c(param1[1],param2[1]))
#   
#   new_param[parset] <- c(param1,param2)
#   
#   output <- solve_model(new_param)
# 
#   total_err <- sum((real_data - output)^2)
# 
#   return(total_err)
# 
# }




