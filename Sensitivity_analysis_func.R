

library(docstring)


base_model_analysis <- function(t,x,p){
  
  #' Base model used solely on analysis (sensitivity and parameter estimation)
  #' 
  #' Inputs:
  #' 
  #' t -> time vector 
  #' x -> essentially our states
  #' p -> parameters used in the ODEs
  #' 
  #' Uses:
  #' 
  #' Used to do sensiivity analysis (sobol, morris, PRCC) or by doing parameter estimation
  #' 
  
  # defining the differential equations and calculating them based on x and p
  
  with(as.list(c(x,p)), {
    
    dN <- rate * N * G - flow * N
    
    dG <- flow * (G_medium - G) - rate * N * G
    
    return(list(c(dN, dG)))
    
  })
}

lactate_model_analysis <- function(t,x,p){
  
  #' lactate_model used solely on analysis (sensitivity and parameter estimation)
  #' 
  #' Inputs:
  #' 
  #' t -> time vector 
  #' x -> essentially our states
  #' p -> parameters used in the ODEs
  #' 
  #' Uses:
  #' 
  #' Used to do sensiivity analysis (sobol, morris, PRCC) or by doing parameter estimation
  #' 
  
  # defining the differential equations and calculating them based on x and p
  
  with(as.list(c(x,p)), {
    
    logL <- (1 / (1 + exp(L_log_growth*(G-L_log_mid))))
    
    logN <- (0.7 / (1 + exp(N_log_growth*(L-N_log_mid)))) + 0.3
    
    dN <- rate * N * G * logN - flow * N
    
    dG <- flow * (G_medium - G) - rate * N * G + L * logL
    
    dL <- G * N * l_rate - flow * L - L * logL
    
    return(list(c(dN, dG, dL)))
    
    
  })
  
}

final_model_analysis <- function(t,x,p){
  
  #' final_model used solely on analysis (sensitivity and parameter estimation)
  #' 
  #' Inputs:
  #' 
  #' t -> time vector 
  #' x -> essentially our states
  #' p -> parameters used in the ODEs
  #' 
  #' Uses:
  #' 
  #' Used to do sensiivity analysis (sobol, morris, PRCC) or by doing parameter estimation
  #' 
  
  # defining the differential equations and calculating them based on x and p
  
  
  with(as.list(c(x,p)), {
    
    # defining the logistic regression functions used as inhibitors
    
    N_rate_inhib <- (N_rate_inhib_max / (1 + exp(N_rate_inhib_growth*(L-N_rate_inhib_mid)))) + (1-N_rate_inhib_max)
    
    lac_con <- (lac_con_max / (1 + exp(lac_con_growth*(G-lac_con_mid)))) + (1-lac_con_max)
    
    lac_prod <- (lac_prod_max / (1 + exp(-lac_prod_growth*(G-lac_prod_mid)))) + (1-lac_prod_max)
    
    
    
    dN <- rate * N * G / (1 + G / G50) * N_rate_inhib + N * L * lac_con - flow * N
    
    dG <- - rate * N * G / (1 + G / G50) * N_rate_inhib + flow * (G_medium - G) - N * lac_prod
    
    dL <- N * lac_prod - flow * L - N * L * lac_con
    
    
    return(list(c(dN, dG, dL), c(N_rate_inhib = N_rate_inhib, lac_con = lac_con, lac_prod = lac_prod)))
    
    
  })
  
}





sobol_sensitivity <- function(model, var_pars, x0_init, var_min, var_max, time_val, n_iterations = 3000){
  
  #' sensitivity analysis maker based on Sobol's method
  #' 
  #' Inputs:
  #' 
  #' model -> the ODEs which you want to analyze
  #' var_pars -> parameters whose sensitivity will be calculated
  #' x0_init -> initital value for each state in the system of ODEs
  #' var_min -> the smallest value the parameters can have (written as a vector)
  #' var_max -> the largest value the parameters can have (written as a vector)
  #' time_val -> the time interval, which the sensitivity analysis will take place.
  #' n_iterations -> the number of iterations that will be run.
  #' 
  #' Output:
  #' 
  #' sobol_result -> 
  #' 
  #' Uses:
  #' 
  #' Used to do sobol sensitivity analysis
  #' 
  
  # calculating the sobol_result data.frame
  
  sobol_result <- ODEsobol(mod = model, # the model that is used
                           pars = var_pars, # the parameters present
                           state_init = x0_init, # initial state value
                           times = time_val, # time value
                           n = n_iterations, # number of iterations
                           rfuncs = "runif", # how to distribute the parameter values
                           rargs = paste0("min = ", var_min, # minimum value the parameters can have
                                          ", max = ", var_max), # maximum value the parameters can have
                           sobol_method = "Martinez", # which sobol method to be used
                           ode_method = "lsoda", # which ode method to be used (this is default in ode)
                           parallel_eval = TRUE, # parallel evaluation? yes or no.
                           parallel_eval_ncores = 2) # number of cores to be used to evaluate
  
  return(sobol_result)
  
}

morris_sensitivity <- function(model, var_pars, x0_init, var_min, var_max, time_val){
  
  #' sensitivity analysis maker based on Morris' method
  #' 
  #' Inputs:
  #' 
  #' model -> the ODEs which you want to analyze
  #' var_pars -> parameters whose sensitivity will be calculated
  #' x0_init -> initital value for each state in the system of ODEs
  #' var_min -> the smallest value the parameters can have (written as a vector)
  #' var_max -> the largest value the parameters can have (written as a vector)
  #' time_val -> the time interval, which the sensitivity analysis will take place.
  #' 
  #' 
  #' Output:
  #' 
  #' morris_result -> 
  #' 
  #' Uses:
  #' 
  #' Used to do morris sensitivity analysis
  #' 
  
  # calculating the morris_result data.frame
  
  
  ODEmorris(model,
            pars = var_pars, # the parameter
            state_init = x0_init, # intiial states
            times = time_val, # time value
            binf = var_min, # minimum value the parameters can have
            bsup = var_max, # maximum value the parameters can have
            r = 3000, # number of "iterations"
            design = list(type = "oat", # oat = one at a time
                          levels = 10,
                          grid.jump =1), # grid jump should be levels / 2
            scale = TRUE,
            ode_method = "lsoda", # which ode-solver that will be used
            parallel_eval = TRUE, # parallel evaluation? Yes or no?
            parallel_eval_ncores = 2) # number of cores to be used to evaluate
}


PRCC_calc <- function(model, state_init, state_name, param_name, time_val, param_data_frame, n_iterations){
  
  #' Doing the PRCC calculation
  #' 
  #' Inputs:
  #' 
  #' model -> the ODEs which you want to analyze
  #' state_init -> initial state value
  #' state_name -> Name of the states
  #' param_name -> vector of the param names.
  #' time_val -> time value vector
  #' param_data_frame -> The data frame you get after PRCC simulation
  #' 
  #' 
  #' Output:
  #' 
  #' morris_result -> 
  #' 
  #' Uses:
  #' 
  #' Used to do morris sensitivity analysis
  #' 
  
  # calculating the morris_result data.frame
  
  
  
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
      
      #state_value[i,k] <- mean(tail(output[,k+1],n=3))
      
      state_value[i,k] <- mean(output[,k+1])
      
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


param_simulator <- function(param_lower, param_upper,param_list,param_name_list,param_index, n_iterations){
  
  parameter_list <- list()
  
  param_ran_list <- seq(from = param_lower, to = param_upper, length.out = n_iterations)
  
  for (i in 1:n_iterations){
    
    #param_ran <- runif(1, min = param_lower, max = param_upper)
    
    param_ran <- param_ran_list[i]
    
    
    for (l in 1:length(param_list)){
      
      if (l == param_index){
        
        parameter_list[l] <- param_ran
        
      } else {
        
        parameter_list[l] <- param_list[[l]]
        
      }
      
    }
    
    
    names(parameter_list) <- param_name_list
    
    sol <- ode(x0,time,final_protein_model,parameter_list)
    
    if (i == 1){
      
      output_final_list <- data.frame(sol)
      
      parameter <- c(rep(param_ran,length(time)))
      
    } else {
      
      parameter_conc <- rep(param_ran,length(time))
      
      parameter <- c(parameter, parameter_conc)
      
      output_list <- data.frame(sol)
      
      output_final_list <- rbind(output_final_list,output_list)
      
      
    }
    
    
  }
  
  output_total <- cbind(output_final_list, parameter)
  
  colnames(output_total)[ncol(output_total)] <- param_name_list[param_index]
  
  return(output_total)
  
}


param_simulator_plotter <- function(param_lower_list,
                                    param_upper_list,
                                    param_list,
                                    param_name_list,
                                    n_iterations){
  
  plot_list <- list()
  
  for (i in 1:length(param_name_list)){
    
    output_iter <- param_simulator(param_lower = param_lower_list[i],
                                   param_upper = param_upper_list[i],
                                   param_list,
                                   param_name_list,
                                   param_index = i,
                                   n_iterations)
    
    #par(mfrow=c(2,2))
    
    plot1 <- ggplot(data = output_iter, aes(x = time, y = P)) + geom_point(aes(color = output_iter[,ncol(output_iter)]),size = 1) +
      labs(title = "P levels", x = "time", y = "P")+
      scale_color_gradientn(colours = rainbow(10), name = param_name_list[i])
    
    plot2 <- ggplot(data = output_iter, aes(x = time, y = N)) + geom_point(aes(color = output_iter[,ncol(output_iter)]),size = 1) +
      labs(title = "N levels", x = "time", y = "N")+
      scale_color_gradientn(colours = rainbow(10), name = param_name_list[i])
    
    plot3 <- ggplot(data = output_iter, aes(x = time, y = G)) + geom_point(aes(color = output_iter[,ncol(output_iter)]),size = 1) +
      labs(title = "G levels", x = "time", y = "G")+
      scale_color_gradientn(colours = rainbow(10), name = param_name_list[i])
    
    plot4 <- ggplot(data = output_iter, aes(x = time, y = L)) + geom_point(aes(color = output_iter[,ncol(output_iter)]),size = 1) +
      labs(title = "L levels", x = "time", y = "L")+
      scale_color_gradientn(colours = rainbow(10), name = param_name_list[i])
    
    plot_total <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2, align = 'h', common.legend = TRUE)
    
    
    plot_list[[i]] <- plot_total
    
  }
  
  names(plot_list) <- param_name_list
  
  
  return(plot_list)
  
  
}




