


plane_plot_base <- function(t, y, parameters) {
  I_0 <- parameters
  dy <- numeric(2)
  dy[1] <- 2 *  (y[2] + y[1] - 1/3 * y[1]^3) + I_0
  dy[2] <- (1 - y[1] - y[2])/2   
  return(list(dy))
}


phase_plane_plot_base <- function(t, func_var, parameters){
  
  rate <- parameters[1]
  
  flow <- parameters[2]
  
  G_medium <- parameters[3]
  
  
  diff_eq <- numeric(2)
  
  
  diff_eq[1] <- rate * func_var[1] * func_var[2] - flow * func_var[1]
  
  diff_eq[2] <- flow * (G_medium - func_var[2]) - rate * func_var[1] * func_var[2]
  
  return(list(diff_eq))
  
}


phase_plane_plot_lactate <- function(t, func_var, parameters){
  
  rate <- parameters[1]
  
  flow <- parameters[2]
  
  G_medium <- parameters[3]
  
  l_rate <- parameters[4]
  
  L_log_growth <- parameters[5]
  
  N_log_growth <- parameters[5]
  
  L_log_mid <- parameters[6]
  
  N_log_mid <- parameters[7]
  
  logL <- (1 / (1 + exp(L_log_growth*(G-L_log_mid))))
  
  logN <- (0.7 / (1 + exp(N_log_growth*(L-N_log_mid)))) + 0.3
  
  
  diff_eq <- numeric(3)
  
  
  diff_eq[1] <- rate * func_var[1] * func_var[2] * logN - flow * func_var[1]
  
  diff_eq[2] <- flow * (G_medium - func_var[2]) - rate * func_var[1] * func_var[2] + func_var[3] * log_L
  
  diff_eq[3] <- func_var[1] * func_var[2] * l_rate - flow * func_var[3] - func_var[3] * log_L
  
  return(list(diff_eq))
  
}



phase_plan_analysis <- function(model, xrange, yrange, t_stop = 300, param_list, state_names){
  
  
  phaseplan_result <- phasePlaneAnalysis(model,
                                         xlim = xrange,
                                         ylim = yrange,
                                         tend = 300,
                                         parameters = param_list,
                                         add= FALSE,
                                         state.names = state_names)
  
  
  return(phaseplan_result)
  
}


phase_plane_data <- function(n_iterations){
  
  
  N0_list <- runif(n_iterations, min = 1, max = 60)
  
  G0_list <- runif(n_iterations, min = 1, max = 340)
  
  L0_list <- p$G_medium - N0_list - G0_list
  
  for (i in 1:n_iterations){
    
    x0_list <- c(N = N0_list[i], G = G0_list[i], L = L0_list[i])
    
    sol_list <- ode(x0_list,time,final_model,p)
    
    
    if (i == 1){
      
      output_final_list <- data.frame(sol_list)
      
      sim_number <- c(rep(i,length(time)))
      
    } else {
      
      next_sim <- rep(i,length(time))
      
      sim_number <- c(sim_number,next_sim)
      
      output_list <- data.frame(sol_list)
      
      output_final_list <- rbind(output_final_list,output_list)
      
    }
    
  }
  
  output_total <- cbind(output_final_list, sim_number)
  
  
  return(output_total)
  
}









