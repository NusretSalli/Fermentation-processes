

rm(list = ls())

require(deSolve)

require(ggplot2)

require(tidyr)

require(phaseR)

require(ODEsensitivity)

require(epiR)

require(FME)

require(dplyr)

require(tidyverse)

require(numDeriv)

require(plotly)

require(minpack.lm)

require(numDeriv)

source("Models.R")

source("Parameter_estimation_func.R")


### Parameter estimation ###

parameters <- c(rate = 0.199,
                flow = 0.75,
                G_medium = 800,
                G50 = 50,
                N_rate_inhib_growth = 0.2,
                lac_con_growth = 0.2,
                lac_prod_growth = 0.2,
                N_rate_inhib_mid = 50,
                lac_con_mid = 60,
                lac_prod_mid = 50)

N0 <- 5
G0 <- 400
L0 <- 0

init <- c(N = N0, G = G0, L = L0)

time <- seq(0,30,0.1)

new_param <- c(rate = 0.6,
               flow = 0.75,
               G_medium = 600,
               G50 = 30,
               N_rate_inhib_growth = 0.2,
               lac_con_growth = 0.2,
               lac_prod_growth = 0.2,
               N_rate_inhib_mid = 50,
               lac_con_mid = 60,
               lac_prod_mid = 50)

sol_real <- ode(init,time,final_model_estimation,new_param)

noised_data <- add_noise(sol_real, mean = 10, sd = 6)

noised_data_frame <- data.frame(noised_data)

ggplot(data = noised_data_frame, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) + 
  labs(title = "Number of cells", x = "time", y = "number of cells")


# noise <- as.data.frame(matrix(rnorm(dim(output_real)[1]*(dim(output_real)[2]-1), mean = 2, sd = 4), ncol = dim(output_real)[2]-1))
# 
# 
# what <- output_real[,2:dim(output_real)[2]] + noise


# adding noise to our data with rnorm with a self-chosen mean and sd value


# function that solves the model given the parameters as input

solve_model <- function(pars) {
  
  # initial values
  
  state <- c(N = N0, G = G0, L = L0)
  
  # time vector
  
  time <- seq(0,30,0.1)
  
  # store the output and returning it
  
  output <- ode(y = state, time, final_model_estimation, parms = pars)
  
  return(output)
}



# our objective function given the parameters

objective_func <- function(x, parset = names(x)) {
  
  # we substitute the parameter values with x ()
  
  parameters[parset] <- x
  
  output <- solve_model(parameters) # solving the model
  
  modCost(output, noised_data, x = "time") # calculating the loss
  
}

# parameters to fit 

param_to_fit <- c(rate = 0.7,
                  G_medium = 60,
                  G50 = 50)

# the minimum and maximum range in which the fit takes place

bound_min_var <- c(0.001,
                   30,
                   5)

bound_max_var <- c(0.7,
                   700,
                   100)

# fitting the model

# good optim algorithms -> bobyqa, Pseudo, L-BFGS-B

fit <- modFit(objective_func,
              p = param_to_fit,
              lower = bound_min_var,
              upper = bound_max_var,
              method = "L-BFGS-B")

fit

parameters[c(1,3,4)] = fit$par 

fit_output <- ode(init,time,final_model_estimation,parameters)




#################### parameter estimation simulation ##################

# we try rate and g50

parameters <- c(rate = 0.199,
                flow = 0.75,
                G_medium = 600,
                G50 = 30,
                N_rate_inhib_growth = 0.2,
                lac_con_growth = 0.2,
                lac_prod_growth = 0.2,
                N_rate_inhib_mid = 50,
                lac_con_mid = 60,
                lac_prod_mid = 50)

N0 <- 5
G0 <- 400
L0 <- 0

init <- c(N = N0, G = G0, L = L0)

time <- seq(0,30,0.1)

new_param <- c(rate = 0.6,
               flow = 0.75,
               G_medium = 600,
               G50 = 60,
               N_rate_inhib_growth = 0.2,
               lac_con_growth = 0.2,
               lac_prod_growth = 0.2,
               N_rate_inhib_mid = 50,
               lac_con_mid = 60,
               lac_prod_mid = 50,
               N_rate_inhib_max = 0.9,
               lac_con_max = 0.9,
               lac_prod_max = 0.9)

sol_real <- ode(init,time,final_model_estimation,new_param)

# parameters to fit 

param_to_fit <- c(rate = 0.199,
                  G50 = 30)

# the minimum and maximum range in which the fit takes place

bound_min_var <- c(0.01,
                   10)

bound_max_var <- c(0.8,
                   100)


results <- fit_simulation(sol_real, param_to_fit, bound_min_var, bound_max_var,"L-BFGS-B", 1000)

# change this so we have a vector of the parameters and a column that explains which parameter is what.

param_list <- results[[1]]

error_sim <- results[[2]]


## plotting ## 

histogram_sim_maker(param_list,c(0.6,60))


noise_test_data <- add_noise(sol_real, mean_val = 0, sd_val = 1)


#### TESTING IF OUR ALGORITHM WORKS AND OTHER STUFF ####

test_results <- param_initial_state_tester(noise_test_data,c("rate","G50"),bound_min_var, bound_max_var,"L-BFGS-B",100)

# only for plotting 

number_vector <- rep(NA,200)

number_vector[c(66,66+100)] <- 1

number_vector[101] <- 99

final_test_results <- cbind(test_results,number_vector)

# plotting the results

ggplot(data = final_test_results, aes(x = rate, y = G50, label = number_vector))+
  geom_point(aes(colour = state), size = 4)+
  #geom_bin_2d(bins = 60)+
  geom_text(size = 4)+
  scale_fill_continuous(low = "black", high = "red")+
  ggtitle("Initial values and their end destination")


### Constructing the contour / grid stuff ###

param1 <- seq(bound_min_var[1],bound_max_var[1],0.01)

names(param1) <- "rate"

#rep("rate", length(param1))


param2 <- seq(bound_min_var[2],bound_max_var[2], length = 80)

names(param2) <- "G50"


error_function(new_param,minima_point, real_data = sol_real)


results_contour <- param_plot_contour(new_param,
                                     c("rate","G50"),
                                     c(0.3,35), # min range
                                     c(0.8,85), # max range
                                     c(500,500),
                                     true_data = sol_real)

results_contour_noise <- param_plot_contour(new_param,
                                            c("rate","G50"),
                                            c(0.3,35), # min range
                                            c(0.8,85), # max range
                                            c(300,300),
                                            true_data = noise_test_data)

parameter_1 <- results_contour[[1]]

parameter_1_noise <- results_contour_noise[[1]]

parameter_2 <- results_contour[[2]]

parameter_2_noise <- results_contour_noise[[2]]

matrix <- results_contour[[3]]

matrix_noise <- results_contour_noise[[3]]

minima_point <- c(parameter_1[which(matrix == min(matrix), arr.ind = TRUE)[2]],
                  parameter_2[which(matrix == min(matrix), arr.ind = TRUE)[1]])

minima_point_noise <- c(parameter_1_noise[which(matrix_noise == min(matrix_noise), arr.ind = TRUE)[2]],
                        parameter_2_noise[which(matrix_noise == min(matrix_noise), arr.ind = TRUE)[1]])

fig <- plot_ly(
  x = parameter_1, 
  y = parameter_2, 
  z = matrix,
  type = "contour",
  #colors = colorRamp(c("dark green", "blue")),
  autocontour = F,
  contours = list(
    start = 0,
    end = 1000,
    size = 10
  ),
  line = list(smoothing = 1),
  colorscale = "Jet"
  
)

fig %>% add_trace(x = ~minima_point[1], y = ~minima_point[2], name = paste0("minimum: ",as.character(min(matrix))), type = "scatter")


fig_noise <- plot_ly(
  x = parameter_1_noise, 
  y = parameter_2_noise, 
  z = matrix_noise,
  type = "contour",
  #colors = colorRamp(c("dark green", "blue")),
  autocontour = F,
  contours = list(
    start = 800,
    end = 1800,
    size = 10
  ),
  line = list(smoothing = 1),
  colorscale = "Jet"
  
)

fig_noise %>% add_trace(x = ~minima_point_noise[1], y = ~minima_point_noise[2], name = paste0("minimum: ",as.character(min(matrix))), type = "scatter")


# which(matrix == min(matrix), arr.ind = TRUE)

### GGPLOT SECTION ###

# rownames(matrix) <- parameter_1
# 
# colnames(matrix) <- parameter_2
# 
# as.data.frame(matrix) %>% 
#   rownames_to_column() %>% 
#   gather(key, value, -rowname) %>% 
#   mutate(key = as.numeric(key), 
#          rowname = as.numeric(rowname)) %>%
#   ggplot() +
#   #geom_contour(aes(x = rowname, y = key, z = value, colour = value), breaks = seq(0,10000,100))+
#   geom_contour_filled(aes(x = rowname, y = key, z = value), alpha = 0.6, breaks = c(0,10,100,1000, 10000,1000000,10000000000))+
#   geom_point(aes(x = (as.numeric(rownames(matrix)[(which(matrix == min(matrix), arr.ind = TRUE)[2])])),
#                  y = (as.numeric(colnames(matrix)[(which(matrix == min(matrix), arr.ind = TRUE)[1])]))),
#                  color = "black", size = 2)+
#   geom_point(data = test_results, aes(x = rate, y = G50, colour = state), size = 3)+
#   #geom_bin_2d(bins = 60)+
#   #geom_jitter(width = 0.01)+
#   #scale_fill_continuous(low = "black", high = "red")+
#   ggtitle("Initial values and their end destination")
#   

## Calculating the covariance and correlation matrix for param_list

cov_emp_matrix <- cov(param_list)

corr_emp_matrix <- cor(param_list)



## calculating the covariance and correlation matrix with numderiv ##

error_function_hessian <- function(x){
  
  parset = names(x)
  
  new_param[parset] <- x
  
  output <- solve_model(new_param)
  
  total_err <- 1/2*sum((noise_test_data - output)^2)
  
  return(total_err)
  
}

calc_hessian_matrix <- hessian(func = error_function_hessian, x = c(minima_point_noise[1],
                                                                    minima_point_noise[2]))

solved_hessian_matrix <- solve(calc_hessian_matrix)

#equal_hessian_matrix <-  2.691218 * solved_hessian_matrix

#corr_theory_matrix <- cov2cor(equal_hessian_matrix) 



### PARAMETER ESTIMATION WITH lactate switches ### 

parameters_lactate <- c(rate = 0.04,
                flow = 0.75,
                G_medium = 400,
                G50 = 50,
                N_rate_inhib_growth = 0.5,
                lac_con_growth = 0.5,
                lac_prod_growth = 4,
                N_rate_inhib_mid = 120,
                lac_con_mid = 20,
                lac_prod_mid = 40,
                N_rate_inhib_max = 0.9,
                lac_con_max = 0.9,
                lac_prod_max = 0.15)

N0 <- 5
G0 <- 395
L0 <- 0

init <- c(N = N0, G = G0, L = L0)

time <- seq(0,30,0.1)

new_param_lactate <- c(rate = 0.04,
                        flow = 0.75,
                        G_medium = 400,
                        G50 = 50,
                        N_rate_inhib_growth = 0.5,
                        lac_con_growth = 0.5,
                        lac_prod_growth = 0.5,
                        N_rate_inhib_mid = 120, # change to 90 for oscillation 120 for non-oscillation
                        lac_con_mid = 20,
                        lac_prod_mid = 120,
                        N_rate_inhib_max = 0.9,
                        lac_con_max = 0.9,
                        lac_prod_max = 0.9)


sol_real <- ode(init,time,final_model_estimation_lactate_switch,new_param_lactate)

sol_real_dataframe <- data.frame(sol_real)

ggplot(data = sol_real_dataframe, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = sol_real_dataframe, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Glucose levels", x = "time", y = "glucose levels") +
  ylim(0, max(sol_real_dataframe$G)+5)

ggplot(data = sol_real_dataframe, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Lactate levels", x = "time", y = "Lactate") +
  ylim(0, max(sol_real_dataframe$L)+5)

param_to_fit_lactate <- c(lac_prod_growth = 4,
                          lac_prod_mid = 40,
                          lac_prod_max = 0.15)

bound_min_var_lactate <- c(0.1,
                           30,
                           0.1)

bound_max_var_lactate <- c(5,
                           170,
                           1)

results_lactate <- fit_simulation_lactate(sol_real,
                          param_to_fit_lactate,
                          bound_min_var_lactate,
                          bound_max_var_lactate,
                          "L-BFGS-B", # L-BFGS-B
                          1000)


param_list_lactate <- results_lactate[[1]]

error_sim_lactate <- results_lactate[[2]]


histogram_sim_maker(param_list_lactate,
                    c(lac_prod_growth = 0.5,
                      lac_prod_mid = 120,
                      lac_prod_max = 0.9))


cov_lactate <- cov(param_list_lactate)

corr_lactate <- cor(param_list_lactate)


solve_model <- function(pars) {
  
  # initial values
  
  state <- c(N = N0, G = G0, L = L0)
  
  # time vector
  
  time <- seq(0,30,0.1)
  
  # store the output and returning it
  
  output <- ode(y = state, time, final_model_estimation_lactate_switch, parms = pars)
  
  return(output)
}


error_function_hessian_lactate <- function(x){
  
  parset = names(x)
  
  new_param_lactate[parset] <- x
  
  output <- solve_model(new_param_lactate)
  
  total_err <- 1/2*sum((sol_real - output)^2)
  
  return(total_err)
  
}

error_function_hessian_lactate(c(lac_prod_growth = 0.5,
                                 lac_prod_mid = 120,
                                 lac_prod_max = 0.9))

calc_hessian_matrix <- hessian(func = error_function_hessian_lactate, c(lac_prod_growth = 0.5,
                                                                lac_prod_mid = 120,
                                                                lac_prod_max = 0.9))

solved_hessian_matrix <- solve(calc_hessian_matrix)



## another example with lac_con ##

parameters <- c(rate = 0.04,
                        flow = 0.75,
                        G_medium = 400,
                        G50 = 50,
                        N_rate_inhib_growth = 0.5,
                        lac_con_growth = 4,
                        lac_prod_growth = 0.5,
                        N_rate_inhib_mid = 120,
                        lac_con_mid = 210,
                        lac_prod_mid = 120,
                        N_rate_inhib_max = 0.9,
                        lac_con_max = 0.3,
                        lac_prod_max = 0.9)

N0 <- 5
G0 <- 395
L0 <- 0

init <- c(N = N0, G = G0, L = L0)

time <- seq(0,30,0.1)

new_param_lactate <- c(rate = 0.04,
                       flow = 0.75,
                       G_medium = 400,
                       G50 = 50,
                       N_rate_inhib_growth = 0.5,
                       lac_con_growth = 0.5,
                       lac_prod_growth = 0.5,
                       N_rate_inhib_mid = 120, # change to 90 for oscillation 120 for non-oscillation
                       lac_con_mid = 20,
                       lac_prod_mid = 120,
                       N_rate_inhib_max = 0.9,
                       lac_con_max = 0.9,
                       lac_prod_max = 0.9)


sol_real <- ode(init,time,final_model_estimation_lactate_switch,new_param_lactate)

sol_real_dataframe <- data.frame(sol_real)

ggplot(data = sol_real_dataframe, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = sol_real_dataframe, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Glucose levels", x = "time", y = "glucose levels") +
  ylim(0, max(sol_real_dataframe$G)+5)

ggplot(data = sol_real_dataframe, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Lactate levels", x = "time", y = "Lactate") +
  ylim(0, max(sol_real_dataframe$L)+5)

param_to_fit_total <- c(lac_con_growth = 4,
                          lac_con_mid = 210,
                          lac_con_max = 0.3,
                          lac_prod_growth = 4,
                          lac_prod_mid = 40,
                          lac_prod_max = 0.15)

bound_min_var_total <- c(0.1,
                           1,
                         0.01,
                         0.1,
                         30,
                         0.1)

bound_max_var_total <- c(5,
                           230,
                           1,
                         5,
                         170,
                         1)

results_lactate <- fit_simulation_lactate(sol_real_total,
                                          param_to_fit_total,
                                          bound_min_var_total,
                                          bound_max_var_total,
                                          "L-BFGS-B", # L-BFGS-B
                                          1000)


param_list_lactate <- results_lactate[[1]]

error_sim_lactate <- results_lactate[[2]]


histogram_sim_maker(param_list_lactate,
                    c(lac_con_growth = 0.5,
                      lac_con_mid = 20,
                      lac_con_max = 0.9))


cov_lactate2 <- cov(param_list_lactate)

corr_lactate2 <- cor(param_list_lactate)


error_function_hessian_lactate <- function(x){
  
  parset = names(x)
  
  new_param_lactate[parset] <- x
  
  output <- solve_model(new_param_lactate)
  
  total_err <- 1/2*sum((sol_real - output)^2)
  
  return(total_err)
  
}

error_function_hessian_lactate(c(lac_con_growth = 0.5,
                                 lac_con_mid = 20,
                                 lac_con_max = 0.9))

calc_hessian_matrix <- hessian(func = error_function_hessian_lactate, c(lac_con_growth = 0.5,
                                                                        lac_con_mid = 20,
                                                                        lac_con_max = 0.9))

solved_hessian_matrix <- solve(calc_hessian_matrix)




### total lactate switch estimation

parameters <- c(rate = 0.04,
                flow = 0.75,
                G_medium = 400,
                G50 = 50,
                N_rate_inhib_growth = 0.5,
                lac_con_growth = 4,
                lac_prod_growth = 4,
                N_rate_inhib_mid = 120,
                lac_con_mid = 210,
                lac_prod_mid = 40,
                N_rate_inhib_max = 0.9,
                lac_con_max = 0.3,
                lac_prod_max = 0.15)

N0 <- 5
G0 <- 395
L0 <- 0

init <- c(N = N0, G = G0, L = L0)

time <- seq(0,30,0.1)

new_param_total <- c(rate = 0.04,
                       flow = 0.75,
                       G_medium = 400,
                       G50 = 50,
                       N_rate_inhib_growth = 0.5,
                       lac_con_growth = 0.5,
                       lac_prod_growth = 0.5,
                       N_rate_inhib_mid = 120, # change to 90 for oscillation 120 for non-oscillation
                       lac_con_mid = 20,
                       lac_prod_mid = 120,
                       N_rate_inhib_max = 0.9,
                       lac_con_max = 0.9,
                       lac_prod_max = 0.9)


sol_real <- ode(init,time,final_model_estimation_lactate_switch,new_param_total)

sol_real_dataframe <- data.frame(sol_real)

ggplot(data = sol_real_dataframe, aes(x = time, y = N)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = sol_real_dataframe, aes(x = time, y = G)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Glucose levels", x = "time", y = "glucose levels") +
  ylim(0, max(sol_real_dataframe$G)+5)

ggplot(data = sol_real_dataframe, aes(x = time, y = L)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Lactate levels", x = "time", y = "Lactate") +
  ylim(0, max(sol_real_dataframe$L)+5)

param_to_fit_total <- c(lac_con_growth = 4,
                        lac_con_mid = 210,
                        lac_con_max = 0.3,
                        lac_prod_growth = 4,
                        lac_prod_mid = 40,
                        lac_prod_max = 0.15)

bound_min_var_total <- c(0.1,
                         1,
                         0.01,
                         0.1,
                         30,
                         0.1)

bound_max_var_total <- c(5,
                         230,
                         1,
                         5,
                         170,
                         1)

results_total <- fit_simulation_lactate(sol_real,
                                          param_to_fit_total,
                                          bound_min_var_total,
                                          bound_max_var_total,
                                          "L-BFGS-B", # L-BFGS-B
                                          500)


param_list_total <- results_total[[1]]

error_sim_total <- results_total[[2]]


histogram_sim_maker(param_list_total,
                    c(lac_con_growth = 0.5,
                      lac_con_mid = 20,
                      lac_con_max = 0.9,
                      lac_prod_growth = 0.5,
                      lac_prod_mid = 120,
                      lac_prod_max = 0.9))


cov_lactate_total <- cov(param_list_total)

corr_lactate_total <- cor(param_list_total)



## logistic function plots.

G_level <- seq(100,130,1)

lac_con_vec <- (1 / (1 + exp(p$lac_con_growth*(G_level-p$lac_con_mid))))

lac_prod_vec <- (1 / (1 + exp(-p$lac_prod_growth*(G_level-p$lac_prod_mid))))

logistic_data_frame <- data.frame(cbind(lac_con_vec, lac_prod_vec, G_level))


ggplot(data = logistic_data_frame, aes(x = G_level, y = lac_con_vec)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Logistic function for rho(G)", x = "Glucose levels", y = "factor")

ggplot(data = logistic_data_frame, aes(x = G_level, y = lac_prod_vec)) + geom_point(size = 3, color = "blue") + geom_line(color = "red", linewidth = 1.5) +
  labs(title = "Logistic function for xi(G)", x = "Glucose levels", y = "factor")





