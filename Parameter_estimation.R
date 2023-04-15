

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
               lac_prod_mid = 50)

sol_real <- ode(init,time,final_model_estimation,new_param)

# parameters to fit 

param_to_fit <- c(rate = 0.199,
                  G50 = 30)

# the minimum and maximum range in which the fit takes place

bound_min_var <- c(0.01,
                   10)

bound_max_var <- c(0.8,
                   100)


results <- fit_simulation(sol_real, param_to_fit, bound_min_var, bound_max_var,"L-BFGS-B", 100)

# change this so we have a vector of the parameters and a column that explains which parameter is what.

param_list <- results[[1]]

error_sim <- results[[2]]

## plotting ## 

histogram_sim_maker(param_list,c(0.6,60))


noise_test_data <- add_noise(sol_real, mean = 0, sd = 1)


#### TESTING IF OUR ALGORITHM WORKS AND OTHER STUFF ####

test_results <- param_initial_state_tester(noise_test_data,c("rate","G50"),bound_min_var, bound_max_var,"L-BFGS-B",50)

# plotting the results

ggplot(data = test_results, aes(x = rate, y = G50))+
  geom_point(aes(colour = state), size = 5)+
  geom_bin_2d(bins = 35)+
  scale_fill_continuous(low = "dark green", high = "blue")+
  ggtitle("Initial values and their end destination")



error_function(new_param, c(rate = 0.6, G50 = 60), real_data = sol_real)



### Constructing the contour / grid stuff ###

param1 <- seq(bound_min_var[1],bound_max_var[1],0.01)

names(param1) <- "rate"

#rep("rate", length(param1))


param2 <- seq(bound_min_var[2],bound_max_var[2], length = 80)

names(param2) <- "G50"


# has to call it real_data


error_function(new_param,c(rate = 0.6, G50 = 60), real_data = sol_real)


results_contour <- param_plot_contour(new_param,
                                     c("rate","G50"),
                                     bound_min_var,
                                     bound_max_var,
                                     c(100,100),
                                     true_data = sol_real)

parameter_1 <- results_contour[[1]]

parameter_2 <- results_contour[[2]]

matrix <- results_contour[[3]]

image(parameter_1, parameter_2, matrix)


fig <- plot_ly(
  x = parameter_1, 
  y = parameter_2, 
  z = matrix,
  type = "contour",
  autocontour = F,
  contours = list(
    start = 0,
    end = 10000,
    size = 100
  ),
  line = list(smoothing = 2)
)

fig

# which(matrix == min(matrix), arr.ind = TRUE)



