

rm(list=ls())

source("Models.R")

source("helper functions/3dplot.R")

source("helper functions/2dplot.R")

source("Sensitivity_analysis.R")

source("Parameter_estimation.R")

require(deSolve)

require(FME)

## Initial state
N0 <- 5
G0 <- 400
L0 <- 0

x0 <- c(N = N0, G = G0, L = L0)

p <- list()

p$rate <- 0.1

p$flow <- 0.75 # this shouldn't change from 0.75

p$G_medium <- 400

p$G50 <- 50

p$N_rate_inhib_growth <- 0.2

p$lac_con_growth <- 0.2

p$lac_prod_growth <- 0.2

p$N_rate_inhib_mid <- 50

p$lac_con_mid <- 60

p$lac_prod_mid <- 30

time <- seq(0,30,0.1)

sol <- ode(x0,time,final_model,p)

output <- data.frame(sol)


model <- function(t,x,p){
  
  with(as.list(c(x,p)), {
    
    N_rate_inhib <- (1 / (1 + exp(N_rate_inhib_growth*(L-N_rate_inhib_mid))))
    
    lac_con <- (1 / (1 + exp(lac_con_growth*(G-lac_con_mid))))
    
    lac_prod <- (1 / (1 + exp(-lac_prod_growth*(G-lac_prod_mid))))
    
    
    
    dN <- rate * N * G / (1 + G / G50) * N_rate_inhib + N * L * lac_con - flow * N
    
    dG <- - rate * N * G / (1 + G / G50) * N_rate_inhib + flow * (G_medium - G) - N * lac_prod
    
    dL <- N * lac_prod - flow * L - N * L * lac_con
    
    
    return(list(c(dN, dG, dL), c(N_rate_inhib = N_rate_inhib, lac_con = lac_con, lac_prod = lac_prod)))
    
    
  })
  
}



p <- c(rate = 0.8,
       flow = 0.25,
       G_medium = 800,
       G50 = 20,
       N_rate_inhib_growth = 0.2,
       lac_con_growth = 0.2,
       lac_prod_growth = 0.2,
       N_rate_inhib_mid = 50,
       lac_con_mid = 60,
       lac_prod_mid = 50)

N0 <- 5
G0 <- 400
L0 <- 0

s <- c(N = N0, G = G0, L = L0)

data <- output[1:4]

fixed = list(N = N0, G = G0, L = L0)

f_result <- fit()

summary(f_result)




### parameter estimation using modfit()

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

output_est <- ode(init,time,final_model_estimation,parameters)


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

output_real <- ode(init,time,final_model_estimation,new_param)


modcost_test <- modCost(output_est,output_real, x = "time")


error_functional_test <- function(param){
  
  est <- ode(init,time,final_model_estimation,param)
  
  return(modCost(est,output_real, x = "time"))
  
  
}

bound_min_var <- c(0.001,
                   0.20,
                   50,
                   10,
                   0.1,
                   0.1,
                   0.1,
                   10,
                   10,
                   10)

bound_max_var <- c(0.7,
                   0.95,
                   1000,
                   400,
                   2,
                   2,
                   2,
                   120,
                   120,
                   120)


fit <- modFit(error_functional_test, p = parameters, lower = bound_min_var, upper = bound_max_var)


output_estimated <- ode(init, time, final_model_estimation,fit$par)



