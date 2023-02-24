

rm(list=ls())

source("Models.R")

source("helper functions/3dplot.R")

source("helper functions/2dplot.R")

require(deSolve)

## Initial state
N0 <- 5
G0 <- 400
L0 <- 0

x0 <- c(N = N0, G = G0, L = L0)

p <- list()

p$rate <- 0.1

p$flow <- 0.35 # this shouldn't change from 0.75

p$G_medium <- 700

p$G50 <- 100

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



p <- c(rate = 0.1,
       flow = 0.75,
       G_medium = 400,
       G50 = 60,
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

f_result <- fit()

summary(f_result)
