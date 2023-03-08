

rm(list=ls())

source("Phase_analysis_func.R")


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

p$n_rate_inhib_max <- 0.8

p$lac_con_max <- 0.9

p$lac_prod_max <- 0.75

time <- seq(0,30,0.1)

sol <- ode(x0,time,final_model,p)

output <- data.frame(sol)


N0_list <- runif(1000, 0, 200)

G0_list <- runif(1000, 0, 350)

L0_list <- p$G_medium - N0_list - G0_list






