
library(docstring)

base_model <- function(t,x,p)
{
  
  #' function that computes dn and dg with our base model.
  #' 
  #' input:
  #' 
  #' t -> time vector 
  #' x -> essentially our states
  #' p -> parameters used in the ODEs
  #' 
  #' 
  #' Uses:
  #' 
  #' used in Desolve's ode-function to solve these ODEs numerically.
  #' 
  #' Example:
  #' 
  #' num_solution <- ode(initial,time_interval, base_model, parameters)
  
  
  
  ## Unpack state by hand 
  n <- length(x) / 2
  N <- x[1:n]
  G <- x[n + (1:n)]
  
  # defines the derivatives as numeric data structures
  
  dN <- numeric(n)
  dG <- numeric(n)
  
  # defining the ODEs with our parameters.
  
  with(p,
       {
         dN <- rate * N * G - flow * N
         dG <- flow * (G_medium - G) - rate * N * G
         
         return(list(c(dN, dG)))
       }
  )
  
}




lactate_model <- function(t,x,p)
{
  
  #' function that computes dN, dG and dL with our lactate model.
  #' 
  #' input:
  #' 
  #' t -> time vector 
  #' x -> essentially our states
  #' p -> parameters used in the ODEs
  #' 
  #' 
  #' Uses:
  #' 
  #' used in Desolve's ode-function to solve these ODEs numerically.
  #' 
  #' Example:
  #' 
  #' num_solution <- ode(initial,time_interval, lactate_model, parameters)
  
  ## Unpack state by hand 
  n <- length(x) / 3
  N <- x[1:n]
  G <- x[n + (1:n)]
  L <- x[2*n + (1:n)]
  
  # defines the derivatives as numeric data structures
  
  dN <- numeric(n)
  dG <- numeric(n)
  dL <- numeric(n)
  
  # defining the ODEs with our parameters.
  
  with(p,
       {
         
         logL <- (1 / (1 + exp(L_log_growth*(G-L_log_mid))))
         
         logN <- (0.7 / (1 + exp(N_log_growth*(L-N_log_mid)))) + 0.3
         
         # maybe logN shouldn't be here.
         
         # G should maybe be "mætning" - maximum number of rate?
         
         dN <- rate * N * G * logN - flow * N
         
         dG <- flow * (G_medium - G) - rate * N * G + L * N * logL
         
         # introduce "switch" in lactate production. 
         
         dL <- G * N * l_rate - flow * L - L * N * logL
         
         return(list(c(dN, dG, dL), c(logisticN = logN, logisticL = logL)))
       }
  )
  
}

final_model <- function(t,x,p)
{
  #' function that computes dN, dG, dL and dP with our final model.
  #' 
  #' input:
  #' 
  #' t -> time vector 
  #' x -> essentially our states
  #' p -> parameters used in the ODEs
  #' 
  #' 
  #' Uses:
  #' 
  #' used in Desolve's ode-function to solve these ODEs numerically.
  #' 
  #' Example:
  #' 
  #' num_solution <- ode(initial,time_interval, final_model, parameters)
  
  
  ## Unpack state by hand 
  n <- length(x) / 3
  N <- x[1:n]
  G <- x[n + (1:n)]
  L <- x[2*n + (1:n)]
  
  # Defines the derivatives as numeric data structures
  
  dN <- numeric(n)
  dG <- numeric(n)
  dL <- numeric(n)
  
  
  # defining the ODEs with our parameters.
  
  with(p,
       {
         
         #logL <- (1 / (1 + exp(L_log_growth*(G-L_log_mid))))
         
         #logN <- (0.7 / (1 + exp(N_log_growth*(L-N_log_mid)))) + 0.3
         
         # maybe logN shouldn't be here.
         
         # G should maybe be "mætning" - maximum number of rate?
         
         N_rate_inhib <- (n_rate_inhib_max / (1 + exp(N_rate_inhib_growth*(L-N_rate_inhib_mid))))
         
         lac_con <- (lac_con_max / (1 + exp(lac_con_growth*(G-lac_con_mid))))
         
         lac_prod <- (lac_prod_max / (1 + exp(-lac_prod_growth*(G-lac_prod_mid))))
         
         
         dN <- rate * N * G / (1 + G / G50) * N_rate_inhib + N * L * lac_con - flow * N
         
         dG <- - rate * N * G / (1 + G / G50) * N_rate_inhib + flow * (G_medium - G) - N * lac_prod
         
         dL <- N * lac_prod - flow * L - N * L * lac_con
         
         
         return(list(c(dN, dG, dL), c(N_rate_inhib = N_rate_inhib, lac_con = lac_con, lac_prod = lac_prod)))
       }
  )
  
}




