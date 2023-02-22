


# The base model which only have G and N as differential equations

base_model <- function(t,x,p)
{
  ## Unpack state by hand 
  n <- length(x) / 2
  N <- x[1:n]
  G <- x[n + (1:n)]
  
  dN <- numeric(n)
  dG <- numeric(n)
  
  with(p,
       {
         dN <- rate * N * G - flow * N
         dG <- flow * (G_medium - G) - rate * N * G
         
         return(list(c(dN, dG)))
       }
  )
  
}


# The lactate model which consists of system of differential equations of N, G and L

lactate_model <- function(t,x,p)
{
  ## Unpack state by hand 
  n <- length(x) / 3
  N <- x[1:n]
  G <- x[n + (1:n)]
  L <- x[2*n + (1:n)]
  
  dN <- numeric(n)
  dG <- numeric(n)
  dL <- numeric(n)
  
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
  ## Unpack state by hand 
  n <- length(x) / 3
  N <- x[1:n]
  G <- x[n + (1:n)]
  L <- x[2*n + (1:n)]
  
  dN <- numeric(n)
  dG <- numeric(n)
  dL <- numeric(n)
  
  with(p,
       {
         
         #logL <- (1 / (1 + exp(L_log_growth*(G-L_log_mid))))
         
         #logN <- (0.7 / (1 + exp(N_log_growth*(L-N_log_mid)))) + 0.3
         
         # maybe logN shouldn't be here.
         
         # G should maybe be "mætning" - maximum number of rate?
         
         N_rate_inhib <- (1 / (1 + exp(N_rate_inhib_growth*(L-N_rate_inhib_mid))))
         
         lac_con <- (1 / (1 + exp(lac_con_growth*(G-lac_con_mid))))
         
         lac_prod <- (1 / (1 + exp(-lac_prod_growth*(G-lac_prod_mid))))
         
         
         
         dN <- rate * N * G / (1 + G / G50) * N_rate_inhib + N * L * lac_con - flow * N
         
         dG <- - rate * N * G / (1 + G / G50) * N_rate_inhib + flow * (G_medium - G) - N * lac_prod
         
         dL <- N * lac_prod - flow * L - N * L * lac_con
         
         
         return(list(c(dN, dG, dL), c(N_rate_inhib = N_rate_inhib, lac_con = lac_con, lac_prod = lac_prod)))
       }
  )
  
}




