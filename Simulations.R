
rm(list=ls())

install.packages("vctrs")

#require(deSolve)

#require(ggplot2)

#require(tidyr)

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


## Initial state
N0 <- 1
G0 <- 100

x0 <- c(N = N0, G = G0)

p <- list()

p$rate <- 0.06

p$flow <- 0.75 # this shouldn't change from 0.75

p$G_medium <- 400


time <- seq(0,15,1)

sol <- ode(x0,time,base_model,p)

output <- data.frame(sol)



#### plotting section ####



ggplot(data = output, aes(x = time, y = N)) + geom_point() + geom_line(color = "red") + 
  labs(title = "Number of cells", x = "time", y = "number of cells")

ggplot(data = output, aes(x = time, y = G)) + geom_point(aes(color = "Glucose levels")) + geom_line() + 
  labs(title = "Glucose levels", x = "time", y = "number of cells") + 
  ylim(0, max(output$G)+5)










