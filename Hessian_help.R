
library(numDeriv)

# Gamma-tæthedsfunktion med shape=2 og rate=5
par(mfrow=c(1,1))
xs <- seq(0,1,by=0.01)
ys <- dgamma(xs,2,5)
plot(xs,ys,type="l")

# Least squares funktion omkring gamma(2,5)
obj <- function(p){
  
  x <- seq(0,1,by=0.01)
  pred <- dgamma(x,p[1],p[2])
  y <- dgamma(x,2,5)
  
  sum((y - pred)^2)
  
}

# Evaluér least suqares i f.eks. de sande parametre vs. nogle andre
obj(c(2,5))
obj(c(1,5))

# Plot least squares skiftevis med de 2 parametre fixed
par(mfrow=c(1,2))
ps <- seq(0.1,10,by=0.1)
plot(ps,sapply(ps,function(x){obj(c(x,5))}),type="l")
plot(ps,sapply(ps,function(x){obj(c(2,x))}),type="l")

# Fisher Informations matrix
hessian(func = obj,x = c(2,5))

# Covariance matrix
solve(hessian(func = obj,x = c(2,5)))
