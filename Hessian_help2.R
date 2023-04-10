
X <- rnorm(n=1000,mean=1,sd=2)

loss <- function(ms)0.5*sum((X-ms[1])^2)/ms[2]^2 + log(ms[2])*length(X)

mv <- seq(0,2,length=101)
sv <- seq(1,3,length=51)
LL <- outer(mv,sv,Vectorize(function(m,s)loss(c(m,s))))

contour(mv,sv,LL)

require(numDeriv)

I <- numDeriv::hessian(loss,c(1,2))

V <- solve(I)

fit <- lm(X~1)

summary(fit)


