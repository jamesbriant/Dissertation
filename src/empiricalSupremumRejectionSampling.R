# Suppose there exists some x in the support of g such that f(x) > g(x) 
# For example, we know the value of 2.3 < sup{f(x) : x in support of f}
# Suppose we don't know the exact value of c, then we can estimate it.
# We expect c to be just above 1.

source('functions/rejection-sampling-functions.R')

n <- 5000 # number of samples to draw

c_estimates <- numeric(n) # vector of estimates for c
c_estimates[1] <- 1.00001

y <- runif(n, 0, 2.3)
x <- drawFromg(n)
fx <- f(x)

for(i in 2:n){
  c_estimates[i] <- max(c_estimates[i-1], fx[i-1]/2.3)
}

plot(c_estimates)
# From this plot it can be seen c usually converges by around 200 iterations

hist(x[y*c_estimates < fx], freq=FALSE, xlab="x", main="Beta(2,5) Distribution from 5000 Draws")
x <- seq(from=0, to=1, length.out=250)
lines(x, dbeta(x, 2, 5), col="red")








