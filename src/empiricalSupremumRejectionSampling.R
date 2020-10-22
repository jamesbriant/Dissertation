# Suppose we don't know the exact value of c, then we can estimate it.
# For example, we know the value of 2.5 is not most efficient from the previous example.
# Instead, we can estimate c to make the algorithm more efficient.
# We expect c to be just above 1.

source('functions/rejection-sampling-functions.R')

n <- 2500 # number of samples to draw

c_estimates <- numeric(n) # vector of estimates for c
c_estimates[1] <- 1.00001

y <- runif(n, 0, 2.3)
x <- drawFromg(n)
fx <- f(x)

rejectionSample <- x[y < fx]

for(i in 2:n){
  c_estimates[i] <- max(c_estimates[i-1], fx[i-1]/2.3)
}

plot(c_estimates)


hist(x[y*c_estimates < fx], freq=FALSE, xlab="x", main="Beta(2,5) Distribution from 10000 Draws")
x <- seq(from=0, to=1, length.out=250)
lines(x, dbeta(x, 2, 5), col="red")








