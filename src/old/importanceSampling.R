source('functions/rejection-sampling-functions.R')

n <- 10000

y <- runif(n, 0, 2.5)
x <- drawFromg(n)
fx <- f(x)

rejectionSample <- x[y < fx]
importanceSample <- y

rejectionMean <- mean(rejectionSample)
rejectionMean
importanceMean <- sum(fx*x/1)/n
importanceMean

# E[Beta(2,5)] = 2/7 = 0.285714


