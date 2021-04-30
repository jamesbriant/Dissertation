####################
# This script implements rejection sampling algorithm for Beta(2, 5).
####################

####################
# Source files
####################
source("functions/Sampling.R")

####################
# Start of Code
####################

# Beta distribution PDF with default values alpha=2 & beta=5
BetaPDF <- function(x, alpha=2, beta=5){
  # Recall Beta distribution is defined only on [0,1]
  
  numerator <- x^(alpha-1) * (1-x)^(beta-1) * gamma(alpha + beta)
  denominator <- gamma(alpha) * gamma(beta)
  
  return(numerator/denominator)
}

n <- 10000

y <- runif(n, 0, 2.5)
x <- GetUniformSample(n, 0, 1)
fx <- BetaPDF(x)

rejectionSample <- x[y < fx]
importanceSample <- y

rejectionMean <- mean(rejectionSample)
rejectionMean
importanceMean <- sum(fx*x/1)/n
importanceMean

# E[Beta(2,5)] = 2/7 = 0.285714



#########################################################################

N <- 100000 # number of approximations to be made
n <- 500 # sample size (from proposal distribution)

rejection.means <- numeric(N)
importance.means <- numeric(N)

for(i in 1:N){
  y <- runif(n, 0, 2.5)
  x <- GetUniformSample(n, 0, 1)
  fx <- BetaPDF(x)
  
  rejectionSample <- x[y < fx]
  importanceSample <- y
  
  rejection.means[i] <- mean(rejectionSample)
  importance.means[i] <- sum(fx*x/1)/n
}

x <- seq(0.2, 0.4, length=200)
par(mfrow=c(1,2))
hist(rejection.means, freq=FALSE, xlab="Sample Means", main="Rejection Sampling Means")
lines(x, dnorm(x, mean(rejection.means), sqrt(var(rejection.means))), col="red")
hist(importance.means, freq=FALSE, xlab="Sample Means", main="Importance Sampling Means")
lines(x, dnorm(x, mean(importance.means), sqrt(var(importance.means))), col="red")
par(mfrow=c(1,1))


mean(rejection.means)
mean(importance.means)

sqrt(var(rejection.means))
sqrt(var(importance.means))














