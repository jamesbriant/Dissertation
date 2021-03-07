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