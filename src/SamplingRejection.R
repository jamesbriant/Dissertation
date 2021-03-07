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

# Use 1000 attempts
n <- 1000
c <- 2.5

# Run the algorithm
U <- runif(n, 0, 1)
X <- GetUniformSample(n, 0, 1)
Y <- X[U < BetaPDF(X)/c]

####################
# Plot Histogram
####################
hist(Y, freq=FALSE, xlab="x", main="Beta(2,5) Distribution from 1000 Draws of Unif(0,1)", xlim=c(0,1))
lines(seq(0, 1, length=500), BetaPDF(seq(0, 1, length=500)), col="red")
legend(0.6, 2.4, legend=c("Rejection Sampling", "Beta(2,5)"), pch="-", col=c("black", "red"))
rug(Y)




