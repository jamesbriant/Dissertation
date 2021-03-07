####################
# This script implements rejection sampling algorithm for Beta(2, 5).
# The improvement guarantees n samples from the target distribution.
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

# Generate 1000 samples
n <- 1000


c <- 2.5

BetaSample <- numeric(n)


# Run the algorithm
i <- 1
while(i <= n){
  U <- SimulateU()
  X <- GetUniformSample(1, 0, 1)
  Y <- BetaPDF(X)/c

  if(U < Y){
    BetaSample[i] <- X
    i <- i + 1
  }  
}

####################
# Plot Histogram
####################

hist(BetaSample, freq=FALSE, xlab="x", main="1000 Samples from Beta(2,5)", xlim=c(0,1))
lines(seq(0, 1, length=500), BetaPDF(seq(0, 1, length=500)), col="red")
legend(0.6, 2.4, legend=c("Rejection Sampling", "Beta(2,5)"), pch="-", col=c("black", "red"))
rug(BetaSample)




