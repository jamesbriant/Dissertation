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

####################
# Plot the functions
####################
x.axis <- seq(0,1, length=200)
plot(x.axis, dbeta(x.axis, 2, 5), type="l", xlab="x", ylab="Density", main="Distribution Comparison")
lines(x.axis, dunif(x.axis, 0, 1), col="blue")
legend("topright", 
       legend=c("Beta(2,5)",
                "Uniform[0,1]"),
       lty=c(1,1),
       col=c("black", "blue")
       )

plot(x.axis, dbeta(x.axis, 2, 5), type="l", xlab="x", ylab="Density", main="Envolope Function")
lines(x.axis, 2.4576*dunif(x.axis, 0, 1), col="red")
legend("right", 
       legend=c("Beta(2,5)",
                "Envelope Function"),
       lty=c(1,1),
       col=c("black", "red")
)

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




