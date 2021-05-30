####################
# This script implements empirical supremum rejection sampling algorithm for 
# Beta(2, 5).
####################

####################
# Source files
####################
source("functions/Sampling.R")

####################
# Start of Code
####################

################################################################################
# Example 1
####################

# Beta distribution PDF with default values alpha=2 & beta=5
BetaPDF <- function(x, alpha=2, beta=5){
  # Recall Beta distribution is defined only on [0,1]
  
  numerator <- x^(alpha-1) * (1-x)^(beta-1) * gamma(alpha + beta)
  denominator <- gamma(alpha) * gamma(beta)
  
  return(numerator/denominator)
}

n <- 500

## vector of estimates for c
c_estimates <- numeric(n+1)
c_estimates[1] <- 1
correct.classification <- rep(TRUE, n)

X <- runif(n)
U <- runif(n)
  
## algorithm implementation
for(i in 1:n){
  ratio <- BetaPDF(X[i])/GetUniformPDF(X[i], 0, 1)
  
  if(c_estimates[i]*U[i] < BetaPDF(X[i]) && 2.4576*U[i] > BetaPDF(X[i])){
    correct.classification[i] <- FALSE
  }
  
  c_estimates[i+1] <- max(c_estimates[i], ratio)
}

c <- 2.4576
x.axis <- seq(0,1, length=200)
plot(x.axis, dbeta(x.axis, 2, 5), type="l", xlab="x", ylab="Density", main="Empirical Supremum Rejection - Mistakes")
lines(x.axis, c*dunif(x.axis, 0, 1), col="blue")
points(X[correct.classification], c*U[correct.classification], pch=4, col="cadetblue1")
points(X[!correct.classification], c*U[!correct.classification], pch=19, col="darkviolet")
legend("bottomleft", 
       legend=c("Beta(2,5)",
                "Envelope Function",
                "Correctly Classified",
                "Incorrectly Accepted"),
       lty=c(1, 1, 0, 0),
       pch=c(26, 26, 4, 19),
       col=c("black", "blue", "cadetblue1", "darkviolet")
)


## plot the estimates
plot(1:200, c_estimates[1:200], xlab="Iterations", ylab="c estimate", type="o",
     main="Estimating c - Convergence")
## plot the true value of c (can be found analytically)
abline(h=2.4576, col="red")

plot(1:(n+1), log10(2.4576 - c_estimates), ylab="log(error) in base 10", type="l", xlab="Iteration", main="Estimating c - Logarithmic Error")

################################################################################
# Example 2
####################

n <- 500

## vector of estimates for c
c_estimates <- numeric(n)
c_estimates[1] <- 1

## algorithm implementation
for(i in 2:n){
  U <- SimulateU()
  X <- rt(1, 2)
  ratio <- dnorm(X, 0, 1)/dt(X, 2)
  
  ## We don't care if the sample is accepted or rejected for estimating c
  
  c_estimates[i] <- max(c_estimates[i-1], ratio)
}

## plot the estimates
plot(1:n, c_estimates, xlab="Iterations", ylab="c estimate", 
     main="Empirical Supremum Rejection Sampling
     Normal(0,1) sampling from student-t2")
## plot the true value of c (can be found analytically)
#abline(h=0, col="red")



################################################################################

f <- function(X){
  return(dnorm(X, 0, 1))
}
g <- function(X){
  return(dt(X, 2))
}
candidate.distribution <- function(){
  return(rt(1, 2))
}

c_estimates <- FindC(f, g, candidate.distribution, TRUE)
plot(1:length(c_estimates), c_estimates, xlab="Iterations", ylab="c estimate", 
     main="Empirical Supremum Rejection Sampling
     Normal(0,1) sampling from student-t2")

hist(GenerateSample(1000, f, g, candidate.distribution))
