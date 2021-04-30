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
c_estimates <- numeric(n)
c_estimates[1] <- 1
  
## algorithm implementation
for(i in 2:n){
  U <- SimulateU()
  X <- GetUniformSample(1, 0, 1)
  ratio <- BetaPDF(X)/GetUniformPDF(X, 0, 1)
  
  ## We don't care if the sample is accepted or rejected for estimating c
  
  c_estimates[i] <- max(c_estimates[i-1], ratio)
}

## plot the estimates
plot(1:n, c_estimates, xlab="Iterations", ylab="c estimate", type="o",
     main="Empirical Supremum Rejection Sampling
     Beta(2,5) sampling from Unif(0,1)")
## plot the true value of c (can be found analytically)
abline(h=2.4576, col="red")

plot(1:n, log10(2.4576 - c_estimates), ylab="log(error) - base 10", type="l", xlab="Iteration", main="Empirical Supremum Rejection Sampling - Estimating c")

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
