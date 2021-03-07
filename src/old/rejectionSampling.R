#########################################################################
# Sample f=Beta(2,5) from g=Unif(0,1) using Rejection Sampling

##########################
# OLD VERSION

source('functions/rejection-sampling-functions.R')

data <- drawSample(10000)
hist(data, freq=FALSE, xlab="x", main="Beta(2,5) Distribution from 2500 Draws")
x <- seq(from=0, to=1, length.out=250)
lines(x, dbeta(x, 2, 5), col="red")
legend(0.5, 2.4, legend=c("Rejection Sampling", "True distribution"), pch="-", col=c("black", "red"))
rug(data)







###################################
# NEW VERSION


g <- function(X){
  # g ~ Uniform(0,1)
  return(1)
}

DrawFromg <- function(n){
  return(runif(n, 0, 1))
}

f <- function(X){
  # Beta distribution with alpha=2 & beta=5
  # Recall Beta distribution is defined only on [0,1]
  
  alpha <- 2
  beta <- 5
  
  numerator <- X^(alpha-1) * (1-X)^(beta-1) * gamma(alpha + beta)
  denominator <- gamma(alpha) * gamma(beta)
  
  return(numerator/denominator)
}

n <- 1000

U <- runif(n, 0, 1)
X <- DrawFromg(n)
Y <- X[U < f(X)/2.5]
hist(Y, freq=FALSE, xlab="x", main="Beta(2,5) Distribution from 1000 Draws of Unif(0,1)", xlim=c(0,1))
lines(seq(0, 1, length=500), f(seq(0, 1, length=500)), col="red")
legend(0.6, 2.4, legend=c("Rejection Sampling", "Beta(2,5)"), pch="-", col=c("black", "red"))
rug(Y)





















