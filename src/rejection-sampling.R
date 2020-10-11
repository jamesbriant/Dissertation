#########################################################################
# Sample f=Beta(2,5) from g=Unif(0,1) using Rejection Sampling

drawFromg <- function(n){
  # g ~ Uniform[0,2.5]
  return(runif(n, 0 , 2.5))  
  # f(x) < 2.5 for x on [0,1]. Hence this is NOT most computationally efficient!
}

f <- function(x){
  # Beta distribution with alpha=2 & beta=5
  # Recall Beta distribution is only define on [0,1]
  
  alpha <- 2
  beta <- 5
  
  numerator <- x^(alpha-1) * (1-x)^(beta-1) * gamma(alpha + beta)
  denominator <- gamma(alpha) * gamma(beta)
  
  return(numerator/denominator)
}

drawSample <- function(n){
  x <- runif(n, 0, 1)
  
  # This is the g function!
  y <- drawFromg(n)
  
  fx <- f(x)
  
  return(x[y < fx])
}


data <- drawSample(10000)
hist(data, freq=FALSE, xlab="x", main="Beta(2,5) Distribution from 10,000 Draws")
x <- seq(from=0, to=1, length.out=250)
lines(x, dbeta(x, 2, 5), col="red")
legend(0.5, 2.4, legend=c("Rejection Sampling", "True distribution"), pch="-", col=c("black", "red"))











