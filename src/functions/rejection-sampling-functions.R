drawFromg <- function(n){
  # g ~ Uniform[0,1]
  return(runif(n, 0 , 1))  
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
  y <- runif(n, 0, 3)
  # f(x) < 2.5 for x on [0,1]. Hence this is NOT most computationally efficient!
  
  # This is the g function!
  x <- drawFromg(n)
  
  fx <- f(x)
  
  return(x[y < fx])
}

