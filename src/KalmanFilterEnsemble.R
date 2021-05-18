
####################
# Source files
####################
source("functions/Kalman.R")
source("functions/Plotting.R")


#########################################################################
# Example: X-1D; Y-1D;
# cosine example

f <- function(x){
  X <- x[1]
  t <- x[2]
  return(3*cos(pi*t/10))
}

g <- function(x){
  X <- x[1]
  t <- x[2]
  return(X)
}

X0 <- 3

# Set variances
W <- 0.8
V <- 1.2

n <- 30
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y


M <- 50

# Initial estimates
m0 <- t(as.matrix(X0 + rnorm(M, 0, W))) # mean
C0 <- W # variance

ensemble.kalman.solution <- ApplyEnsembleKalmanFilter(Y, f, g, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(ensemble.kalman.solution$mbar), 
                   unlist(ensemble.kalman.solution$C)
)


#########################################################################
# Example: X-1D; Y-1D;
# Population growth example

r <- 0.1
K <- 6000

# f <- function(x){
#   X <- x[1]
#   t <- x[2]
#   
#   dPdt <- r*X*(1-X/K)
#   
#   return(X + dPdt)
# }
# 
# g <- function(x){
#   X <- x[1]
#   t <- x[2]
#   return(X*0.9)
# }

f <- function(x){
  X <- x[1]
  t <- x[2]
  
  return(X)
}

g <- function(x){
  X <- x[1]
  t <- x[2]
  return(X)
}

X0 <- 500

# Set variances
W <- 600
V <- 40

n <- 100
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y



# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

extended.kalman.solution <- ApplyExtendedKalmanFilter(Y, f, g, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(extended.kalman.solution$m), 
                   unlist(extended.kalman.solution$C)
)




M <- 50

# Initial estimates
m0 <- t(as.matrix(X0 + rnorm(M, 0, sqrt(W)))) # mean
C0 <- W # variance

ensemble.kalman.solution <- ApplyEnsembleKalmanFilter(Y, f, g, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(ensemble.kalman.solution$mbar), 
                   unlist(ensemble.kalman.solution$C)
)

confidence=0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * unlist(ensemble.kalman.solution$C)
polygon(c(1:n, n:1), 
        c(unlist(ensemble.kalman.solution$mbar) - confidenceRange, rev(unlist(ensemble.kalman.solution$mbar) + confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:n, unlist(ensemble.kalman.solution$mbar) - confidenceRange, col="blue", lty=2)
lines(1:n, unlist(ensemble.kalman.solution$mbar) + confidenceRange, col="blue", lty=2)
lines(1:n, unlist(ensemble.kalman.solution$mbar), col="blue", lty=1)
