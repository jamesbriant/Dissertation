
####################
# Source files
####################
source("functions/Kalman.R")
source("functions/Plotting.R")



#########################################################################
# Example: X-1D; Y-1D;

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



# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyExtendedKalmanFilter(Y, f, g, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(kalmanSolution$m), 
                   unlist(kalmanSolution$C)
)



#########################################################################
# Example: X-1D; Y-1D; Logisitic equation - population growth

r <- 0.1
K <- 6000

f <- function(x){
  X <- x[1]
  t <- x[2]
  
  dPdt <- r*X*(1-X/K)
  
  return(X + dPdt)
}

g <- function(x){
  X <- x[1]
  t <- x[2]
  return(X*0.9)
}

X0 <- 500

# Set variances
W <- 16000
V <- 400

n <- 100
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y



# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyExtendedKalmanFilter(Y, f, g, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(kalmanSolution$m), 
                   unlist(kalmanSolution$C)
)

################################################################################
# Example: X-1D; Y-1D;
# KALMAN vs EXTENDED
# stationary model

# This obviously doesn't make sense. The linear approximation of a DLM is 
# obviously linear, hence the approximation is exact.

A <- as.matrix(1)
B <- as.matrix(1)

f <- function(x){
  return(A %*% x[1])
}
g <- function(x){
  return(B %*% x[1])
}

X0 <- as.matrix(0)

# Set variances
W <- as.matrix(1)
V <- as.matrix(1)

set.seed(2022)

n <- 15
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y
plot(0:n, c(0, unlist(Y)), xlab="time", ylab="Obersations", main="Obersavtions Time Series")


# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

#SavePlotToPNG("plots/KalmanStationary.png")
kalmanSolution <- ApplyKalmanFilter(Y, A, B, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(kalmanSolution$m), 
                   unlist(kalmanSolution$C)
)


extended.kalman.solution <- ApplyExtendedKalmanFilter(Y, f, g, m0, C0, W, V)

confidence=0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * unlist(extended.kalman.solution$C)
polygon(c(1:n, n:1), 
        c(unlist(extended.kalman.solution$m) - confidenceRange, rev(unlist(extended.kalman.solution$m) + confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:n, unlist(extended.kalman.solution$m) - confidenceRange, col="blue", lty=2)
lines(1:n, unlist(extended.kalman.solution$m) + confidenceRange, col="blue", lty=2)
lines(1:n, unlist(extended.kalman.solution$m), col="blue", lty=1)
