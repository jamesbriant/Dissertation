library(mvtnorm) # for multivariate normal
library(plotly) # for 3d plots

GenerateKalmanData <- function(f, g, X0, W, V, n=20){
  # Generate the unobserved system X(t) and the observed data Y(t)
  # X(t+1) = fX(t) + W
  # Y(t) = gX(t) + V
  #
  # INPUTS:
  #   f:    matrix
  #   g:    matrix
  #   x0:   start value
  #   W:    variance matrix
  #   V:    variance matrix
  #   n:    number of data points to be produced
  # OUTPUTS:
  #   X:    list of system states
  #   Y:    list of observed states
  
  f.matrix <- as.matrix(f)
  g.matrix <- as.matrix(g)
  X0.matrix <- as.matrix(X0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  # Generate system data X and observations Y
  X <- list()
  X[[1]] <- t(rmvnorm(1, f.matrix %*% X0.matrix, W.matrix))
  Y <- list()
  Y[[1]] <- t(rmvnorm(1, g.matrix %*% X[[1]], V.matrix))
  for (i in 2:n){
    X[[i]] <- t(rmvnorm(1, f.matrix %*% X[[i-1]], W.matrix))
    Y[[i]] <- t(rmvnorm(1, g.matrix %*% X[[i]], V.matrix))
  }
  
  return(list(X=X, Y=Y))
}



ApplyKalmanFilter <- function(Y, f, g, m0, C0, W, V){
  # INPUTS:
  #   Y:    list of observed states
  #   f:    matrix
  #   g:    matrix
  #   m0:   initial estimate for mean of distribution of X
  #   C0:   initial estimate for variance of distribution of X
  #   W:    variance matrix
  #   V:    variance matrix
  # OUTPUTS:
  #   a:    list of priors for the distribution mean of X given y at t-1
  #   m:    list of posteriors for the distribution mean of X given y at t
  #   C:    list of posterior variance matrices for distn  of X at time t
  #   R:    list of prior variance matrices for distn of X at time t+1
  #   Q:    list of prior variance matrices for distn of Y at time t+1
  #   K:    list of Kalman gain values for each time 
  
  CalcKalmanGain <- function(R, G, Q){
    if(dim(Q)[1]==1){
      return(1/as.numeric(Q) * R %*% t(G))
    }
    else{
      return(R %*% t(G) %*% solve(Q))
    }
  }
  
  
  f.matrix <- as.matrix(f)
  g.matrix <- as.matrix(g)
  m0.matrix <- as.matrix(m0)
  C0.matrix <- as.matrix(C0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  n <- length(Y)
  p <- nrow(f.matrix)
  
  a <- list()
  m <- list()
  C <- list()
  R <- list()
  Q <- list()
  K <- list()
  
  # Perform initial calculation using initial conditions
  a[[1]] <- f.matrix %*% m0.matrix
  R[[1]] <- f.matrix %*% C0.matrix %*% t(f.matrix) + W.matrix
  Q[[1]] <- g.matrix %*% R[[1]] %*% t(g.matrix) + V.matrix
  
  K[[1]] <- CalcKalmanGain(R[[1]], g.matrix, Q[[1]])
  m[[1]] <- a[[1]] + K[[1]] %*% (Y[[1]] - g.matrix %*% a[[1]])
  C[[1]] <- (diag(p) -  K[[1]] %*% g.matrix) %*% R[[1]]
  
  for(t in 2:n){
    # predict
    a[[t]] <- f.matrix %*% m[[t-1]]
    R[[t]] <- f.matrix %*% C[[t-1]] %*% t(f.matrix) + W.matrix
    Q[[t]] <- g.matrix %*% R[[t]] %*% t(g.matrix) + V.matrix
    
    # update
    K[[t]] <- CalcKalmanGain(R[[t]], g.matrix, Q[[t]])
    m[[t]] <- a[[t]] + K[[t]] %*% (Y[[t]] - g.matrix %*% a[[t]])
    C[[t]] <- (diag(p) - K[[t]] %*% g.matrix) %*% R[[t]]
  }
  
  return(list(a=a, m=m, C=C, R=R, Q=Q, K=K))
}




#########################################################################
# Example: X-1D; Y-1D;

f <- 2
g <- 1

X0 <- 1

# Set variances
W <- 4
V <- 8

n <- 7
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y



# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyKalmanFilter(Y, f, g, m0, C0, W, V)



plot(1:n, unlist(X), type="l", xlab="time", ylab="state")
lines(1:n, unlist(kalmanSolution$m), col="red")
points(1:n, unlist(kalmanSolution$m), col="red")
#points(1:n, unlist(kalmanSolution$a), col="blue")
legend("topleft", 
       legend=c("X(t) - Unobserved System",
                "m(t) - Posterior",
                "a(t) - Prior"), 
       col=c("black", "red", "blue"), 
       lty=c(1, 1, 0), 
       pch=c(26, 1, 1)
       )

confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * unlist(kalmanSolution$C)
segments(1:n, unlist(kalmanSolution$m) - confidenceRange, 1:n, unlist(kalmanSolution$m) + confidenceRange, col="red")







###########################################################################
# Example: X-2D; Y-1D; 

f <- matrix(c(2, 0, 0, 2), nrow=2, byrow=TRUE)
g <- matrix(c(1, 1), nrow=1)

W <- matrix(c(1, 0, 0, 2), nrow=2, byrow=TRUE)
V <- 4

#initial conditions
X0 <- matrix(c(1, 1), nrow=2)

n <- 6
systemData <- GenerateKalmanData(f, g, X0, W, V, n)
X <- systemData$X
Y <- systemData$Y

# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyKalmanFilter(Y, f, g, m0, C0, W, V)

data <- data.frame(time=1:n,
                   X1=unlist(X)[2*(1:n)-1], 
                   X2=unlist(X)[2*(1:n)], 
                   Y=unlist(Y),
                   a1=unlist(kalmanSolution$a)[2*(1:n)-1], 
                   a2=unlist(kalmanSolution$a)[2*(1:n)],
                   m1=unlist(kalmanSolution$m)[2*(1:n)-1],
                   m2=unlist(kalmanSolution$m)[2*(1:n)]
                   )

plot(data$X1, data$X2, type="l", col="blue", xlab="X1", ylab="X2", main="Predictions and Estimates for 2D X using 1D Y")
lines(data$m1, data$m2, type="l", col="orange")
lines(data$a1, data$a2, col="green")
lines(2:100, 2:100, lty=2, col="grey")
legend("topleft", 
       legend=c("X(t) - Unobserved System",
                "m(t) - Estimated state (posterior)",
                "a(t) - Predicted state (prior)",
                "X1=X2"),
       col=c("blue", "orange", "green", "grey"),
       lty=c(1, 1, 1, 2))

fig <- plot_ly(data, x=~X1, y=~X2, z=~1:n, type = 'scatter3d', mode = 'lines', name="X(t) - Unobserved System")
fig <- fig %>% add_trace(data, x=~m1, y=~m2, z=~1:n, type = 'scatter3d', mode = 'lines', name="m(t) - Estimated state (posterior)")
fig <- fig %>% add_trace(data, x=~a1, y=~a2, z=~time, type = 'scatter3d', mode = 'lines', name="a(t) - Predicted state (prior)")
fig <- fig %>% layout(scene = list(zaxis = list(title="time")))
fig




###########################################################################
# Example: X-2D; Y-2D; 

f <- matrix(c(2, -1, 0, 2), nrow=2, byrow=TRUE)
g <- matrix(c(1, 0, 0, 1), nrow=2, byrow=TRUE)

W <- matrix(c(16, 0, 0, 16), nrow=2, byrow=TRUE)
V <- matrix(c(49, 0, 0, 49), nrow=2, byrow=TRUE)

#initial conditions
X0 <- matrix(c(1, 1), nrow=2)

n <- 6
systemData <- GenerateKalmanData(f, g, X0, W, V, n)
X <- systemData$X
Y <- systemData$Y

# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyKalmanFilter(Y, f, g, m0, C0, W, V)

data <- data.frame(time=1:n,
                   X1=unlist(X)[2*(1:n)-1], 
                   X2=unlist(X)[2*(1:n)], 
                   Y1=unlist(Y)[2*(1:n)-1],
                   Y2=unlist(Y)[2*(1:n)],
                   a1=unlist(kalmanSolution$a)[2*(1:n)-1], 
                   a2=unlist(kalmanSolution$a)[2*(1:n)],
                   m1=unlist(kalmanSolution$m)[2*(1:n)-1],
                   m2=unlist(kalmanSolution$m)[2*(1:n)]
)

fig <- plot_ly(data, x=~X1, y=~X2, z=~1:n, type = 'scatter3d', mode = 'lines', name="X(t) - Unobserved System")
fig <- fig %>% add_trace(data, x=~m1, y=~m2, z=~1:n, type = 'scatter3d', mode = 'lines', name="m(t) - Estimated state (posterior)")
fig <- fig %>% add_trace(data, x=~a1, y=~a2, z=~time, type = 'scatter3d', mode = 'lines', name="a(t) - Predicted state (prior)")
fig <- fig %>% layout(scene = list(zaxis = list(title="time")))
fig

















