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

ApplyKalmanSmoother <- function(f, a, m, R, C){
  # INPUTS:
  #   f:    matrix
  #   a:    list of prior means of distribution of X
  #   m:    list of posterior means of distribution of X
  #   R:    list of prior variances of distribution of X
  #   C:    list of posterior variances of distribution of X
  # OUTPUTS:
  #   s:    list of smoothed distribution mean of X given y at t
  #   S:    list of smoothed distribution variance matrices of X at time t
  
  Invert <- function(R){
    if(dim(R)[1]==1){
      return(as.matrix(1/as.numeric(R)))
    }
    else{
      return(solve(R))
    }
  }
  
  f.matrix <- as.matrix(f)
  
  n <- length(a)
  
  s <- list()
  S <- list()
  
  # Perform initial calculation using initial conditions
  s[[1]] <- m[[n]]
  S[[1]] <- C[[n]]
  
  for(i in 2:n){
    t <- n - i + 1
    A <- C[[t]] %*% t(f) %*% Invert(R[[t+1]])
    
    # predict
    s[[i]] <- m[[t]] + A %*% (as.matrix(s[[i-1]] - a[[t+1]]))
    S[[i]] <- C[[t]] + A %*% (as.matrix(S[[i-1]] - R[[t+1]])) %*% t(A)
  }
  
  return(list(s=rev(s), S=rev(S)))
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
kalmanSmoothedSolution <- ApplyKalmanSmoother(f, kalmanSolution$a, kalmanSolution$m, kalmanSolution$R, kalmanSolution$C)



plot(1:n, unlist(X), type="l", xlab="time", ylab="state")
lines(1:n, unlist(kalmanSmoothedSolution$s), col="red")
confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * unlist(kalmanSmoothedSolution$S)
segments(1:n, unlist(kalmanSmoothedSolution$s) - confidenceRange, 1:n, unlist(kalmanSmoothedSolution$s) + confidenceRange, col="red")
legend("topleft",
       legend=c("X(t) - Unobserved System",
                "m(t) - Smoothed Estimate",
                "95% confident interval"), 
       col=c("black", "red", "red"), 
       lty=c(1, 1, NA), 
       pch=c(NA, NA, "|")
)



















