library(mvtnorm)  # for multivariate normal - rmvnorm()
library(pracma)   # for jacobian()
library(plotly)

GenerateExtendedKalmanData <- function(f, g, X0, W, V, n=20){
  # Generate the unobserved system X(t) and the observed data Y(t)
  # X(t+1) = fX(t) + W
  # Y(t) = gX(t) + V
  #
  # INPUTS:
  #   f:    function
  #   g:    function
  #   x0:   start value
  #   W:    variance matrix
  #   V:    variance matrix
  #   n:    number of data points to be produced
  # OUTPUTS:
  #   X:    list of system states
  #   Y:    list of observed states
  
  X0.matrix <- as.matrix(X0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  # Generate system data X and observations Y
  X <- list()
  X[[1]] <- t(rmvnorm(1, f(c(X0.matrix, 1)), W.matrix))
  Y <- list()
  Y[[1]] <- t(rmvnorm(1, g(c(X[[1]], 1)), V.matrix))
  for (i in 2:n){
    X[[i]] <- t(rmvnorm(1, f(c(X[[i-1]], i)), W.matrix))
    Y[[i]] <- t(rmvnorm(1, g(c(X[[i]], i)), V.matrix))
  }
  
  return(list(X=X, Y=Y))
}

ApplyExtendedKalmanFilter <- function(Y, f, g, m0, C0, W, V){
  # INPUTS:
  #   Y:    list of observed states
  #   f:    function
  #   g:    function
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
  #   K:    list of Kalman gain values for each time t
  
  CalcKalmanGain <- function(R, G, Q){
    if(dim(Q)[1]==1){
      return(1/as.numeric(Q) * R %*% t(G))
    }
    else{
      return(R %*% t(G) %*% solve(Q))
    }
  }
  
  RemoveTime <- function(v){
    # assume the last variable is t
    # v must be a 1 by n matrix
    len <- dim(as.matrix(v))[2]
    return(as.matrix(v[, 1:(len-1)]))
  }
  
  m0.matrix <- as.matrix(m0)
  C0.matrix <- as.matrix(C0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  f.m0 <- as.matrix(f(c(m0.matrix, 1)))
  g.fm0 <- as.matrix(g(c(f.m0, 1)))
  df.m0 <- RemoveTime(jacobian(f, c(m0.matrix, 1)))
  dg.fm0 <- t(RemoveTime(jacobian(g, c(f.m0, 1))))
  
  n <- length(Y)
  p <- dim(dg.fm0)[2]
  
  a <- list()
  m <- list()
  C <- list()
  R <- list()
  Q <- list()
  K <- list()
  
  # Perform initial calculation using initial conditions
  a[[1]] <- f.m0
  R[[1]] <- df.m0 %*% C0.matrix %*% t(df.m0) + W.matrix
  Q[[1]] <- dg.fm0 %*% R[[1]] %*% t(dg.fm0) + V.matrix
  
  K[[1]] <- CalcKalmanGain(R[[1]], dg.fm0, Q[[1]])
  m[[1]] <- a[[1]] + K[[1]] %*% (Y[[1]] - g.fm0)
  C[[1]] <- (diag(p) -  K[[1]] %*% dg.fm0) %*% R[[1]]
  
  for(t in 2:n){
    f.m <- as.matrix(f(c(m[[t-1]], t)))
    g.fm <- as.matrix(g(c(f.m, t)))
    df.m <- RemoveTime(jacobian(f, c(m[[t-1]], t)))
    dg.fm <- t(RemoveTime(jacobian(g, c(f.m, t))))
    
    # predict
    a[[t]] <- f.m
    R[[t]] <- df.m %*% C[[t-1]] %*% t(df.m) + W.matrix
    Q[[t]] <- dg.fm %*% R[[t]] %*% t(dg.fm) + V.matrix
    
    # update
    K[[t]] <- CalcKalmanGain(R[[t]], dg.fm, Q[[t]])
    m[[t]] <- a[[t]] + K[[t]] %*% (Y[[t]] - g.fm)
    C[[t]] <- (diag(p) - K[[t]] %*% dg.fm) %*% R[[t]]
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
  
  RemoveTime <- function(v){
    # assume the last variable is t
    # v must be a 1 by n matrix
    len <- dim(as.matrix(v))[2]
    return(as.matrix(v[, 1:(len-1)]))
  }
  
  n <- length(a)
  
  s <- list()
  S <- list()
  
  # Perform initial calculation using initial conditions
  s[[1]] <- m[[n]]
  S[[1]] <- C[[n]]
  
  for(i in 2:n){
    t <- n - i + 1
    df.m <- RemoveTime(jacobian(f, c(m[[t+1]], t)))
    A <- C[[t]] %*% t(df.m) %*% Invert(R[[t+1]])
    
    # predict
    s[[i]] <- m[[t]] + A %*% (as.matrix(s[[i-1]] - a[[t+1]]))
    S[[i]] <- C[[t]] + A %*% (as.matrix(S[[i-1]] - R[[t+1]])) %*% t(A)
  }
  
  return(list(s=rev(s), S=rev(S)))
}



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
data <- GenerateExtendedKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y



# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyExtendedKalmanFilter(Y, f, g, m0, C0, W, V)
kalmanSmoothedSolution <- ApplyKalmanSmoother(f, kalmanSolution$a, kalmanSolution$m, kalmanSolution$R, kalmanSolution$C)

plot(1:n, unlist(X), type="l", xlab="time", ylab="state")
lines(1:n, unlist(kalmanSmoothedSolution$s), col="red")
confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * unlist(kalmanSmoothedSolution$S)
segments(1:n, unlist(kalmanSmoothedSolution$s) - confidenceRange, 1:n, unlist(kalmanSmoothedSolution$s) + confidenceRange, col="red")
legend(5, 4.5,
       legend=c("X(t) - Unobserved System",
                "m(t) - Smoothed Estimate",
                "95% confident interval"), 
       col=c("black", "red", "red"), 
       lty=c(1, 1, NA), 
       pch=c(NA, NA, "|")
)


















