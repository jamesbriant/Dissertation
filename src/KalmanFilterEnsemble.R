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

ApplyEnsembleKalmanFilter <- function(Y, f, g, m0, C0, W, V){
  # INPUTS:
  #   Y:    list of observed states
  #   f:    function
  #   g:    function
  #   m0:   p by M matrix of initial estimate for mean of distributions of Xi's
  #   C0:   initial estimate for variance of distributions of Xi's
  #   W:    variance matrix
  #   V:    variance matrix
  # OUTPUTS:
  #   a:    list of priors for the distribution mean of Xbar given y at t-1
  #   m:    list of posteriors for the distribution mean of Xbar given y at t
  #   C:    list of posterior variance matrices for distn  of Xbar at time t
  #   R:    list of prior variance matrices for distn of Xbar at time t+1
  #   Q:    list of prior variance matrices for distn of Ybar at time t+1
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
  
  p <- dim(m0)[1] # dimension of X
  M <- dim(m0)[2] # size of ensemble

  m0.matrix <- as.matrix(m0)
  C0.matrix <- as.matrix(C0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  n <- length(Y) # number of time steps
  q <- dim(Y[[1]])[1] # dimension of Y
    
  a <- list()     # list of matrices, each matrix is M by p
  abar <- list()  # list of vectors
  m <- list()     # list of matrices, each matrix is M by p
  mbar <- list()  # list of vectors
  C <- list()
  R <- list()
  #Q <- list()
  Qprime <- list()
  K <- list()
  
  f.m0 <- matrix(0, nrow=p, ncol=M)
  g.fm0 <- matrix(0, nrow=q, ncol=M)
  for(i in 1:M){
    f.m0[, i] <- as.matrix(f(c(m0.matrix[, i], 1))) + rmvnorm(1, zeros(p, 1), W.matrix) 
    g.fm0[, i] <- as.matrix(g(c(f.m0[, i], 1)))
  }
  
  # Perform initial calculation using initial conditions
  a[[1]] <- f.m0
  abar[[1]] <- t(as.matrix(apply(a[[1]], 1, mean)))
  Ea <- a[[1]] - matrix(rep(abar[[1]], M), byrow=FALSE, ncol=M)
  R[[1]] <- 1/(M-1) * Ea %*% t(Ea)
  #Q[[1]] <- dg.fm0 %*% R[[1]] %*% t(dg.fm0) + V.matrix
  
  y <- matrix(0, nrow=q, ncol=M)
  for(i in 1:M){
    y[, i] <- rmvnorm(q, Y[[1]], V.matrix) 
  }
  ybar <- t(as.matrix(apply(y, 1, mean)))
  Ey <- y - matrix(rep(ybar, M), byrow=FALSE, ncol=M)
  Qprime[[1]] <- 1/(M-1) * Ey %*% t(Ey)
  
  K[[1]] <- Ea %*% t(Ey) %*% solve(Ey %*% t(Ey))
  
  m[[1]] <- matrix(0, nrow=p, ncol=M)
  for(i in 1:M){
    m[[1]][, i] <- a[[1]][, i] + K[[1]] %*% (y[, i] - g.fm0[, i])
  }
  mbar[[1]] <- t(as.matrix(apply(m[[1]], 1, mean)))
  Em <- m[[1]] - matrix(rep(mbar[[1]], M), byrow=FALSE, ncol=M)
  
  C[[1]] <- 1/(M-1) * Em %*% t(Em)
  
  for(t in 2:n){
    f.m <- matrix(0, nrow=p, ncol=M)
    g.fm <- matrix(0, nrow=q, ncol=M)
    for(i in 1:M){
      f.m[, i] <- as.matrix(f(c(m[[t-1]][, i], t))) + rmvnorm(1, zeros(p, 1), W.matrix) 
      g.fm[, i] <- g(c(f.m[, i], t))
    }
    
    # predict
    a[[t]] <- f.m
    abar[[t]] <- t(as.matrix(apply(a[[t]], 1, mean)))
    Ea <- a[[t]] - matrix(rep(abar[[t]], M), byrow=FALSE, ncol=M)
    R[[t]] <- 1/(M-1) * Ea %*% t(Ea)
    
    # update
    y <- matrix(0, nrow=q, ncol=M)
    for(i in 1:M){
      y[, i] <- rmvnorm(q, Y[[t]], V.matrix) 
    }
    ybar <- t(as.matrix(apply(y, 1, mean)))
    Ey <- y - matrix(rep(ybar, M), byrow=FALSE, ncol=M)
    Qprime[[t]] <- 1/(M-1) * Ey %*% t(Ey)
    
    K[[t]] <- Ea %*% t(Ey) %*% solve(Ey %*% t(Ey))
    
    m[[t]] <- matrix(0, nrow=p, ncol=M)
    for(i in 1:M){
      m[[t]][, i] <- a[[t]][, i] + K[[t]] %*% (y[, i] - g.fm[, i])
    }
    mbar[[t]] <- t(as.matrix(apply(m[[t]], 1, mean)))
    Em <- m[[t]] - matrix(rep(mbar[[t]], M), byrow=FALSE, ncol=M)
    
    C[[t]] <- 1/(M-1) * Em %*% t(Em)
  }
  
  return(list(abar=abar, mbar=mbar, K=K, C=C, R=R, Qprime=Qprime, K=K))
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


M <- 100

# Initial estimates
m0 <- t(as.matrix(X0 + rnorm(M, 0, W))) # mean
C0 <- W # variance

kalmanSolution <- ApplyEnsembleKalmanFilter(Y, f, g, m0, C0, W, V)


plot(1:n, unlist(X), type="l", xlab="time", ylab="state")
lines(1:n, unlist(kalmanSolution$mbar), col="red")
points(1:n, unlist(kalmanSolution$mbar), col="red")
confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * unlist(kalmanSolution$C)
segments(1:n, unlist(kalmanSolution$m) - confidenceRange, 1:n, unlist(kalmanSolution$m) + confidenceRange, col="red")

legend(4, 4.9,
       legend=c("X(t) - Unobserved System",
                "m(t) - Estimated State",
                "95% confident interval"), 
       col=c("black", "red", "red"), 
       lty=c(1, 1, NA), 
       pch=c(NA, "o", "|")
)


#########################################################################
# Example: X-2D; Y-1D;

f <- function(x){
  X1 <- x[1]
  X2 <- x[2]
  t <- x[3]
  return(c(3*cos(pi*t/10), 6*sin(pi*t/7)))
}

g <- function(x){
  X1 <- x[1]
  X2 <- x[2]
  t <- x[3]
  return(X1 + X2)
}

#initial conditions
X0 <- matrix(c(3, 0), nrow=2)

# Set variances
W <- matrix(c(1, 0, 0, 1), nrow=2, byrow=TRUE)
V <- 1.3

n <- 30
data <- GenerateExtendedKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y

# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyExtendedKalmanFilter(Y, f, g, m0, C0, W, V)

data <- data.frame(time=1:n,
                   X1=unlist(X)[2*(1:n)-1], 
                   X2=unlist(X)[2*(1:n)], 
                   Y=unlist(Y),
                   a1=unlist(kalmanSolution$a)[2*(1:n)-1], 
                   a2=unlist(kalmanSolution$a)[2*(1:n)],
                   m1=unlist(kalmanSolution$m)[2*(1:n)-1],
                   m2=unlist(kalmanSolution$m)[2*(1:n)]
)

fig <- plot_ly(data, x=~X1, y=~X2, z=~1:n, type = 'scatter3d', mode = 'lines')
fig <- fig %>% add_trace(data, x=~m1, y=~m2, z=~1:n, type = 'scatter3d', mode = 'lines')
fig <- fig %>% add_trace(data, x=~a1, y=~a2, z=~time, type = 'scatter3d', mode = 'lines')
fig


#########################################################################
# Example: X-2D; Y-2D;

f <- function(x){
  X1 <- x[1]
  X2 <- x[2]
  t <- x[3]
  return(c(3*cos(pi*t/10), 6*sin(pi*t/7)))
}

g <- function(x){
  X1 <- x[1]
  X2 <- x[2]
  t <- x[3]
  return(c(X1, X2))
}

#initial conditions
X0 <- matrix(c(3, 0), nrow=2)

# Set variances
W <- matrix(c(1, 0, 0, 1), nrow=2, byrow=TRUE)
V <- matrix(c(2, 0, 0, 2), nrow=2, byrow=TRUE)

n <- 30
data <- GenerateExtendedKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y

# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalmanSolution <- ApplyExtendedKalmanFilter(Y, f, g, m0, C0, W, V)

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

fig <- plot_ly(data, x=~X1, y=~X2, z=~1:n, type = 'scatter3d', mode = 'lines')
fig <- fig %>% add_trace(data, x=~m1, y=~m2, z=~1:n, type = 'scatter3d', mode = 'lines')
fig <- fig %>% add_trace(data, x=~a1, y=~a2, z=~time, type = 'scatter3d', mode = 'lines')
fig
