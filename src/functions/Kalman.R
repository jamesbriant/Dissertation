library(mvtnorm)

#############
# Defintions of all Kalman functions
#############

GenerateKalmanData <- function(f, g, X0, W, V, n=20){
  # Generate the unobserved system X(t) and the observed data Y(t)
  # X(t+1) = fX(t) + W
  # Y(t) = gX(t) + V
  #
  # INPUTS:
  #   f:    function, must return a vector
  #   g:    function, must return a vector
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
  for (t in 2:n){
    X[[t]] <- t(rmvnorm(1, f(c(X[[t-1]], t)), W.matrix))
    Y[[t]] <- t(rmvnorm(1, g(c(X[[t]], t)), V.matrix))
  }
  
  return(list(X=X, Y=Y))
}

LaplacianPDF <- function(x, mu=0, b=1){
  return(exp(-abs(x-mu)/b)/(2*b))
}

rLaplacian <- function(b=1){
  #mu=0
  return(rexp(1, 1/b) - rexp(1, 1/b))
}

GenerateKalmanDataLaplacian <- function(f, g, X0, W, b, n=20){
  # This only works in 1D
  # Generate the unobserved system X(t) and the observed data Y(t)
  # X(t+1) = fX(t) + W
  # Y(t) = gX(t) + V
  #
  # INPUTS:
  #   f:    function, must return a vector
  #   g:    function, must return a vector
  #   x0:   start value
  #   W:    variance matrix
  #   b:    laplacian scale parameter
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
  Y[[1]] <- g(c(X[[1]], 1)) + rLaplacian(b) #t(rmvnorm(1, g(c(X[[1]], 1)), V.matrix))
  for (t in 2:n){
    X[[t]] <- t(rmvnorm(1, f(c(X[[t-1]], t)), W.matrix))
    Y[[t]] <- g(c(X[[t]], t)) + rLaplacian(b) #t(rmvnorm(1, g(c(X[[t]], t)), V.matrix))
  }
  
  return(list(X=X, Y=Y))
}

## Calculates the Kalman gain
CalcKalmanGain <- function(R, G, Q){
  if(dim(Q)[1]==1){
    return(1/as.numeric(Q) * R %*% t(G))
  }else{
    return(R %*% t(G) %*% solve(Q))
  }
}

## Removes a dimension (time) from the jacobian in Extended/Ensemble filters
RemoveTime <- function(v){
  # assume the last variable is t
  # v must be a 1 by n matrix
  len <- dim(as.matrix(v))[2]
  return(as.matrix(v[, 1:(len-1)]))
}

## Returns the row-mean vector of x
CalculateMeanVector <- function(x){
  return(t(as.matrix(apply(x, 1, mean))))
}

## Returns x - xbar as a matrix
CalculateCenteredVector <- function(x, xbar, M){
  return(x - matrix(rep(xbar, M), byrow=FALSE, ncol=M))
}

## Returns the unbiased sample covariance matrix from X
CalculateSampleCovariance <- function(X, M){
  return(X %*% t(X) / (M-1))
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
  
  f.matrix <- as.matrix(f)
  g.matrix <- as.matrix(g)
  m0.matrix <- as.matrix(m0)
  C0.matrix <- as.matrix(C0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  n <- length(Y)
  p <- dim(m0)[1]
  
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
  
  m0.matrix <- as.matrix(m0)
  C0.matrix <- as.matrix(C0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  f.m0 <- as.matrix(f(c(m0.matrix, 1)))
  g.fm0 <- as.matrix(g(c(f.m0, 1)))
  df.m0 <- RemoveTime(pracma::jacobian(f, c(m0.matrix, 1)))
  dg.fm0 <- t(RemoveTime(pracma::jacobian(g, c(f.m0, 1))))
  
  n <- length(Y)
  p <- dim(dg.fm0)[2] # dim(m0)[1]
  
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
    df.m <- RemoveTime(pracma::jacobian(f, c(m[[t-1]], t)))
    dg.fm <- t(RemoveTime(pracma::jacobian(g, c(f.m, t))))
    
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
  
  p <- dim(m0)[1] # dimension of X
  size.ensemble <- dim(m0)[2] # size of ensemble
  q <- dim(Y[[1]])[1] # dimension of Y
  max.time <- length(Y) # number of time steps

  m0.matrix <- as.matrix(m0)
  C0.matrix <- as.matrix(C0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)

  ## define output lists
  a <- list()     # list of matrices, each matrix is M by p
  abar <- list()  # list of vectors
  m <- list()     # list of matrices, each matrix is M by p
  mbar <- list()  # list of vectors
  C <- list()
  R <- list()
  #Q <- list()
  Qprime <- list()
  K <- list()
  
  f.m0 <- matrix(0, nrow=p, ncol=size.ensemble)
  g.fm0 <- matrix(0, nrow=q, ncol=size.ensemble)
  for(i in 1:size.ensemble){
    f.m0[, i] <- as.matrix(f(c(m0.matrix[, i], 1))) + rmvnorm(1, t(t(numeric(p))), W.matrix) 
    g.fm0[, i] <- as.matrix(g(c(f.m0[, i], 1)))
  }
  
  # Perform initial calculation using initial conditions
  a[[1]] <- f.m0
  abar[[1]] <- CalculateMeanVector(a[[1]])
  Ea <- CalculateCenteredVector(a[[1]], abar[[1]], size.ensemble)
  R[[1]] <- CalculateSampleCovariance(Ea, size.ensemble)
  #Q[[1]] <- dg.fm0 %*% R[[1]] %*% t(dg.fm0) + V.matrix
  
  y <- matrix(0, nrow=q, ncol=size.ensemble)
  for(i in 1:size.ensemble){
    y[, i] <- rmvnorm(q, Y[[1]], V.matrix) 
  }
  #ybar <- CalculateMeanVector(y)
  #Ey <- CalculateCenteredVector(y, ybar, size.ensemble)
  #Qprime[[1]] <- CalculateSampleCovariance(Ey, size.ensemble)
  
  #K[[1]] <- Ea %*% t(Ey) %*% solve(Ey %*% t(Ey))
  #dg.fm <- t(RemoveTime(pracma::jacobian(g, c(a[[1]], 1))))
  #K[[1]] <- Ea %*% t(Ea) %*% dg.fm %*% solve(Qprime[[1]])/(size.ensemble-1)
  
  Ga <- matrix(0, nrow=q, ncol=size.ensemble)
  for(i in 1:size.ensemble){
    Ga[1, i] <- g(c(a[[1]][i], 1))
  }
  Ga.bar <- CalculateMeanVector(Ga)
  E.Ga <- CalculateCenteredVector(Ga, Ga.bar, size.ensemble)
  
  K[[1]] <- Ea %*% t(E.Ga) %*% solve(E.Ga %*% t(E.Ga)/(size.ensemble-1))/(size.ensemble-1)
  
  
  
  m[[1]] <- matrix(0, nrow=p, ncol=size.ensemble)
  for(i in 1:size.ensemble){
    m[[1]][, i] <- a[[1]][, i] + K[[1]] %*% (y[, i] - g.fm0[, i])
  }
  mbar[[1]] <- CalculateMeanVector(m[[1]])
  Em <- CalculateCenteredVector(m[[1]], mbar[[1]], size.ensemble)
  C[[1]] <- CalculateSampleCovariance(Em, size.ensemble)
  
  for(t in 2:max.time){
    f.m <- matrix(0, nrow=p, ncol=size.ensemble)
    g.fm <- matrix(0, nrow=q, ncol=size.ensemble)
    for(i in 1:size.ensemble){
      f.m[, i] <- as.matrix(f(c(m[[t-1]][, i], t))) + rmvnorm(1, t(t(numeric(p))), W.matrix) 
      g.fm[, i] <- g(c(f.m[, i], t))
    }
    
    # predict
    a[[t]] <- f.m
    abar[[t]] <- CalculateMeanVector(a[[t]])
    Ea <- CalculateCenteredVector(a[[t]], abar[[t]], size.ensemble)
    R[[t]] <- CalculateSampleCovariance(Ea, size.ensemble)
    
    # update
    y <- matrix(0, nrow=q, ncol=size.ensemble)
    for(i in 1:size.ensemble){
      y[, i] <- rmvnorm(q, Y[[t]], V.matrix) 
    }
    #ybar <- CalculateMeanVector(y)
    #Ey <- CalculateCenteredVector(y, ybar, size.ensemble)
    #Qprime[[t]] <- CalculateSampleCovariance(Ey, size.ensemble)
    
    #K[[t]] <- Ea %*% t(Ey) %*% solve(Ey %*% t(Ey))
    #dg.fm <- t(RemoveTime(pracma::jacobian(g, c(a[[t]], t))))
    #K[[t]] <- Ea %*% t(Ea) %*% dg.fm %*% solve(Qprime[[t]])/(size.ensemble-1)
    
    
    Ga <- matrix(0, nrow=q, ncol=size.ensemble)
    for(i in 1:size.ensemble){
      Ga[1, i] <- g(c(a[[t]][i], t))
    }
    Ga.bar <- CalculateMeanVector(Ga)
    E.Ga <- CalculateCenteredVector(Ga, Ga.bar, size.ensemble)
    
    K[[t]] <- Ea %*% t(E.Ga) %*% solve(E.Ga %*% t(E.Ga)/(size.ensemble-1))/(size.ensemble-1)
    
    m[[t]] <- matrix(0, nrow=p, ncol=size.ensemble)
    for(i in 1:size.ensemble){
      m[[t]][, i] <- a[[t]][, i] + K[[t]] %*% (y[, i] - g.fm[, i])
    }
    mbar[[t]] <- CalculateMeanVector(m[[t]])
    Em <- CalculateCenteredVector(m[[t]], mbar[[t]], size.ensemble)
    C[[t]] <- CalculateSampleCovariance(Em, size.ensemble)
  }
  
  return(list(abar=abar, mbar=mbar, C=C, R=R, K=K))#, Qprime=Qprime))
}

ApplyKalmanFilterThinned <- function(Y, f, g, m0, C0, W, V){
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
  
  f.matrix <- as.matrix(f)
  g.matrix <- as.matrix(g)
  m0.matrix <- as.matrix(m0)
  C0.matrix <- as.matrix(C0)
  W.matrix <- as.matrix(W)
  V.matrix <- as.matrix(V)
  
  n <- length(Y)
  p <- dim(m0)[1]
  
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
  if(!is.na(Y[[1]])){
    m[[1]] <- a[[1]] + K[[1]] %*% (Y[[1]] - g.matrix %*% a[[1]])
    C[[1]] <- (diag(p) -  K[[1]] %*% g.matrix) %*% R[[1]]
  }else{
    m[[1]] <- a[[1]]
    C[[1]] <- R[[1]]
  }

  
  for(t in 2:n){
    # predict
    a[[t]] <- f.matrix %*% m[[t-1]]
    R[[t]] <- f.matrix %*% C[[t-1]] %*% t(f.matrix) + W.matrix
    Q[[t]] <- g.matrix %*% R[[t]] %*% t(g.matrix) + V.matrix
    
    # update
    K[[t]] <- CalcKalmanGain(R[[t]], g.matrix, Q[[t]])
    if(!is.na(Y[[t]])){
      m[[t]] <- a[[t]] + K[[t]] %*% (Y[[t]] - g.matrix %*% a[[t]])
      C[[t]] <- (diag(p) -  K[[t]] %*% g.matrix) %*% R[[t]]
    }else{
      m[[t]] <- a[[t]]
      C[[t]] <- R[[t]]
    }
  }
  
  return(list(a=a, m=m, C=C, R=R, Q=Q, K=K))
}

ThinData <- function(data, p=0.95){
  # INPUTS:
  #   data: list of observed states
  #   p:    probability of keeping results
  # OUTPUTS:
  #   thinned.data: list of observed states with missing values
  
  thinned.data.list <- list()
  for(i in 1:length(data)){
    if(runif(1,0,1) < p){
      thinned.data.list[[i]] <- data[[i]]
    }else{
      thinned.data.list[[i]] <- NA
    }
  }
  
  return(thinned.data.list)
}



###########################################################

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
