# Kalman Filter Example in one dimension
# For simplicity, let G(X(t))=1 for all t, hence the observed Y is just noisy X.
# Let W=1, V=2.
# 

library(mvtnorm)
library(expm)
library(plotly)

f <- matrix(c(2, 0, 0, 2), nrow=2, byrow=TRUE)

g <- matrix(c(1, 1), nrow=1)

ObserveY <- function(m, Q){
  return(rmvnorm(1, m, sqrtm(Q)[1]))
}

# Set variances
W <- matrix(c(0.05, 0, 0, 0.05), nrow=2, byrow=TRUE)
V <- 0.1

#initial conditions
m0 <- matrix(c(2, 2), nrow=2) # mean
C0 <- matrix(c(0.05, 0, 0, 0.05), nrow=2, byrow=TRUE) # variance

maxt <- 6
m <- list() #numeric(maxt+1)
m[[1]] <- m0
C <- list() #numeric(maxt+1)
C[[1]] <- C0
X <- list() #numeric(maxt)
Y <- list() #numeric(maxt)
R <- list() #numeric(maxt)
Q <- list() #numeric(maxt)
K <- list() #numeric(maxt)
prior <- list() #numeric(maxt)

for(t in 1:maxt){
  # predict
  R[[t]] <- f %*% C[[t]] %*% t(f) + W
  X[[t]] <- rmvnorm(1, f %*% m[[t]], R[[t]])
  prior[[t]] <- X[[t]]
  #R[t] <- f*C[t]*f + W
  Q[[t]] <- g %*% R[[t]] %*% t(g) + V
  
  # observe Y
  Y[[t]] <- rmvnorm(1, g %*% f %*% m[[t]], Q[[t]])
  
  # make observation and update
  K[[t]] <- g %*% R[[t]] / Q[[t]][1]
  m[[t+1]] <- f %*% m[[t]] + t(K[[t]] * (Y[[t]] - g %*% f %*% m[[t]])[1])
  C[[t+1]] <- (diag(2) - t(K[[t]]) %*% g) %*% R[[t]]
  X[[t]] <- rmvnorm(1, m[[t+1]], C[[t+1]])
}


data <- data.frame(x.axis=0:maxt, y.axis=0:maxt, true=2*2^(0:maxt), X1=c(m0[1], unlist(X)[2*(1:maxt)-1]),
                   X2=c(m0[2], unlist(X)[2*(1:maxt)]), Y=c(0, unlist(Y)))

fig <- plot_ly(data, x=~x.axis, y=~y.axis, z=~true, type = 'scatter3d', mode = 'lines')
fig <- fig %>% add_trace(data, x=~X1, y=~X2, z=~Y, type = 'scatter3d', mode = 'lines')
fig









