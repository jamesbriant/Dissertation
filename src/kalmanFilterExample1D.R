# Kalman Filter Example in one dimension
# For simplicity, let G(X(t))=1 for all t, hence the observed Y is just noisy X.
# Let W=1, V=2.
# 

f <- 2
g <- 0.5

# Set variances
W <- 0.2
V <- 0.4

#initial conditions
m0 <- 1 # mean
C0 <- 0.2 # variance

maxt <- 6
m <- numeric(maxt+1)
m[1] <- m0
C <- numeric(maxt+1)
C[1] <- C0
X <- numeric(maxt)
Y <- numeric(maxt)
R <- numeric(maxt)
Q <- numeric(maxt)
K <- numeric(maxt)
prior <- numeric(maxt)

for(t in 1:maxt){
  # predict
  R[t] <- f*C[t]*f + W
  X[t] <- rnorm(1, f*m[t], R[t])
  prior[t] <- X[t]
  #R[t] <- f*C[t]*f + W
  Q[t] <- g*R[t]*g + V
  
  # observe Y
  Y[t] <- rnorm(1, g*f*m[t], sqrt(Q[t]))
  
  # make observation and update
  K[t] <- g*R[t]/Q[t]
  m[t+1] <- f*m[t] + K[t]*(Y[t] - g*f*m[t])
  C[t+1] <- (1 - g*K[t])*R[t]
  X[t] <- rnorm(1, m[t+1], C[t+1])
}

plot(0:maxt, 2^(0:maxt), type="l", xlab="time", ylab="state")
lines(1:maxt, X, col="red")
lines(1:maxt, Y, type="l", col="blue")
points(1:maxt, prior)
legend(0, 60, 
       legend=c("True Solution", 
                "X(t) - Unobserved System (posterior)",
                "Y(t) - Observed States", 
                "X(t) - Unobserved System (prior)"), 
       col=c("black", "red", "blue", "black"), 
       lty=c(1, 1, 1, 0), 
       pch=c(26, 26, 26, 1))





###########################################################################

# Kalman Filter Example in one dimension
# For simplicity, let G(X(t))=1 for all t, hence the observed Y is just noisy X.
# Let W=1, V=2.
# 

f <- function(t){
  return(2*cos(pi*t/6))
}
fprime <- function(t){
  return(-2*pi*(sin(pi*t/5)/5))
}
g <- 1

# Set variances
W <- 0.05
V <- 0.06

#initial conditions
m0 <- f(0) # mean
C0 <- 0.05 # variance

maxt <- 20
m <- numeric(maxt+1)
m[1] <- m0
C <- numeric(maxt+1)
C[1] <- C0
X <- numeric(maxt)
Y <- numeric(maxt)
R <- numeric(maxt)
Q <- numeric(maxt)
K <- numeric(maxt)
prior <- numeric(maxt)

for(t in 1:maxt){
  # predict
  ft <- f(t)
  fprimet <- fprime(t)
  R[t] <- ft*C[t]*ft + W
  X[t] <- rnorm(1, ft*m[t], R[t])
  prior[t] <- X[t]
  #R[t] <- f*C[t]*f + W
  Q[t] <- g*R[t]*g + V
  
  # observe Y
  Y[t] <- rnorm(1, g*ft*m[t], sqrt(Q[t]))
  
  # make observation and update
  K[t] <- g*R[t]/Q[t]
  m[t+1] <- ft*m[t] + K[t]*(Y[t] - g*ft*m[t])
  C[t+1] <- (1 - g*K[t])*R[t]
  X[t] <- rnorm(1, m[t+1], C[t+1])
}

plot(0:maxt, f(0:maxt), type="l", xlab="time", ylab="state")
lines(1:maxt, X, col="red")
lines(1:maxt, Y, type="l", col="blue")
points(1:maxt, prior)
legend(0, 60, 
       legend=c("True Solution", 
                "X(t) - Unobserved System (posterior)", 
                "Y(t) - Observed States", 
                "X(t) - Unobserved System (prior)"), 
       col=c("black", "red", "blue", "black"), 
       lty=c(1, 1, 1, 0), 
       pch=c(26, 26, 26, 1))













