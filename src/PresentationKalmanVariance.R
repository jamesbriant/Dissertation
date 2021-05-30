####################
# Source files
####################
source("functions/Kalman.R")
source("functions/Plotting.R")


################################################################################
# Example: X-1D; Y-1D;
# stationary model

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

n <- 25
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y
#plot(0:n, c(0, unlist(Y)), xlab="time", ylab="Observations", main="Obseravtions Time Series")


# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

#SavePlotToPNG("plots/KalmanStationary.png")
kalmanSolution <- ApplyKalmanFilter(Y, A, B, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(kalmanSolution$m), 
                   unlist(kalmanSolution$C)
)
#dev.off()


test.plot.data <- function(i, data){
  plot(1:i, data[1:i], 
       xlim=c(1, length(data)), 
       ylim=1.1*c(min(data), max(data)), 
       pch=16, 
       xlab="Time",
       ylab="Observation",
       main="Observations collected in time")
}

AnimatePlot("animations/test1.mp4", test.plot.data, 15, unlist(Y)[1:15], auto.play=TRUE, width=700, height=400, res=100)





test.plot.solution <- function(i, data1){
  # INPUT
  #   X:  true DLM system states at times 1,...,T, vector
  #   m:  estimated system states, vector
  #   C:  estimated variance, vector
  #   confidence: size of confidence bands on estimate
  
  confidence <- 0.95
  X <- data1$X
  m <- data1$m
  C <- data1$C
  
  n <- length(X)
  plot(1:i, X[1:i], 
       xlab="Time", 
       ylab="System States",
       main="System State Estimates Updated with Time",
       xlim=c(1, n),
       ylim=c(min(X, m)-1, max(X, m)+1)
       )
  
  confidenceRange <- qnorm(1-(1-confidence)/2) * C[1:i]
  polygon(c(1:i, i:1), 
          c(m[1:i] - confidenceRange, rev(m[1:i] + confidenceRange)),
          col=rgb(1, 0, 0, 0.075),
          border=NA)
  lines(1:i, m[1:i] - confidenceRange, col="red", lty=2)
  lines(1:i, m[1:i] + confidenceRange, col="red", lty=2)
  lines(1:i, m[1:i], col="red", lty=1)
  
  #plot(1:length(X), X, add=TRUE)
  
  suppressWarnings(legend("topleft", 
                          legend=c("X(t) - Unobserved System",
                                   "m(t) - Posterior",
                                   paste0(100*confidence, "% confidence region")),
                          col=c("black", "red", "red"), 
                          lty=c(0, 1, 2), 
                          pch=c(1, 26, 26)
  ))
}

AnimatePlot("animations/test2.mp4", test.plot.solution, 15, 
            data.frame(X=unlist(X)[1:15], 
                       m=unlist(kalmanSolution$m)[1:15], 
                       C=unlist(kalmanSolution$C)[1:15]), 
            auto.play=TRUE, width=800, height=600, res=100)

################################################################################
# Example: X-1D; Y-1D;
# stationary model

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
W <- as.matrix(5)
V <- as.matrix(0.5)

set.seed(2022)

n <- 25
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y
#plot(0:n, c(0, unlist(Y)), xlab="time", ylab="Obersations", main="Obersavtions Time Series")


# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

#SavePlotToPNG("plots/KalmanStationary.png")
kalmanSolution <- ApplyKalmanFilter(Y, A, B, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(kalmanSolution$m), 
                   unlist(kalmanSolution$C)
)
#dev.off()