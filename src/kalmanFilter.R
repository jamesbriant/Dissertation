####################
# Libraries
####################
library(plotly) # for 3d plots

####################
# Source files
####################
source("functions/Kalman.R")
source("functions/Plotting.R")

################################################################################
# Example: X-1D; Y-1D;

A <- as.matrix(2)
B <- as.matrix(1)

f <- function(x){
  return(A %*% x[1])
}
g <- function(x){
  return(B %*% x[1])
}

X0 <- as.matrix(1)

# Set variances
W <- as.matrix(3)
V <- as.matrix(2)

n <- 7
data <- GenerateKalmanData(f, g, X0, W, V, n)
X <- data$X
Y <- data$Y



# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

#SavePlotToPNG("test1.png")
kalmanSolution <- ApplyKalmanFilter(Y, A, B, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                  unlist(kalmanSolution$m), 
                  unlist(kalmanSolution$C)
)
#dev.off()

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
#dev.off()

################################################################################
# Example: X-2D; Y-1D; 

f <- function(x){
  return(matrix(c(2, 0, 0, 2), nrow=2, byrow=TRUE))
}
g <- function(x){
  return(matrix(c(1, 1), nrow=1))
}

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




################################################################################
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

















