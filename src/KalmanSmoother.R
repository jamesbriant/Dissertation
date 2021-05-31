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
#plot(0:n, c(0, unlist(Y)), xlab="time", ylab="Obersations", main="Obersavtions Time Series")


# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

#SavePlotToPNG("plots/KalmanStationary.png")
kalman.solution <- ApplyKalmanFilter(Y, A, B, m0, C0, W, V)

kalman.smoothed.sol <- ApplyKalmanSmoother(A, 
                                           kalman.solution$a,
                                           kalman.solution$m,
                                           kalman.solution$R,
                                           kalman.solution$C)


PlotKalmanSolution(unlist(X), 
                   unlist(kalman.solution$m), 
                   unlist(kalman.solution$C)
)
lines(1:n, kalman.smoothed.sol$s, col="blue")
confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * sqrt(unlist(kalman.smoothed.sol$S))
polygon(c(1:n, n:1), 
        c(unlist(kalman.smoothed.sol$s) - confidenceRange, rev(unlist(kalman.smoothed.sol$s) + confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:n, unlist(kalman.smoothed.sol$s) - confidenceRange, col="blue", lty=2)
lines(1:n, unlist(kalman.smoothed.sol$s) + confidenceRange, col="blue", lty=2)

suppressWarnings(legend("topleft", 
                        legend=c("X(t) - Unobserved System",
                                 "Kalman Filtered Mean",
                                 "Kalman Smoothed Mean"),
                        col=c("black", "red", "blue"), 
                        lty=c(0, 1, 1), 
                        pch=c(1, 26, 26)
))

MAE.KS <- mean(abs(unlist(X) - unlist(kalman.smoothed.sol$s)))
MAE.KS
MAE.KF <- mean(abs(unlist(X) - unlist(kalman.solution$m)))
MAE.KF


##############################################################################
# Different variance smoother
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
kalman.solution <- ApplyKalmanFilter(Y, A, B, m0, C0, W, V)

kalman.smoothed.sol <- ApplyKalmanSmoother(A, 
                                           kalman.solution$a,
                                           kalman.solution$m,
                                           kalman.solution$R,
                                           kalman.solution$C)


PlotKalmanSolution(unlist(X), 
                   unlist(kalman.solution$m), 
                   unlist(kalman.solution$C)
)
lines(1:n, kalman.smoothed.sol$s, col="blue")
confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * sqrt(unlist(kalman.smoothed.sol$S))
polygon(c(1:n, n:1), 
        c(unlist(kalman.smoothed.sol$s) - confidenceRange, rev(unlist(kalman.smoothed.sol$s) + confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:n, unlist(kalman.smoothed.sol$s) - confidenceRange, col="blue", lty=2)
lines(1:n, unlist(kalman.smoothed.sol$s) + confidenceRange, col="blue", lty=2)

suppressWarnings(legend("topleft", 
                        legend=c("X(t) - Unobserved System",
                                 "Kalman Filtered Mean",
                                 "Kalman Smoothed Mean"),
                        col=c("black", "red", "blue"), 
                        lty=c(0, 1, 1), 
                        pch=c(1, 26, 26)
))

MAE.KS <- mean(abs(unlist(X) - unlist(kalman.smoothed.sol$s)))
MAE.KS
MAE.KF <- mean(abs(unlist(X) - unlist(kalman.solution$m)))
MAE.KF

################################################################################
# Missing data example

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
Y.thinned <- ThinData(Y, 0.85)
#plot(0:n, c(0, unlist(Y)), xlab="time", ylab="Obersations", main="Obersavtions Time Series")


# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

kalman.solution <- ApplyKalmanFilterThinned(Y.thinned, A, B, m0, C0, W, V)

kalman.smoothed.sol <- ApplyKalmanSmoother(A, 
                                           kalman.solution$a,
                                           kalman.solution$m,
                                           kalman.solution$R,
                                           kalman.solution$C)

PlotKalmanSolution(unlist(X), 
                   unlist(kalman.solution$m), 
                   unlist(kalman.solution$C),
                   location = "bottomright"
)
lines(1:n, kalman.smoothed.sol$s, col="blue")
confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * sqrt(unlist(kalman.smoothed.sol$S))
polygon(c(1:n, n:1), 
        c(unlist(kalman.smoothed.sol$s) - confidenceRange, rev(unlist(kalman.smoothed.sol$s) + confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:n, unlist(kalman.smoothed.sol$s) - confidenceRange, col="blue", lty=2)
lines(1:n, unlist(kalman.smoothed.sol$s) + confidenceRange, col="blue", lty=2)

points(which(is.na(unlist(Y.thinned))), unlist(X)[which(is.na(unlist(Y.thinned)))], col="magenta", pch=15)
suppressWarnings(legend("bottomright", 
                        legend=c("X(t) - Unobserved System",
                                 "Kalman Filtered Mean",
                                 "Kalman Smoothed Mean",
                                 "Missing observation"),
                        col=c("black", "red", "blue", "magenta"), 
                        lty=c(0, 1, 1, 0), 
                        pch=c(1, 26, 26, 15)
))
