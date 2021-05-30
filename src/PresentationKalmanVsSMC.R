####################
# Source files
####################
source("functions/Kalman.R")
source("functions/Plotting.R")
source("functions/SMC.R")

#########################################################################
# Kalman

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
plot(0:n, c(0, unlist(Y)), xlab="time", ylab="Obersations", main="Obersavtions Time Series")


# Initial estimates
m0 <- X0 # mean
C0 <- W # variance

#SavePlotToPNG("plots/KalmanStationary.png")
kalman.solution <- ApplyKalmanFilter(Y, A, B, m0, C0, W, V)
PlotKalmanSolution(unlist(X), 
                   unlist(kalman.solution$m), 
                   unlist(kalman.solution$C),
                   location = "topleft"
)
#dev.off()




########################################################################
# SMC

set.seed(2022)

EvolutionEquation <- function(x, t){
  return(x + rnorm(1, 0, 1))
}

# ObservationEquation <- function(x){
#   return(x + rnorm(1, 0, 1.5))
# }

ObservationLikelihood <- function(y, x, observation.variance=1.5^2){
  return(dnorm(y, x, sqrt(observation.variance)))
}


## GENERATE DATA

f <- function(x){
  return(x[1])
}

g <- function(x){
  return(x[1])
}

W <- 1
V <- 1

X0 <- 0

MaxT <- 25

data <- GenerateKalmanData(f, g, X0, W, V, MaxT)




## ANALYSIS

N <- 1000 # Number of particles

x0 <- rnorm(N, X0, sqrt(W))

solution.SMC <- ApplyBootstrapFilter(unlist(data$Y), 
                                     x0,
                                     EvolutionEquation, 
                                     ObservationLikelihood,
                                     N)

lines(1:MaxT, rowMeans(solution.SMC$particles), col="blue")
suppressWarnings(legend("topleft", 
                        legend=c("X(t) - Unobserved System",
                                 "Kalman Posterior Mean",
                                 "Particle Filter Mean"),
                        col=c("black", "red", "blue"), 
                        lty=c(0, 1, 1), 
                        pch=c(1, 26, 26)
))

pf.variance.estimates <- numeric(n)
pf.mean <- numeric(n)
for(i in 1:n){
  pf.variance.estimates[i] <- var(solution.SMC$particles[i, ])
  pf.mean[i] <- mean(solution.SMC$particles[i, ])
}
confidence <- 0.95
confidenceRange <- qnorm(1-(1-confidence)/2) * sqrt(pf.variance.estimates)
polygon(c(1:n, n:1), 
        c(pf.mean - confidenceRange, rev(pf.mean + confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:n, pf.mean - confidenceRange, col="blue", lty=2)
lines(1:n, pf.mean + confidenceRange, col="blue", lty=2)

MAE.PF <- mean(abs(unlist(X) - pf.mean))
MAE.PF
MAE.KF <- mean(abs(unlist(X) - unlist(kalman.solution$m)))
MAE.KF

##############
hist(solution.SMC$particles[20, ], xlab="Particles", freq=FALSE, main="Distribution of Particles at t=20")
lines(seq(-4, 4, length=300), dnorm(seq(-4, 4, length=300), mean(solution.SMC$particles[20, ]), sqrt(var(solution.SMC$particles[20, ]))), col="red")
