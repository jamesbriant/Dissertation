source("functions/SMC.R")
source("functions/Kalman.R")
source("functions/Plotting.R")

set.seed(2022)

EvolutionEquation <- function(x){
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

MaxT <- 8

data <- GenerateKalmanData(f, g, X0, W, V, MaxT)

  
  
  
## ANALYSIS

N <- 1000 # Number of particles

x0 <- rnorm(N, X0, sqrt(W))

test <- ApplyBootstrapFilter(unlist(data$Y), 
                             x0,
                             EvolutionEquation, 
                             ObservationLikelihood,
                             N)

PlotSMCSolution(test$particles.tilde, test$weights, 1)

#########################################################
# TESTING FOR PLOTTING FUNCTION

###############
# Create function parameter variables
###############

particles <- test$particles.tilde
weights <- test$weights
plotting.interval <- 2

###############
# Function implementation
###############

M <- dim(particles)[1]
N <- dim(particles)[2]

## weights do not have to be provided.
## Assume equal weights if not provided

if(identical(weights, FALSE)){
  weights = matrix(1/N, nrow=M, ncol=N)
}

# the times which are to be plotted
times <- seq(1, M, by=plotting.interval)
m <- length(times)

# collect the data at the designated times for plotting
particles.plot <- particles[times, ]
weights.plot <- weights[times, ]
weights.plot <- t(apply(weights.plot, 1, StandardiseWeights))

dens <- sapply(1:m, function(i) density(particles.plot[i, ], weights=weights.plot[i, ]))

data <- data.frame(
          x = rev(unlist(lapply(1:m, function(i) dens[, i]$x))),
          y = rev(unlist(lapply(1:m, function(i) dens[, i]$y))),
          z = factor(rep(times, each=length(dens[, 1]$x)))
        )

fig <- plot_ly(data, x=~x, y=~z, z=~y, type = 'scatter3d', mode = 'lines', color=~z)
fig


################################################################

#hist(test$particles[7, ])


t <- 7

dens.tilde <- density(test$particles.tilde[t, ])
plot(dens.tilde$x, dens.tilde$y, type="l", col="orange", ylim=c(0, 0.45))

abline(v=unlist(data$Y)[t], col="green")

dens.tilde <- density(test$particles.tilde[t, ], weights=StandardiseWeights(test$weights[t, ]))
lines(dens.tilde$x, dens.tilde$y, col="blue")

dens <- density(test$particles[t, ])
lines(dens$x, dens$y)

legend("topleft",
       legend=c(
         "prediction",
         "observation",
         "importance sample",
         "resampled"
       ),
       col=c("orange", "green", "blue", "black"),
       lty=c(1, 1, 1, 1)
       )



dens.tilde <- density(test$particles.tilde[t+1, ])
lines(dens.tilde$x, dens.tilde$y, col="orange", lty=2)

abline(v=unlist(data$Y)[t+1], col="green", lty=2)

dens <- density(test$particles[t+1, ])
lines(dens$x, dens$y, lty=2)



dens.tilde <- density(test$particles.tilde[t, ]*test$weights[t, ])
lines(dens.tilde$x, N*dens.tilde$y, col="blue")

abline(v=unlist(data$Y)[t-1], col="green")
abline(v=unlist(data$Y)[t], col="orange")
abline(v=unlist(data$Y)[t+1], col="red")











