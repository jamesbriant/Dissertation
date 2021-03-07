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

W <- 5
V <- 0.5

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
