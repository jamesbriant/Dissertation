source("functions/SMC.R")
source("functions/Kalman.R")
source("functions/Plotting.R")

set.seed(2022)

EvolutionEquation <- function(x, t){
  return(x)
}

SystemNoise <- function(N=1){
  return(rnorm(N, 0, 1))
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

######
# TRY FORCING EXTREME DATA
#data$Y[[19]] <- data$X[[19]]-4


## ANALYSIS

N <- 20 # Number of particles

x0 <- rnorm(N, X0, sqrt(W))

set.seed(2022)
solution.SMC.aux <- ApplyAuxiliaryParticleFilter(unlist(data$Y), 
                                                 x0,
                                                 EvolutionEquation, 
                                                 SystemNoise,
                                                 ObservationLikelihood,
                                                 N=N)


apf.variance.estimates <- numeric(MaxT)
apf.mean <- numeric(MaxT)
# for(i in 1:MaxT){
#   pf.mean[i] <- mean(solution.SMC$particles[i, ])
#   pf.variance.estimates[i] <- var(solution.SMC$particles[i, ])
# }
for(i in 1:MaxT){
  apf.mean[i] <- sum(solution.SMC.aux$particles.tilde[i, ] * StandardiseWeights(solution.SMC.aux$weights[i, ]))
  apf.variance.estimates[i] <- sum(solution.SMC.aux$particles.tilde[i, ]^2 * StandardiseWeights(solution.SMC.aux$weights[i, ])) + apf.mean[i]^2
}
apf.confidence <- 0.95
apf.confidenceRange <- qnorm(1-(1-apf.confidence)/2) * sqrt(apf.variance.estimates)



#################
# Normal bootstrap

EvolutionEquation <- function(x, t){
  return(x + rnorm(1, 0, 1))
}

# ObservationEquation <- function(x){
#   return(x + rnorm(1, 0, 1.5))
# }

ObservationLikelihood <- function(y, x, observation.variance=1.5^2){
  return(dnorm(y, x, sqrt(observation.variance)))
}

set.seed(2022)
solution.SMC <- ApplyBootstrapFilter(unlist(data$Y), 
                                     x0,
                                     EvolutionEquation, 
                                     ObservationLikelihood,
                                     N)

pf.variance.estimates <- numeric(MaxT)
pf.mean <- numeric(MaxT)
# for(i in 1:MaxT){
#   pf.mean[i] <- mean(solution.SMC$particles[i, ])
#   pf.variance.estimates[i] <- var(solution.SMC$particles[i, ])
# }
for(i in 1:MaxT){
  pf.mean[i] <- sum(solution.SMC$particles.tilde[i, ] * StandardiseWeights(solution.SMC$weights[i, ]))
  pf.variance.estimates[i] <- sum(solution.SMC$particles.tilde[i, ]^2 * StandardiseWeights(solution.SMC$weights[i, ])) + pf.mean[i]^2
}
pf.confidence <- 0.95
pf.confidenceRange <- qnorm(1-(1-pf.confidence)/2) * sqrt(pf.variance.estimates)






################## 
# Plotting
plot(1:MaxT, 
     unlist(data$X), 
     xlab="time", 
     ylab="system state", 
     ylim=c(min(unlist(data$X), 
                pf.mean - pf.confidenceRange, 
                apf.mean - apf.confidenceRange
                ), 
            max(unlist(data$X), 
                pf.mean + pf.confidenceRange, 
                apf.mean + apf.confidenceRange
                )
            ),
     main="Auxiliary Particle Filter vs Normal Particle Filter"
     )

lines(1:MaxT, rowMeans(solution.SMC$particles), col="blue")
polygon(c(1:MaxT, MaxT:1), 
        c(pf.mean - pf.confidenceRange, rev(pf.mean + pf.confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:MaxT, pf.mean - pf.confidenceRange, col="blue", lty=2)
lines(1:MaxT, pf.mean + pf.confidenceRange, col="blue", lty=2)

lines(1:MaxT, apf.mean, col="magenta")
polygon(c(1:MaxT, MaxT:1), 
        c(apf.mean - apf.confidenceRange, rev(apf.mean + apf.confidenceRange)),
        col=rgb(1, 0.1, 0.8, 0.075),
        border=NA)
lines(1:MaxT, apf.mean - apf.confidenceRange, col="magenta", lty=2)
lines(1:MaxT, apf.mean + apf.confidenceRange, col="magenta", lty=2)
suppressWarnings(legend("topleft", 
                        legend=c("X(t) - Unobserved System",
                                 "Normal Particle Filter",
                                 "Auxiliary Particle Filter"),
                        col=c("black", "blue", "magenta"), 
                        lty=c(0, 1, 1), 
                        pch=c(1, 26, 26)
))



MAE.PF <- mean(abs(unlist(data$X) - pf.mean))
MAE.PF
MAE.APF <- mean(abs(unlist(data$X) - unlist(apf.mean)))
MAE.APF

hist(unlist(data$X) - unlist(data$Y))






#PlotSMCSolution(test$particles.tilde, test$weights, 1)
