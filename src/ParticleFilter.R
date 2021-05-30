source("functions/SMC.R")
source("functions/Kalman.R")
source("functions/Plotting.R")

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
PlotSMCSolution(test$particles)

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


################################################################################
# Effects of resampling on the posterior distribution

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

N <- 100 # Number of particles

x0 <- rnorm(N, X0, sqrt(W))

# PF solution
set.seed(2022)
solution.SMC <- ApplyBootstrapFilter(unlist(data$Y), 
                                     x0,
                                     EvolutionEquation, 
                                     ObservationLikelihood,
                                     N)
par(mfrow=c(2,2))
for(i in 1:4){
  t <- i+10
  plot(density(solution.SMC$particles[t, ]), col="red", xlab="x", main=paste0("t=",t))
  lines(density(solution.SMC$particles.tilde[t, ], weights=StandardiseWeights(solution.SMC$weights[t, ])))
}
par(mfrow=c(1,1))

################################################################################
# Effects of sample size on statistical performance of the filter

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


N <- c(10, 80, 150, 250)
mean.estimates <- list()
variance.estimates <- list()
for(i in 1:length(N)){
  x0 <- rnorm(N[i], X0, sqrt(W))

  solution.SMC <- ApplyBootstrapFilter(unlist(data$Y), 
                             x0,
                             EvolutionEquation, 
                             ObservationLikelihood,
                             N[i])
  mean.estimates[[i]] <- numeric(MaxT)
  variance.estimates[[i]] <- numeric(MaxT)
  for(t in 1:MaxT){
    mean.estimates[[i]][t] <- mean(solution.SMC$particles[t, ])
    variance.estimates[[i]][t] <- var(solution.SMC$particles[t, ])
  }
}

for(i in 1:length(N)){
  if(i==1){
    plot(1:n, variance.estimates[[i]], col=i, type="o", xlab="time", ylab="sample variance", main="Sample Variance for Different Sized Samples")
    leg <- c(paste0("N=", N[[i]]))
  }else{
    lines(1:n, variance.estimates[[i]], col=i, type="o")
    leg <- c(leg, paste0("N=", N[i]))
  }
}
legend("topright",
       legend = leg,
       pch=26,
       lty=1,
       col=1:length(N))

############################################################################
# Resampling methods - Naive


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

N <- 100 # Number of particles

x0 <- rnorm(N, X0, sqrt(W))

# PF solution
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


# Naive solution
set.seed(2022)
solution.SMC.naive <- ApplyBootstrapFilterNaive(unlist(data$Y), 
                             x0,
                             EvolutionEquation, 
                             ObservationLikelihood,
                             N, K=5)

pfn.variance.estimates <- numeric(MaxT)
pfn.mean <- numeric(MaxT)
# for(i in 1:MaxT){
#   pfn.mean[i] <- sum(solution.SMC.naive$particles.posterior[i, ] * StandardiseWeights(solution.SMC.naive$weights.posterior))
#   pfn.variance.estimates[i] <- sum(solution.SMC.naive$particles.posterior[i, ]^2 * StandardiseWeights(solution.SMC.naive$weights.posterior)) + pfn.mean[i]^2
# }
for(i in 1:MaxT){
  pfn.mean[i] <- sum(solution.SMC.naive$particles.prior[i+1, ] * StandardiseWeights(solution.SMC.naive$weights.prior[i+1, ]))
  pfn.variance.estimates[i] <- sum(solution.SMC.naive$particles.prior[i+1, ]^2 * StandardiseWeights(solution.SMC.naive$weights.prior[i+1, ])) + pfn.mean[i]^2
}
pfn.confidence <- 0.95
pfn.confidenceRange <- qnorm(1-(1-pfn.confidence)/2) * sqrt(pfn.variance.estimates)

### Plotting
plot(1:MaxT, 
     unlist(data$X), 
     xlab="time", 
     ylab="system state", 
     ylim=c(min(unlist(data$X), 
                pf.mean - pf.confidenceRange, 
                pfn.mean - pfn.confidenceRange
                ), 
            max(unlist(data$X), 
                pf.mean + pf.confidenceRange, 
                pfn.mean + pfn.confidenceRange
                )
            ),
     main="Naive Resampling"
     )

lines(1:MaxT, rowMeans(solution.SMC$particles), col="blue")
polygon(c(1:MaxT, MaxT:1), 
        c(pf.mean - pf.confidenceRange, rev(pf.mean + pf.confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:MaxT, pf.mean - pf.confidenceRange, col="blue", lty=2)
lines(1:MaxT, pf.mean + pf.confidenceRange, col="blue", lty=2)

lines(1:MaxT, pfn.mean, col="magenta")
polygon(c(1:MaxT, MaxT:1), 
        c(pfn.mean - pfn.confidenceRange, rev(pfn.mean + pfn.confidenceRange)),
        col=rgb(1, 0.1, 0.8, 0.075),
        border=NA)
lines(1:MaxT, pfn.mean - pfn.confidenceRange, col="magenta", lty=2)
lines(1:MaxT, pfn.mean + pfn.confidenceRange, col="magenta", lty=2)
suppressWarnings(legend("topleft", 
                        legend=c("X(t) - Unobserved System",
                                 "Resampling K=1",
                                 "Resampling K=5"),
                        col=c("black", "blue", "magenta"), 
                        lty=c(0, 1, 1), 
                        pch=c(1, 26, 26)
))

############################################################################
# Resampling methods - Effective Sample size


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

N <- 20 # Number of particles

x0 <- rnorm(N, X0, sqrt(W))

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

set.seed(2022)
solution.SMC.eff <- ApplyBootstrapFilterEff(unlist(data$Y), 
                                               x0,
                                               EvolutionEquation, 
                                               ObservationLikelihood,
                                               N, N0=7)

pfe.variance.estimates <- numeric(MaxT)
pfe.mean <- numeric(MaxT)
# for(i in 1:MaxT){
#   pfn.mean[i] <- sum(solution.SMC.eff$particles.posterior[i, ] * StandardiseWeights(solution.SMC.eff$weights.posterior))
#   pfn.variance.estimates[i] <- sum(solution.SMC.eff$particles.posterior[i, ]^2 * StandardiseWeights(solution.SMC.eff$weights.posterior)) + pfe.mean[i]^2
# }
for(i in 1:MaxT){
  pfe.mean[i] <- sum(solution.SMC.eff$particles.prior[i+1, ] * StandardiseWeights(solution.SMC.eff$weights.prior[i+1, ]))
  pfe.variance.estimates[i] <- sum(solution.SMC.eff$particles.prior[i+1, ]^2 * StandardiseWeights(solution.SMC.eff$weights.prior[i+1, ])) + pfe.mean[i]^2
}
pfe.confidence <- 0.95
pfe.confidenceRange <- qnorm(1-(1-pfe.confidence)/2) * sqrt(pfe.variance.estimates)

### Plotting
plot(1:MaxT, 
     unlist(data$X), 
     xlab="time", 
     ylab="system state", 
     ylim=c(min(unlist(data$X), 
                pf.mean - pf.confidenceRange, 
                pfe.mean - pfe.confidenceRange
     ), 
     max(unlist(data$X), 
         pf.mean + pf.confidenceRange, 
         pfe.mean + pfe.confidenceRange
     )
     ),
     main="Resampling - Effective Sample Size, N=20"
)

lines(1:MaxT, rowMeans(solution.SMC$particles), col="blue")
polygon(c(1:MaxT, MaxT:1), 
        c(pf.mean - pf.confidenceRange, rev(pf.mean + pf.confidenceRange)),
        col=rgb(0, 0, 1, 0.075),
        border=NA)
lines(1:MaxT, pf.mean - pf.confidenceRange, col="blue", lty=2)
lines(1:MaxT, pf.mean + pf.confidenceRange, col="blue", lty=2)

lines(1:MaxT, pfe.mean, col="magenta")
polygon(c(1:MaxT, MaxT:1), 
        c(pfe.mean - pfe.confidenceRange, rev(pfe.mean + pfe.confidenceRange)),
        col=rgb(1, 0.1, 0.8, 0.075),
        border=NA)
lines(1:MaxT, pfe.mean - pfe.confidenceRange, col="magenta", lty=2)
lines(1:MaxT, pfe.mean + pfe.confidenceRange, col="magenta", lty=2)
abline(v=7, col="green")
abline(v=13, col="green")
abline(v=22, col="green")
suppressWarnings(legend("topleft", 
                        legend=c("X(t) - Unobserved System",
                                 "Naive Resampling, K=1",
                                 "Minimum Effective Size = 7"),
                        col=c("black", "blue", "magenta"), 
                        lty=c(0, 1, 1), 
                        pch=c(1, 26, 26)
))