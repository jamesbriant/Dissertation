# Returns the standardised vector of vec, ie. sum(|vec|) = 1
StandardiseWeights <- function(vec){
  return(vec/sum(sqrt(vec^2)))
}

# The resample step required for the Bootstrap Filter.
# Returns the new sample, each particle has equal weighting
Resample <- function(x.vec, weights){
  N <- length(x.vec)
  
  sample.counts <- rmultinom(1, N, weights)[, 1]
  
  return(rep(x.vec, sample.counts))
}

# input weights can be non-standardised
EffectiveSampleSize <- function(weights){
  return(1/sum(StandardiseWeights(weights)^2))
}

ApplyBootstrapFilter <- function(data, x0, Evolve, ObservationLikelihood, 
                                 N=1000){
  # data:                   observations made. Vector of length MaxT
  # x0:                     initial sample
  # Evolve:                 evolution equation, INCLUDING noise, MUST accept x and t
  # ObservationLikelihood:  observation likelihood function
  # N:                      Number of particles to use
  
  MaxT <- length(data)

  weights <- matrix(0, nrow=MaxT, ncol=N) # not standardised!!

  particles       <- matrix(0, nrow=MaxT+1, ncol=N) # each particle has weight 1/N, includes data from t=0,...,N
  particles.tilde <- matrix(0, nrow=MaxT, ncol=N) # each particle has weight weights[t,i], includes data from t=1,...,N
  
  particles[1, ] <- x0

  for(t in 1:MaxT){
    # Evolution
    
    particles.tilde[t, ] <- sapply(particles[t, ], Evolve, t)
    
    # Calculate the weights
    weights[t, ] <- mapply(ObservationLikelihood, data[t], particles.tilde[t, ])
    
    # Re-sample particles
    particles[t+1, ] <- Resample(particles.tilde[t, ], weights[t, ])
  }
  
  return(list("particles" = particles[2:(MaxT+1), ], 
              "particles.tilde" = particles.tilde,
              "weights" = weights))
}

ApplyBootstrapFilterNaive <- function(data, x0, Evolve, ObservationLikelihood, 
                                 N=1000, K=4){
  # data:                   observations made. Vector of length MaxT
  # x0:                     initial sample
  # Evolve:                 evolution equation, INCLUDING noise, MUST accept x and t
  # ObservationLikelihood:  observation likelihood function
  # N:                      Number of particles to use
  # K:                      Resample every Kth step
  
  MaxT <- length(data)
  
  weights.prior     <- matrix(0, nrow=MaxT+1, ncol=N) # t=1,...,N
  weights.posterior <- matrix(0, nrow=MaxT, ncol=N) # t=1,...,N
  
  particles.prior       <- matrix(0, nrow=MaxT+1, ncol=N) # each particle has weight w_t, includes data from t=1,...,N
  particles.posterior   <- matrix(0, nrow=MaxT, ncol=N) # each particle has weight w_t, includes data from t=1,...,N
  
  particles.prior[1, ]  <- x0
  weights.prior[1, ]    <- rep(1/(length(x0)), times=length(x0))
  
  for(t in 1:MaxT){
    # Evolution
    particles.posterior[t, ]  <- sapply(particles.prior[t, ], Evolve, t)
    weights.posterior[t, ]    <- weights.prior[t, ] * mapply(ObservationLikelihood, data[t], particles.prior[t, ])
    
    if(t %% K == 0){
      #resample
      particles.prior[t+1, ]  <- Resample(particles.posterior[t, ], weights.posterior[t, ])
      weights.prior[t+1, ]    <- rep(1/(length(x0)), times=length(x0))
    }else{
      #don't resample
      particles.prior[t+1, ]  <- particles.posterior[t, ]
      weights.prior[t+1, ]    <- weights.posterior[t, ]
    }
  }
  
  return(list("particles.prior" = particles.prior, 
              "particles.posterior" = particles.posterior,
              "weights.prior" = weights.prior,
              "weights.posterior" = weights.posterior)
         )
}



ApplyBootstrapFilterEff <- function(data, x0, Evolve, ObservationLikelihood, 
                                      N=1000, N0=FALSE){
  # data:                   observations made. Vector of length MaxT
  # x0:                     initial sample
  # Evolve:                 evolution equation, INCLUDING noise, MUST accept x and t
  # ObservationLikelihood:  observation likelihood function
  # N:                      Number of particles to use
  # N0:                     Minimum effective sample size allowed before resampling
  
  if(identical(N0, FALSE)){
    N0 = sqrt(length(x0))
  }
  MaxT <- length(data)
  
  weights.prior     <- matrix(0, nrow=MaxT+1, ncol=N) # t=1,...,N
  weights.posterior <- matrix(0, nrow=MaxT, ncol=N) # t=1,...,N
  
  particles.prior       <- matrix(0, nrow=MaxT+1, ncol=N) # each particle has weight w_t, includes data from t=1,...,N
  particles.posterior   <- matrix(0, nrow=MaxT, ncol=N) # each particle has weight w_t, includes data from t=1,...,N
  
  particles.prior[1, ]  <- x0
  weights.prior[1, ]    <- rep(1/(length(x0)), times=length(x0))
  
  effective.size <- numeric(MaxT)
  
  for(t in 1:MaxT){
    # Evolution
    particles.posterior[t, ]  <- sapply(particles.prior[t, ], Evolve, t)
    weights.posterior[t, ]    <- weights.prior[t, ] * mapply(ObservationLikelihood, data[t], particles.prior[t, ])
    
    effective.size[t] <- EffectiveSampleSize(weights.posterior[t, ])
    
    if(effective.size[t] <= N0){
      #resample
      particles.prior[t+1, ]  <- Resample(particles.posterior[t, ], weights.posterior[t, ])
      weights.prior[t+1, ]    <- rep(1/(length(x0)), times=length(x0))
    }else{
      #don't resample
      particles.prior[t+1, ]  <- particles.posterior[t, ]
      weights.prior[t+1, ]    <- weights.posterior[t, ]
    }
  }
  
  return(list("particles.prior" = particles.prior, 
              "particles.posterior" = particles.posterior,
              "weights.prior" = weights.prior,
              "weights.posterior" = weights.posterior,
              "effective.sample.sizes" = effective.size)
         )
}


ApplyAuxiliaryParticleFilter <- function(data, x0, Evolve, SystemNoise,
                                         ObservationLikelihood, N=1000){
  # data:                   observations made. Vector of length MaxT
  # x0:                     initial sample
  # Evolve:                 evolution function, EXCLUDING noise, MUST accept x and t
  # SystemNoise:            evolution noise function
  # ObservationLikelihood:  observation likelihood function
  # N:                      Number of particles to use
  
  MaxT <- length(data)
  
  weights <- matrix(0, nrow=MaxT, ncol=N) # not standardised!!
  
  particles       <- matrix(0, nrow=MaxT+1, ncol=N) # each particle has weight 1/N, includes data from t=0,...,N
  particles.tilde <- matrix(0, nrow=MaxT, ncol=N) # each particle has weight weights[t,i], includes data from t=1,...,N
  
  particles[1, ] <- x0
  
  for(t in 1:MaxT){
    # Evolution
    particles.evolved <- sapply(particles[t, ], Evolve, t) # noise not applied
    
    expected.likelihood <- mapply(ObservationLikelihood, data[t], particles.evolved)
    
    #particles.k <- rmultinom(particles[t, ], 1:N, expected.likelihood)[, 1]
    particles.k <- Resample(particles[t, ], expected.likelihood) # probability changes when weights are not equal
    particles.tilde[t, ] <- particles.k + SystemNoise(N=N)
    
    
    
    # Calculate the weights
    weights[t, ] <- mapply(ObservationLikelihood, data[t], particles.tilde[t, ])/expected.likelihood
    
    # Re-sample particles
    particles[t+1, ] <- Resample(particles.tilde[t, ], weights[t, ])
  }
  
  return(list("particles" = particles[2:(MaxT+1), ], 
              "particles.tilde" = particles.tilde,
              "weights" = weights)
         )
}



ApplyParamLearnAuxFilter <- function(data, x0, psi0, theta, Evolve, SystemNoise,
                                     ObservationLikelihood, a=0.99){
  # data:                   observations made. list of length MaxT
  # x0:                     initial sample, matrix 
  # psi0:                   initial parameter sample, matrix
  # theta:                  initial value for theta
  # Evolve:                 evolution function, EXCLUDING noise, MUST accept x and t
  # SystemNoise:            evolution noise function
  # ObservationLikelihood:  observation likelihood function
  # a:                      `a` parameter where a^2 + h^2 = 1
  
  if(!(dim(x0)[2] == dim(psi0)[2])){
    print("Number of particles in psi0 and x0 must be the same!")
    return(0)
  }
  
  h2 <- 1- a^2
  
  no.particles <- dim(x0)[2]
  no.states <- dim(x0)[1]
  no.params <- dim(psi0)[1]
  MaxT <- length(data)
  
  weights <- matrix(0, nrow=MaxT, ncol=no.particles) # not standardised!!
  
  particles       <- list() # each particle has weight 1/N, includes data for t=0,...,T
  particles.tilde <- list() # each particle has weight weights[t,i], includes data for t=1,...,T
  psi             <- list() # params at time t=0,...,T
  m               <- list() # list of adjust psi particles
  theta.vec       <- numeric()
  
  #particles       <- matrix(0, nrow=MaxT+1, ncol=no.particles) # each particle has weight 1/N, includes data from t=0,...,N
  #particles.tilde <- matrix(0, nrow=MaxT, ncol=no.particles) # each particle has weight weights[t,i], includes data from t=1,...,N
  
  particles[[1]]  <- x0
  psi[[1]]        <- psi0
  theta.vec[1]    <- theta
  
  for(t in 1:MaxT){
    
    psi.bar <- mean(psi[[t]][1, ])
    Psi     <- var(psi[[t]][1,])
    m[[t]]  <- matrix(0, nrow=no.params, ncol=no.particles)
    for(i in 1:no.params){
      m[[t]][i, ] <- a*psi[[t]][i, ] + (1-a)*psi.bar
    }
    
    # Evolution
    particles.evolved <- matrix(0, nrow=no.states, ncol=no.particles)
    expected.likelihood <- numeric(no.particles)
    for(i in 1:no.particles){
      particles.evolved[, i] <- Evolve(particles[[t]][, i], t, m[[t]][, i], theta.vec[t])
      expected.likelihood[i] <- ObservationLikelihood(data[[t]], particles.evolved[, i])
    }
    expected.likelihood <- StandardiseWeights(expected.likelihood)
    #particles.evolved <- mapply(Evolve, particles[[t]], t, m[[t]], theta.vec[t]) # noise not applied
    #expected.likelihood <- mapply(ObservationLikelihood, data[[t]], particles.evolved)
    
    particles.tilde[[t]] <- matrix(0, nrow=no.states, ncol=no.particles)
    
    psi[[t+1]] <- matrix(0, nrow=no.params, ncol=no.particles)
    I <- rmultinom(1, no.particles, expected.likelihood)
    i <- 0
    for(k in 1:no.particles){
      if(I[k] > 0){
        for(j in 1:I[k]){
          i <- i + 1
          psi[[t+1]][1, i] <- rnorm(1, m[[t]][1, k], sqrt(h2*Psi))
          particles.tilde[[t]][, i] <- Evolve(particles[[t]][, k], t, psi[[t+1]][1, k], theta.vec[t]) + SystemNoise()
          
          # Calculate the weights
          weights[t, i] <- ObservationLikelihood(data[[t]], particles.tilde[[t]][, k])/expected.likelihood[k]
        }
      }
    }
    
    #weights[t, ] <- mapply(ObservationLikelihood, data[t], particles.tilde[t, ])/expected.likelihood
    
    # Re-sample particles and psi
    I <- rmultinom(1, no.particles, weights[t, ])
    i <- 0
    particles[[t+1]] <- matrix(0, nrow=no.states, ncol=no.particles)
    for(k in 1:no.particles){
      if(I[k] > 0){
        for(j in 1:I[k]){
          i <- i + 1
          particles[[t+1]][, i] <- particles.tilde[[t]][, k]
          psi[[t+1]][1, i] <- psi[[t+1]][1, k]
        }
      }
      
    }
    #particles[[t+1]] <- rep(particles.tilde[[t]], I[, 1])
    #psi[[t+1]] <- rep(psi[[t+1]], I[, 1])
    
    #calculate the next theta
    current.posterior <- rowMeans(particles[[t+1]])
    previous.posterior <- rowMeans(particles[[t]])
    theta.vec[t+1] <- CalculateNewTheta(current.posterior, previous.posterior)
  }
  
  return(list("particles" = particles[2:(MaxT+1)], 
              "particles.tilde" = particles.tilde,
              "psi" = psi,
              "weights" = weights,
              "theta" = theta.vec)
  )
}













