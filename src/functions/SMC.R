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

# ApplyParticleFilter <- function(data, evolve, observe, N=1000){
#   MaxT <- length(data)
#   
#   weights <- matrix(0, nrow=MaxT, ncol=N) # not standardised!!
#   weights[1, ] <- rep(1/N, N)
#   
#   particles <- numeric(N)
#   
#   for(i in 1:MaxT){
#     
#   }
# }


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
    particles.k <- Resample(particles[t, ], expected.likelihood)
    particles.tilde[t, ] <- particles.k + SystemNoise(N=N)
    
    
    
    # Calculate the weights
    weights[t, ] <- mapply(ObservationLikelihood, data[t], particles.tilde[t, ])/expected.likelihood
    
    # Re-sample particles
    particles[t+1, ] <- Resample(particles.tilde[t, ], weights[t, ])
  }
  
  return(list("particles" = particles[2:(MaxT+1), ], 
              "particles.tilde" = particles.tilde,
              "weights" = weights))
  
  
  
  
  
  
}
















