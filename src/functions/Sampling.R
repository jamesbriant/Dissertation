##### U ~ Unif(0,1) implementation

## returns a single sample drawn from Uniform(0,1)
SimulateU <- function(){
  return(runif(1, 0, 1))
}



##### sampling distribution (g in the literature) implementation

## returns PDF of Uniform(a,b) at x
## x can be a vector
GetUniformPDF <- function(x, a, b){
  if(x > a && x < b){
    return(1/(b - a))
  }
  
  return(0)
}

## returns sample of size n drawn from Uniform(a,b)
GetUniformSample <- function(n, a, b){
  return(runif(n, a, b))
}


##### Empirical Supremum Rejection Sampling

## Implements the Empirical Supremum Rejection Sampling algorithm to find an
## estimate for c.
FindC <- function(target.PDF, candidate.PDF, candidate.sampler, 
                    return.all.candidates=FALSE){
  # input
  #   target.PDF:  function - target distribution PDF
  #   candidate.PDF:  function - candidate distribution PDF
  #   candidate.sampler:  function - Draws sample of size 1 from candidate distribution
  
  c_estimates <- c(1)
  
  ## time (in iterations) since the c_estimate changed
  last.change <- 0
  i <- 1
  
  while(last.change < 200){
    last.change <- last.change + 1
    i <- i + 1
    
    U <- SimulateU()
    X <- candidate.sampler()
    ratio <- target.PDF(X)/candidate.PDF(X)
    
    ## We don't care if the sample is accepted or rejected for estimating c
    
    c_estimates[i] <- max(c_estimates[i-1], ratio)
    
    if(c_estimates[i] != c_estimates[i-1]){
      last.change <- 0
    }
  }
  
  if(return.all.candidates == TRUE){
    return(c_estimates)
  }else{
    return(max(c_estimates))
  }
}

###### Sample Generator

## Generates a sample of size n from target.PDF using candidates from 
## candidate.PDF

GenerateSample <- function(n, target.PDF, candidate.PDF, candidate.sampler, 
                            c=0){
  # input
  #   n:    size of target sample
  #   target.PDF:  function - target distribution PDF
  #   candidate.PDF:  function - candidate distribution PDF
  #   candidate.sampler:  function - Draws sample of size 1 from candidate distribution
  #   c:  
  if(c == 0){
    c <- FindC(target.PDF, candidate.PDF, candidate.sampler)
  }
  
  target.sample <- numeric(n)
  
  for(i in 1:n){
    accept <- FALSE
    
    while(accept == FALSE){
      U <- SimulateU()
      X <- candidate.sampler()
      ratio <- target.PDF(X)/candidate.PDF(X)
      
      if(U*c < ratio){
        accept <- TRUE
      }
    }
    
    target.sample[i] <- X
  }
  
  return(target.sample)
}


