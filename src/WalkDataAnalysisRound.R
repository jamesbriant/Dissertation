
library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(tidyverse)
library(ggmap)
library(mvtnorm)

source("functions/SMC.R")
source("functions/Plotting.R")


#####################
# Code
#####################

StandardiseTimes <- function(raw.times){
  # remove the dates and normalises the times. starts from 0.
  # units are seconds
  
  times1 <- mapply(sub, ".*T", "", raw.times)
  times.matrix <- matrix(0, nrow=length(times1), ncol=3)
  for(i in 1:length(times1)){
    times.matrix[i, ] <- as.vector(unlist(strsplit(times1[i], ":")))
    times.matrix[i, 3] <- mapply(sub, "Z.*", "", times.matrix[i, 3])
  }
  times.seconds.final <- 60*60*as.numeric(times.matrix[, 1])
  times.seconds.final <- times.seconds.final + 60*as.numeric(times.matrix[, 2])
  times.seconds.final <- times.seconds.final + round(as.numeric(times.matrix[, 3]))
  times.seconds.final <- times.seconds.final - times.seconds.final[1]
  
  return(times.seconds.final)
}

MaxS <- 50

# read in the raw data
#raw.data <- read.csv("data/20210320135655.csv") # 1st walk
raw.data <- read.csv("data/20210322154634.csv")#[1:(MaxS+1), ] # circle walk
#raw.data <- read.csv("data/20210322162522.csv")[1:(MaxS+1), ] # straight walk

MaxS <- length(raw.data$lat)-1

times <- StandardiseTimes(raw.data$time)

df <- data.frame(x = raw.data$lon,
                 y = raw.data$lat)

cov.df <- cov(df)

## Work out the initial bearing here


####################################################
# Dynamics functions

#sorted
EvolutionEquation <- function(X, s, psi, theta){
  x     <- X[1]
  y     <- X[2]
  r     <- exp(psi)*(times[s+1] - times[s])

  if(theta < -pi/2 || theta > pi/2){
    dx <- r*cos(theta - pi/2)
    dy <- r*sin(theta - pi/2)
    
    if(theta < -pi/2){
      x.next <- x + dx
      y.next <- y - dy
    }else{
      x.next <- x + dx
      y.next <- y - dx
    }
  }else{
    dx <- r*sin(theta)
    dy <- r*cos(theta)
    
    if(theta < 0){
      x.next <- x + dx
      y.next <- y + dy
    }else{
      x.next <- x + dx
      y.next <- y + dx
    }
  }
  
  return(c(x.next, y.next))# + SystemNoise())
    # should there be system noise here????
}

#sorted
SystemNoise <- function(){
  return(c(rmvnorm(1, c(0,0), sigma=cov.df/100)))
}

ThetaCheckDomain <- function(theta){
  # check theta is still defined on [-pi, pi).
  # Correct any points that are not on this domain.
  
  theta.checked <- numeric(length(theta))
  
  for(i in 1:length(theta)){
    if(theta[i] >= pi){
      theta.checked[i] <- theta[i] - 2*pi
    }else if(theta[i] < -pi){
      theta.checked[i] <- theta[i] + 2*pi
    }else{
      theta.checked[i] <- theta[i]
    }
  }
  
  return(theta.checked)
}

CalculateNewTheta <- function(X.new, X.old){
  if(X.new[2] > X.old[2]){
    # if we moved upwards
    theta <- atan((X.new[1] - X.old[1])/(X.new[2] - X.old[2]))
  }else{
    theta <- atan((X.new[1] - X.old[1])/(X.new[2] - X.old[2])) + pi/2
  }
  
  # if(X.new[1] < X.old[1]){
  #   theta <- -theta
  # }
  
  return(ThetaCheckDomain(theta))
  #return(theta)
}

#sorted
# ObservationLikelihood <- function(y, x){
#   return(pmvnorm(lower=y,
#                  upper=Inf,
#                  mean=x,
#                  sigma=cov.df/1000)[1]
#          )
# }
ObservationLikelihood <- function(y, x){
  return(dmvnorm(y,
                 mean=x,
                 sigma=cov.df)[1]
  )
}



###############
# ANLAYSIS
###############

N <- 100 # Number of particles

data.list <- list()
for(t in 1:(MaxS+1)){
  data.list[[t]] <- c(df$x[t], df$y[t])
}

x0 <- t(rmvnorm(N, c(df$x[1], df$y[1]), (cov.df/1000)*diag(2)))
psi0 <- log(t(rnorm(N,
                    sqrt((data.list[[2]][1]-data.list[[1]][1])^2 + (data.list[[2]][2]-data.list[[1]][2])^2)/times[2],
                    0.001^2)))
theta0 <- CalculateNewTheta(data.list[[2]], data.list[[1]])

sol <- ApplyParamLearnAuxFilter(data.list[2:(MaxS+1)], 
                              x0, psi0, theta0,
                              EvolutionEquation, 
                              SystemNoise,
                              ObservationLikelihood
                              )



x.sol <- numeric(MaxS)
y.sol <- numeric(MaxS)
for(t in 1:MaxS){
  x.sol[t] <- mean(sol$particles[[t]][1, ])
  y.sol[t] <- mean(sol$particles[[t]][2, ])
}
df.sol <- data.frame(x = x.sol,
                     y = y.sol)

# get google cloud key
source("functions/GoogleKeyInfo.R")
MY_KEY <- GetKey()
ggmap::register_google(key = MY_KEY)

ggmap(get_googlemap(center = c(mean(raw.data$lon), 
                               mean(raw.data$lat)),
                    zoom = 16, 
                    size = c(640, 320),
                    scale = 2,
                    maptype ='terrain',
                    color = 'color',
                    markers = df.sol,
                    path = df.sol)
)



