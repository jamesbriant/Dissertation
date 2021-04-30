#####################
# Load Libraries
#####################

library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(tidyverse)
library(ggmap)

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

MaxS <- 36

# read in the raw data
#raw.data <- read.csv("data/20210320135655.csv") # 1st walk
#raw.data <- read.csv("data/20210322154634.csv") # circle walk
raw.data <- read.csv("data/20210322162522.csv")[1:(MaxS+1), ] # straight walk

times <- StandardiseTimes(raw.data$time)

df <- data.frame(x = raw.data$lon,
                 y = raw.data$lat)

velocity.x <- (df$x[MaxS+1] - df$x[1])/times[MaxS]
velocity.y <- (df$y[MaxS+1] - df$y[1])/times[MaxS]

W.x <- 1e-9
W.y <- W.x







#############
# x-direction
EvolutionEquation.x <- function(x.previous, s){
  return(velocity.x*(times[s+1] - times[s]) + x.previous + SystemNoise.x())
}

SystemNoise.x <- function(N=1){
  return(rnorm(N, 0, sqrt(W.x)))
}

ObservationLikelihood.x <- function(y, x, observation.variance=var(df$x)){
  return(dnorm(y, x, sqrt(observation.variance)))
}

#############
# y-direction
EvolutionEquation.y <- function(y.previous, s){
  return(velocity.y*(times[s+1] - times[s]) + y.previous + SystemNoise.y())
}

SystemNoise.y <- function(N=1){
  return(rnorm(N, 0, sqrt(W.y)))
}

ObservationLikelihood.y <- function(y, x, observation.variance=var(df$y)){
  return(dnorm(y, x, sqrt(observation.variance)))
}


## ANALYSIS

N <- 1000 # Number of particles

x0 <- rnorm(N, df$x[1], sqrt(1e-9))
y0 <- rnorm(N, df$y[1], sqrt(1e-9))



sol.x <- ApplyBootstrapFilter(df$x[2:(MaxS+1)], 
                             x0,
                             EvolutionEquation.x, 
                             ObservationLikelihood.x,
                             N=N)

sol.y <- ApplyBootstrapFilter(df$y[2:(MaxS+1)],
                             y0,
                             EvolutionEquation.y,
                             ObservationLikelihood.y,
                             N=N)

PlotSMCSolution(sol.y$particles.tilde, sol.y$weights, 3)


df.sol <- data.frame(x = rowMeans(sol.x$particles),
                     y = rowMeans(sol.y$particles))

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




