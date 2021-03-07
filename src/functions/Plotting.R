#############
# LOADING PACKAGES
#############

LoadPlotly <- function(){
  if(!("plotly" %in% (.packages()))){
    library(plotly)
  }
}

LoadAV <- function(){
  if(!("av" %in% (.packages()))){
    library(av)
  }
}

#############
# Function definitions for plotting DLM data
#############

PlotKalmanSolution <- function(X, m, C, confidence=0.95){
  # INPUT
  #   X:  true DLM system states at times 1,...,T, vector
  #   m:  estimated system states, vector
  #   C:  estimated variance, vector
  #   confidence: size of confidence bands on estimate

  n <- length(X)
  plot(1:n, X, xlab="time", ylab="system state", ylim=c(min(X, m)-1, max(X, m)+1))

  confidenceRange <- qnorm(1-(1-confidence)/2) * C
  polygon(c(1:n, n:1), 
        c(m - confidenceRange, rev(m + confidenceRange)),
        col=rgb(1, 0, 0, 0.075),
        border=NA)
  lines(1:n, m - confidenceRange, col="red", lty=2)
  lines(1:n, m + confidenceRange, col="red", lty=2)
  lines(1:n, m, col="red", lty=1)

  #plot(1:length(X), X, add=TRUE)

  suppressWarnings(legend("topleft", 
        legend=c("X(t) - Unobserved System",
                  "m(t) - Posterior",
                  paste0(100*confidence, "% confidence region")),
        col=c("black", "red", "red"), 
        lty=c(0, 1, 2), 
        pch=c(1, 26, 26)
        ))
}

PlotSMCSolution <- function(particles, weights=FALSE, plotting.interval=1){
  LoadPlotly()
  
  m <- dim(particles)[1]
  N <- dim(particles)[2]
  
  ## weights do not have to be provided.
  ## Assume equal weights if not provided
  
  if(identical(weights, FALSE)){
    weights = matrix(1/N, nrow=m, ncol=N)
  }
  
  # the times which are to be plotted
  times <- seq(1, m, by=plotting.interval)
  
  # collect the data at the designated times for plotting
  particles.plot <- particles[times, ]
  weights.plot <- weights[times, ]
  weights.plot <- t(apply(weights.plot, 1, StandardiseWeights))
  m <- dim(particles.plot)[1]
  
  dens <- sapply(1:m, function(i) density(particles.plot[i, ], weights=weights.plot[i, ]))
  
  data <- data.frame(
    x = -unlist(lapply(1:m, function(i) dens[, i]$x)),
    y = unlist(lapply(1:m, function(i) dens[, i]$y)),
    z = factor(rep(times, each=length(dens[, 1]$x)))
  )
  
  fig <- plot_ly(data, x=~z, y=~x, z=~y, type = 'scatter3d', mode = 'lines', color=~z)
  fig <- fig %>% layout(
    title = "Test",
    scene = list(
      yaxis = list(title = "System State"), 
      zaxis = list(title = "Density"),
      xaxis = list(title = "Time",
                   autorange = "reversed")
    ))
  #fig <- fig %>% add_trace(x = 1:m, y=sapply(1:m, function(i) mean(data$x[which(data$z == i)])), z=rep(0, m), type="scatter3d", mode="lines")
  fig
}

#############
# SAVING IMAGES
#############

# Save the current plot. Units in mm
# smaller file, lower quality
SavePlotToJPG <- function(filename, w=150, h=90){
  jpeg(filename, width=w, height=h, units="mm", res=100)
}

# Save the current plot. Units in mm
# larger file, higher quality
#### USE THIS FUNCTION ####
SavePlotToPNG <- function(filename, w=150, h=90){
  png(filename, width=w, height=h, units="mm", res=100)
}


##############
# SAVING VIDEO
##############

AnimatePlot <- function(output, user.func, MaxT, data=FALSE, width=1280, height=720, res=108, auto.play=FALSE){
  LoadAV()
  
  png("plots/temp/temp%03d.png", width=width, height=height, res=res)
  for(i in 1:MaxT){
    user.func(i, data)
  }
  dev.off()
  
  png_files <- sprintf("plots/temp/temp%03d.png", 1:MaxT)
  av::av_encode_video(png_files, output, framerate = 1)
  
  if(auto.play == TRUE){
    utils::browseURL(output)
  }
}
