raw.data <- read.csv("data/dataexport_20210311T102318.csv")
data <- as.numeric(as.character(raw.data$Basel[10:dim(raw.data)[1]]))

ts.plot(data)

acf(data)
pacf(data)

hist(data, freq=FALSE)
x <- seq(-5, 40, length=300)
lines(x, dnorm(x, mean(data), sd=sqrt(var(data))))






mu <- mean(data)
ts.plot(data)

t <- 0:length(data)
lines(t, mu*(sin((t-100)*2*pi/365)+1), col="red")

