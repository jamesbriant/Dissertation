#########################################################################
# Sample f=Beta(2,5) from g=Unif(0,1) using Rejection Sampling

source('functions/rejection-sampling-functions.R')

data <- drawSample(10000)
hist(data, freq=FALSE, xlab="x", main="Beta(2,5) Distribution from 2500 Draws")
x <- seq(from=0, to=1, length.out=250)
lines(x, dbeta(x, 2, 5), col="red")
legend(0.5, 2.4, legend=c("Rejection Sampling", "True distribution"), pch="-", col=c("black", "red"))
rug(data)


