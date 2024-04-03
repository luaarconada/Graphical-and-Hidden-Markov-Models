rm(list=ls())
#
# Earthquake data.
#
if (!require("MixtureInf"))
  install.packages("MixtureInf")
library(MixtureInf)
data(earthquake)
help(earthquake)
#
# Bar chart of data.  
#
y <- earthquake$number
counts <- table(y)
plot(sort(unique(y)),counts/sum(counts),type='h',xlim=c(0,45),ylim=c(0,0.12),
     xlab="number of earthquakes",ylab="frequency")
ytick <- seq(0,0.12,by=0.02)
axis(side=2, at=ytick, labels = ytick)
#
# Does a single Poisson distribution look reasonable?
#
lines(sort(unique(y)),dpois(sort(unique(y)),mean(y)),
      type='h',lwd=2,lty=2,col="purple")
#
# Use MxtureInf to fit a mixture model.  We could try with k=1 or 2 also.
#
k <- 3
out <- pmle.pois(earthquake,k,1)
aic <- -2*out$"Penalized log-likelihood:"
aic
w <- out$"PMLE of mixing proportions:"
lambda <- out$"PMLE of component parameters:"
w
lambda
#
# Add the weighted component estimates to the plot.
#
plot(sort(unique(y)),counts/sum(counts),type='h',xlim=c(0,45),ylim=c(0,0.12),
     xlab="number of earthquakes",ylab="frequency")
ytick <- seq(0,0.12,by=0.02)
axis(side=2, at=ytick, labels = ytick)
for (j in 1:k){
  lines(c(0:45),w[j]*dpois(c(0:45),lambda[j]),type='h',col=j,lty=2,lwd=2)
}
#
# Calculate the full, fitted mixture distribution.
#
pred <- rep(0,46)
for (j in 1:k){
  pred <- pred+w[j]*dpois(c(0:45),lambda[j])
}  
#
# Plot with the real data again.
#
plot(sort(unique(y)),counts/sum(counts),type='h',xlim=c(0,45),ylim=c(0,0.12),
     xlab="number of earthquakes",ylab="frequency")
ytick <- seq(0,0.12,by=0.02)
axis(side=2, at=ytick, labels = ytick)
lines(c(0:45),pred,type='h',lwd=2,lty=2,col='purple')
#
# Time series plot of data.
#
x <- c(1900:2006)
plot(x,y,type='l',lwd=2,xlab="Year",ylab="Number of earthquakes")
#
# Look at the ACF of the earthquake data.  There is some autocorrelation.  A simple mixture isn't appropriate.
#
acf(earthquake,main="")
#
# Fit an HMM instead
#
if (!require("depmixS4"))
  install.packages("depmixS4")
library(depmixS4)
#
hmminit <- depmix(number~1,nstates=k,family=poisson(),data=earthquake)
set.seed(1)
hmmfit <- fit(hmminit)
hmmfit
#
# Extract parameters
#
params <- getpars(hmmfit)
pinit <- params[1:k]
pinit
trans <- t(matrix(params[(k+1):(k+k^2)],nrow=3))
trans
lambda <- exp(params[(k+k^2+1):(2*k+k^2)])
states <- depmixS4::posterior(hmmfit)
states
states <- states[,1]
#
# Plot with predicted mean values added.  The state changes are the point where the heights of the predictions change
#
plot(x,y,type='l',lwd=2,xlab="Year",ylab="Number of earthquakes",ylim=c(0,50))
lines(x,lambda[states],type='l',lwd=2,col='purple')
#
# Now calculate the steady state distribution of the MC.  The steady state weights are fairly close to the 
# independent mixture estimates.
#
library(markovchain)
ergodicmc <- new("markovchain",transitionMatrix=trans)
w <- steadyStates(ergodicmc)
w
out$"PMLE of mixing proportions:"
#
# Also compare the lambdas.  Again they are fairly close.
#
lambda
out$"PMLE of component parameters:"
#
# Now plot the marginal distribution again with the steady state component and predictive probabilities.  
# These are similar to the results for the independent model. 
#
plot(sort(unique(y)),counts/sum(counts),type='h',xlim=c(0,45),ylim=c(0,0.12),
     xlab="number of earthquakes",ylab="frequency")
ytick <- seq(0,0.12,by=0.02)
axis(side=2, at=ytick, labels = ytick)
for (j in 1:k){
  lines(c(0:45),w[j]*dpois(c(0:45),lambda[j]),type='h',col=j,lty=2,lwd=2)
}
pred <- rep(0,46)
for (j in 1:k){
  pred <- pred+w[j]*dpois(c(0:45),lambda[j])
}  
plot(sort(unique(y)),counts/sum(counts),type='h',xlim=c(0,45),ylim=c(0,0.12),
     xlab="number of earthquakes",ylab="frequency")
ytick <- seq(0,0.12,by=0.02)
axis(side=2, at=ytick, labels = ytick)
lines(c(0:45),pred,type='h',lwd=2,lty=2,col='purple')
