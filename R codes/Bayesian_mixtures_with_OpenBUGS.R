#
# Bayesian mixtures 
#
rm(list=ls())
#
# Earthquake data.
#
if (!require("MixtureInf"))
  install.packages("MixtureInf")
library(MixtureInf)
data(earthquake)
y <- earthquake$number
#
# Use MxtureInf to fit a classical mixture model.  Use these for initial values for Bayesian mixture.
#
k <- 3
out <- pmle.pois(earthquake,k,1)
aic <- -2*out$"Penalized log-likelihood:"
aic
w <- out$"PMLE of mixing proportions:"
lambda <- out$"PMLE of component parameters:"
w
lambda
cc <- sort(lambda,decreasing=TRUE,index.return=TRUE)
lambda <- cc$x
w <- w[cc$ix]
#
# Bayesian mixture model BUGS code
#
bugsmixture <- function(){
  for (j in 1:n){
    y[j] ~ dpois(lambda[s[j]])
    s[j] ~ dcat(P[])
  }
  P[1:k] ~ ddirich(alpha[])
  lambda1 ~ dgamma(0.001,0.001)
  lambda[1] <- lambda1
  beta[1] ~ dbeta(0.001,0.001)
  for (i in 2:k){
    beta[i]   ~ dbeta(0.001,0.001)
    lambda[i] <-  beta[i]*lambda[i-1]
  }
}  
#
# Define the data for the BUGS code
#
n <- length(y)
data <- list("k"=k,"n"=n,"y"=y,"alpha"=rep(1,k))
#
# Define the initial values.
#
lambda[2]/lambda[1]
lambda[3]/lambda[2]
inits <- function(){
  list(lambda1=31.07,beta=c(0.5,0.64,0.66),P=c(0.147,0.552,0.301))
}
#
# Setting the parameters to be monitored
#
params <- c("lambda","P")
#
# Run OpenBUGS
#
if (!require(R2OpenBUGS)){ 
  install.packages(R2OpenBUGS) 
}
library(R2OpenBUGS)
res <- bugs(data,inits,parameters.to.save=params,n.iter=10000,n.burnin=2000,
            n.chain=3,model.file=bugsmixture,codaPkg=FALSE,
            OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe")
#
chain <- res$sims.list
lambda <- chain$lambda
P <- chain$P
#
# Convergence checks
#
iters <- length(lambda[,1])
ts.plot(lambda[,1])
ts.plot(cumsum(lambda[,1])/c(1:iters))
acf(lambda[,1])
#
# Posterior summary statistics and density plots.
#
print(res)
plot(density(lambda[,1]),xlab="lambda[1]",main="Posterior density of lambda[1]")
plot(density(lambda[,2]),xlab="lambda[2]",main="Posterior density of lambda[2]")
plot(density(lambda[,3]),xlab="lambda[3]",main="Posterior density of lambda[3]")
#
# Posterior predictive density.
#
pred <- rep(0,46)
for (i in 1:iters){
  for (j in 1:k){
    pred <- pred+P[i,j]*dpois(c(0:45),lambda[i,j])
  }
}
pred <- pred/iters
#
# Plot data and predictive density.
#
counts <- table(y)
plot(sort(unique(y)),counts/sum(counts),type='h',xlim=c(0,45),ylim=c(0,0.12),
     xlab="number of earthquakes",ylab="frequency")
ytick <- seq(0,0.12,by=0.02)
axis(side=2, at=ytick, labels = ytick)
lines(c(0:45),pred,type='h',lwd=2,lty=2,col='purple')
