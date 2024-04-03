#
# Bayesian HMMs
#
rm(list=ls())
if (!require("MixtureInf"))
  install.packages("MixtureInf")
library(MixtureInf)
data(earthquake)
y <- earthquake$number
#
# Use MxtureInf to fit a classical mixture model.  Use the initial values of lambda for the Bayesian HMM.
#
k <- 3
out <- pmle.pois(earthquake,k,1)
lambda <- out$"PMLE of component parameters:"
lambda <- sort(lambda,decreasing=TRUE)
#
# Function of OpenBUGS code.
#
bugsHMM <- function(){
  s[1] ~ dcat(P0[])
  y[1] ~ dpois(lambda[s[1]])
  for (j in 2:n){
    s[j] ~ dcat(P.mat[s[j-1],])
    y[j] ~ dpois(lambda[s[j]])
  }
  for(i in 1:k){
    alpha[i] <- 1
    P0[i] <- 1/k
    P.mat[i,1:k] ~ ddirich(alpha[])
  }
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
data <- list("k"=k,"n"=n,"y"=y)
#
# Define the initial values.
#
lambda[1]
lambda[2]/lambda[1]
lambda[3]/lambda[2]
inits <- function(){
  list(lambda1=31.07,beta=c(0.5,0.64,0.66))
}
#
# Setting the parameters to be monitored
#
params <- c("lambda","P.mat","s")
#
# Run OpenBUGS
#
if (!require(R2OpenBUGS)){ 
  install.packages(R2OpenBUGS) 
}
library(R2OpenBUGS)
res <- bugs(data,inits,parameters.to.save=params,n.iter=10000,n.burnin=2000,
            n.chain=3,model.file=bugsHMM,codaPkg=FALSE,
            OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe")
cc <- res$sims.list
lambda <- cc$lambda
#
# Convergence checks
#
ts.plot(lambda[,1],ylab="lambda[1]")
ts.plot(cumsum(lambda[,1])/c(1:length(lambda[,1])),ylab="E[lambda[1]]")
acf(lambda[,1],main="lambda[1]")
#
# Posterior mean
#
elambda <- apply(lambda,2,mean)
elambda
#
# Posterior density plots
#
hist(lambda[,1],probability=TRUE,xlab="lambda[1]",ylab='f',main="Posterior density of lambda[1]")
lines(density(lambda[,1]),lwd=2,col='red')
hist(lambda[,2],probability=TRUE,xlab="lambda[2]",ylab='f',main="Posterior density of lambda[2]")
lines(density(lambda[,2]),lwd=2,col='red')
hist(lambda[,3],probability=TRUE,xlab="lambda[1]",ylab='f',main="Posterior density of lambda[3]")
lines(density(lambda[,3]),lwd=2,col='red')
#
# Also calculate the posterior mean transition matrix
#
P.mat <- cc$P.mat
eP.mat <- apply(P.mat,2:3,mean)
eP.mat
#
# Posterior mean values for each state
#
s <- cc$s
es <- apply(s,2,mean)
es
#
# Calculate predictive state probabilities and individual posterior state estimates.
#
ps <- matrix(rep(0,107*k),nrow=107)
maxstate <- rep(0,107)
for (i in 1:107){
  for (j in 1:k){
    ps[i,j] <- sum(s[,i]==j)
  }
  maxstate[i] <- which(ps[i,]==max(ps[i,]))
}
ps <- ps/dim(s)[1]
ps
maxstate
