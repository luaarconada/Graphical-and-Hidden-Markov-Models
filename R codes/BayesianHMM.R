#
# Bayesian HMMs
#
rm(list=ls())
#
# Zucchini earthquake data showing the number of (magnitude 7) earthquakes every  
# year from 1900 to 2006.
#
y <- c(13, 14, 8, 10, 16, 26, 32, 27, 18, 32, 36, 24, 22, 23, 22, 18,
       25, 21, 21, 14, 8, 11, 14, 23, 18, 17, 19, 20, 22, 19, 13, 26,
       13, 14, 22, 24, 21, 22, 26, 21, 23, 24, 27, 41, 31, 27, 35, 26,
       28, 36, 39, 21, 17, 22, 17, 19, 15, 34, 10, 15, 22, 18, 15, 20,
       15, 22, 19, 16, 30, 27, 29, 23, 20, 16, 21, 21, 25, 16, 18, 15,
       18, 14, 10, 15, 8, 15, 6, 11, 8, 7, 18, 16, 13, 12, 13, 20,
       15, 16, 12, 18, 15, 16, 13, 15, 16, 11, 11)
x <- c(1900:2006)
earthquake <- data.frame(x = x, y = y)
n <- length(y)
#
# Fit a classical mixture model.  Use the initial values of lambda for the 
# Bayesian HMM.
#
if (!require("depmixS4")){install.packages("depmixS4")}
library(depmixS4)
earthquake <- data.frame(x = x, y = y)
temp <- mix(y ~ 1, nstates = 3, family = poisson(), data = earthquake)
fit3 <- fit(temp)
params <- getpars(fit3)
p3 <- unname(params[1:3])
lambda <- unname(exp(params[4:6]))
#
# Function for an HMM in OpenBUGS code.
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
k <- 3
data <- list("k" = k, "n" = n, "y" = y)
#
# Define the initial values.
#
lambda[1]
lambda[2]/lambda[1]
lambda[3]/lambda[2]
inits <- function(){
  list(lambda1 = 31.07, beta = c(0.5, 0.64, 0.66))
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
