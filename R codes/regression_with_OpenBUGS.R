#
# Multiple linear regression via OpenBUGS
#
rm(list=ls())
if (!require(R2OpenBUGS)){ 
  install.packages(R2OpenBUGS) 
}
library(R2OpenBUGS)
#
# Regression data.
#
if (!require(dplyr)){ 
  install.packages(dplyr) 
}
library(dplyr)
data(mtcars)
help(mtcars)
head(mtcars)
#
# Just use a subset of the data.
#
df <- mtcars[,c(1,3,5,8,9)]
#
# Look at plots of mpg against the different regressors.
#
plot(df$disp,df$mpg,xlab="disp",ylab="mpg")
plot(df$drat,df$mpg,xlab="drat",ylab="mpg")
plot(df$vs,df$mpg,xlab="vs",ylab="mpg")
plot(df$am,df$mpg,xlab="am",ylab="mpg")
#
# Perform classical regression.
#
model <- mpg ~ disp + drat + vs + am
class_reg <- lm(model,data=df)
summary(class_reg)
#
# Bayesian regression: define a function with the BUGS regression model
#
bugsregression <- function(){
    for( i in 1 : n ) {
		  mpg[i] ~ dnorm(mu[i], tau)
      mu[i] <- beta0 + beta1 * disp[i] + beta2 * drat[i] + beta3 * vs[i] + beta4 * am[i]
    }
    beta0 ~ dnorm(0.0, 1.0E-6)
    beta1 ~ dnorm(0.0, 1.0E-6)
    beta2 ~ dnorm(0.0, 1.0E-6)
    beta3 ~ dnorm(0.0, 1.0E-6)
    beta4 ~ dnorm(0.0, 1.0E-6)
    tau ~ dgamma(0.001, 0.001)
    sigma <- 1 / sqrt(tau)
}
#
# Define the data for Bayesian regression.
# 
n <- length(df$mpg)
data <- list("n"=n,"mpg"=df$mpg,"disp"=df$disp,"drat"=df$drat,"vs"=df$vs,"am"=df$am)
#
# Define the initial values.
#
inits <- function(){
  list(beta0=rnorm(1,20,1),beta1=rnorm(1,0,1),beta2=rnorm(1,0,1),beta3=rnorm(1,0,1),beta4=rnorm(1,0,1),tau=rgamma(1,1,1))
}
#
# Setting the parameters to be monitored
#
params <- c("beta0","beta1","beta2","beta3","beta4","tau","sigma")
#
# Running the code: We may have to change the OpenBUGS.pgm (path to OpenBUGS) to do this.  You can do this by looking at the 
# properties of the OpenBUGS icon and copying the path.
#
res <- bugs(data,inits,parameters.to.save=params,n.iter=10000,n.burnin=2000,
            n.chain=3,model.file=bugsregression,codaPkg=FALSE,
            OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe")
#
# Looking at the output
#
print(res)
plot(res)
chain <- res$sims.list
#
# Checking convergence
#
ts.plot(chain$beta0,ylab="beta0")
ts.plot(cumsum(chain$beta0)/c(1:length(chain$beta0)),ylab="beta0")
acf(chain$beta0,main="beta0")
#
# Posterior density plot
#
hist(chain$beta0,prob=TRUE,xlab="beta0",main="")
lines(density(chain$beta0),col="red",lwd=2)
#
# Prediction.
#
# disp = 150, drat = 3.8, vs = 1, am = 0.
#
mpgpred <- chain$beta0+chain$beta1*150+chain$beta2*3.8+chain$beta3+rnorm(length(chain$sigma),0,chain$sigma)
hist(mpgpred,prob=TRUE,xlab="mpg",main="")
lines(density(mpgpred),col="green",lwd=2)
mean(mpgpred)
#
# Classical model selection.
# 
summary(class_reg)
#
# Classical stepwise regression using AIC.
#
if (!require(MASS)){ 
  install.packages(MASS) 
}
library(MASS)
step <- stepAIC(class_reg,direction="both")
#
# Check the DIC of the Bayesian full model fit
#
res$DIC
#
# Now fit the model without drac.
#
bugsregression2 <- function(){
  for( i in 1 : n ) {
    mpg[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 * disp[i] + beta3 * vs[i] + beta4 * am[i]
  }
  beta0 ~ dnorm(0.0, 1.0E-6)
  beta1 ~ dnorm(0.0, 1.0E-6)
  beta3 ~ dnorm(0.0, 1.0E-6)
  beta4 ~ dnorm(0.0, 1.0E-6)
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1 / sqrt(tau)
}
data2 <- list("n"=n,"mpg"=df$mpg,"disp"=df$disp,"vs"=df$vs,"am"=df$am)
inits2 <- function(){
  list(beta0=rnorm(1,20,1),beta1=rnorm(1,0,1),beta3=rnorm(1,0,1),beta4=rnorm(1,0,1),tau=rgamma(1,1,1))
}
params2 <- c("beta0","beta1","beta3","beta4","tau","sigma")
res2 <- bugs(data2,inits2,parameters.to.save=params2,n.iter=10000,n.burnin=2000,
            n.chain=3,model.file=bugsregression2,codaPkg=FALSE,
            OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe")
res2$DIC