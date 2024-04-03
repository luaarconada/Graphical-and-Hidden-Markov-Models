#
# Bayesian Stochastic volatility using MCMC & particle filtering.
#
# SV model:
#
# y_t ~ N(0,h_t)
# h_t = mu + phi(h_{t-1}-mu), sigma^2)
# h_0 ~ N(0,sigma^2/(1-phi^2))
#
rm(list=ls())
if(!require(stochvol)){install.packages("stochvol")}
if(!require(nimble)){install.packages("nimble")}
if(!require(nimbleSMC)){install.packages("nimbleSMC")}
#
# 1) MCMC approach
#
library(stochvol)
#
# Load data
#
data("exrates")
#
# Plot time series of Euro dollar prices and calculate and plot (demeaned) log returns.
#
ret <- logret(exrates$USD, demean = TRUE)
par(mfrow = c(1, 2), mar = c(1.9, 1.9, 1.9, 0.5), mgp = c(2, 0.6, 0))
plot(exrates$date, exrates$USD, type = "l",main = "Price of 1 EUR in USD")
plot(exrates$date[-1], ret, type = "l", main = "Demeaned log returns")
#
# Call MCMC sampler
#
res <- svsample(ret, priormu = c(-10, 1), priorphi = c(20, 1.1),priorsigma = 0.1)
#
# Look at trace plots.
#
par(mfrow = c(3, 1))
paratraceplot(res)
#
# Look at summary statistics and plot parameter estimates.
#
summary(res, showlatent = FALSE)
par(mfrow=c(1,3))
paradensplot(res)
#
# Fitted volatility estimates.
#
par(mfrow=c(1,1))
volplot(res, dates = exrates$date[-1])
#
# Adding out of sample predictions.
#
volplot(res, forecast = 365, dates = exrates$date[-1])
#
# 2) Fitting using sequential Monte Carlo
#
# You need to have Rtools installed to be able to run this.
#
library(nimble)
library(nimbleSMC)
stochVCode <- nimbleCode({
  x[1] ~ dnorm(phi * x0, sigmaSquaredInv)
  y[1] ~ dnorm(0, var = betaSquared * exp(x[1]))
  for(t in 2:T){
    x[t] ~ dnorm(phi * x[t-1], sigmaSquaredInv)
    y[t] ~ dnorm(0, var = betaSquared * exp(x[t]))
  }
  x0 ~ dnorm(1, sigmaSquaredInv)
  phi <- 2 * phiStar - 1
  phiStar ~ dbeta(18, 1)
  sigmaSquaredInv ~ dgamma(5, 20)
  betaSquared <- 1 / betaSquaredInv
  betaSquaredInv ~ dgamma(5, 20)
  })
stochVolModel <- nimbleModel(code = stochVCode, name = 'stochVol',constants = list(T = 67), data = list(y = ret),
  inits = list(betaSquaredInv = 2.785, phi = .9702,sigmaSquaredInv = 31.561))
stochVolModel$getNodeNames(includeData=FALSE)
CstochVolModel <- compileNimble(stochVolModel)
stochVolLiuWestFilter <- buildLiuWestFilter(model = stochVolModel,nodes = 'x', params = c( 'betaSquaredInv', 'phiStar',
  'sigmaSquaredInv'),control = list(saveAll=TRUE))
CstochVolLiuWestFilter <- compileNimble(stochVolLiuWestFilter,project = stochVolModel)
CstochVolLiuWestFilter$run(10000)
#
# Estimated parameter distributions
#
sigmaSquaredSamples <- 1 / as.matrix(CstochVolLiuWestFilter$mvEWSamples,'sigmaSquaredInv')
hist(sigmaSquaredSamples, main ='', xlab = 'sigma2')
betaSquaredSamples <- 1/as.matrix(CstochVolLiuWestFilter$mvEWSamples,'betaSquaredInv')
hist(betaSquaredSamples, main ='', xlab = 'beta2')
phiSamples <- 2*as.matrix(CstochVolLiuWestFilter$mvEWSamples,'phiStar')-1
hist(phiSamples, main ='', xlab = 'phi')


