#
# Continuous output HMM.
#
rm(list=ls())
if (!require("MixtureInf")){install.packages("MixtureInf")}
if (!require("depmixS4")){install.packages("depmixS4")}
data(speed)
help(speed)
head(speed)
#
# Draw a histogram of the response time data.
#
y <- speed$rt
hist(y, freq = FALSE, xlab = "response times", ylab = "f", 
     main = "Histogram of response time data")
#
# Fit different sized mixtures of normal distributions using Bayesian methods.
#
if (!require("mixAK")){install.packages("mixAK")}
library(mixAK)
#
par(mfrow = c(2,2))
cols <- c("green", "orange", "blue", "red")
dic <- rep(NA, 4)
nMCMC <- c(burn=5000, keep=10000, thin=5, info=1000)
for (k in 1:4){
  kprior.fixed <- list(priorK = "fixed", Kmax = k)
  fit.bayes <- NMixMCMC(y0 = y, prior = kprior.fixed,
                      scale = list(shift = 0, scale = 1), nMCMC = nMCMC, 
                      PED = F)
  pdens <- NMixPredDensMarg(fit.bayes, lgrid=300)
  hist(y, freq = FALSE, xlab = "response times", ylab = "f", 
       main = paste(k, "component fit"))
  lines(pdens$x$x1, pdens$dens$"1", lwd = 2, col = cols[k])
  dic[k] <- fit.bayes$DIC[1]
}
#
# The selected single model is the 3 component mixture.
#
par(mfrow = c(1,1))
#
# Use reversible jump to run with a prior on k.
#
kprior.unif <- list(priorK = "uniform", Kmax=30)
RJModel <- NMixMCMC(y0 = y, prior = kprior.unif, nMCMC = nMCMC, 
                    scale=list(shift = 0, scale = 1),PED = F)
pdensRJ <- NMixPredDensMarg(RJModel, lgrid = 300)
#
# Now replot the histogram with the fitted density.
#
hist(y, freq = FALSE, xlab = "response times", 
     ylab = "f", main = "Response time data with fitted density", 
     ylim = c(0, max(pdensRJ$dens$"1")))
lines(pdensRJ$x$x1, pdensRJ$dens$"1", lwd = 2, col = "blue")
#
# Look at the mixing of k and the posterior probability distribution of k.
#
ts.plot(RJModel$K, xlab = "iteration", ylab = "k")
barplot(prop.table(table(RJModel$K)), xlab = "k", ylab = "P(k|data)", 
        ylim = c(0,1))
#
# Look at the series over time.
#
plot(y, type = "l", xlab = "test number", ylab = "response time")
acf(y, xlab = "lag", ylab = "ACF")
#
# This looks like a time series!
#

















