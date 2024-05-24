rm(list = ls())
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
# ------------------------------------------------------------------------------
#
# MIXTURE MODELS  
#
# Draw a bar chart of the data.
#
counts <- table(y)
plot(sort(unique(y)), counts/sum(counts), type = 'h', xlim = c(0,45), 
     ylim = c(0,0.12), xlab = "number of earthquakes", ylab = "frequency", 
     main = "Histogram of yearly earthquake numbers between 1900 and 2006")
ytick <- seq(0, 0.12, by = 0.02)
axis(side = 2, at = ytick, labels = ytick)
#
# Fit a Poisson distribution to the data.
#
grid.y <- seq(min(y), max(y), 1)
mean.y <- mean(y)
lines(grid.y, dpois(grid.y, mean.y), type = 'h', lwd = 2, lty = 2, col = "blue")
#
# Check that this does not fit.
#
if (!require("energy")){install.packages("energy")}
energy::poisson.mtest(as.integer(y), R = 199)
#
# Now fit mixture models of different sizes.
#
aic <- rep(NA, 5)
bic <- rep(NA, 5)
dev <- rep(NA, 5)
dev[1] <- -2*sum(dpois(y, mean.y,log = TRUE))
aic[1] <- dev[1]+2
bic[1] <- dev[1]+2*log(n)
if (!require("depmixS4")){install.packages("depmixS4")}
library(depmixS4)
earthquake <- data.frame(x = x, y = y)
#
# 2 components.
#
temp <- mix(y ~ 1, nstates = 2, family=poisson(), data=earthquake)
fit2 <- fit(temp)
params <- getpars(fit2)
p2 <- unname(params[1:2])
lambda2 <- unname(exp(params[3:4]))
aic[2] <- AIC(fit2)
bic[2] <- BIC(fit2)
#
# Plot the weighted components and the full probability distribution fit.
#
par(mfrow = c(2,2))
temp <- c(-0.2, 0.2)
cols <- c("green", "yellow")
dens2 <- rep(0, length(grid.y))
plot(sort(unique(y)), counts/sum(counts), type = 'h', xlim = c(0,45), 
     ylim = c(0,0.12), xlab = "number of earthquakes", ylab = "frequency",
     main = "mixture size 2")
ytick <- seq(0, 0.12, by = 0.02)
axis(side = 2, at = ytick, labels = ytick)
for (i in 1:2){
  dens <- p2[i]*dpois(grid.y, lambda2[i])
  lines(grid.y+temp[i], dens, type = "h", lty = 1, lwd = 2, col = cols[i])
  dens2 <- dens2+dens
}
lines(grid.y, dens2, type = "h", lty = 3, lwd = 2, col = "blue")
#
# 3 components.  
#
temp <- mix(y ~ 1, nstates = 3, family=poisson(), data=earthquake)
fit3 <- fit(temp)
params <- getpars(fit3)
p3 <- unname(params[1:3])
lambda3 <- unname(exp(params[4:6]))
aic[3] <- AIC(fit3)
bic[3] <- BIC(fit3)
temp <- c(-0.2, 0, 0.2)
cols <- c("green", "red", "yellow")
dens3 <- rep(0, length(grid.y))
plot(sort(unique(y)), counts/sum(counts), type = 'h', xlim = c(0,45),
     ylim = c(0,0.12), xlab = "number of earthquakes", ylab = "frequency",
     main = "mixture size 3")
ytick <- seq(0, 0.12, by = 0.02)
axis(side = 2, at = ytick, labels = ytick)
for (i in 1:3){
  dens <- p3[i]*dpois(grid.y, lambda3[i])
  lines(grid.y+temp[i], dens, type = "h", lty = 1, lwd = 2, col = cols[i])
  dens3 <- dens3+dens
}
lines(grid.y, dens3, type = "h", lty = 3, lwd = 2, col = "blue")
#
# 4 components.  
#
temp <- mix(y ~ 1, nstates = 4, family = poisson(), data = earthquake)
fit4 <- fit(temp)
params <- getpars(fit4)
p4 <- unname(params[1:4])
lambda4 <- unname(exp(params[5:8]))
aic[4] <- AIC(fit4)
bic[4] <- BIC(fit4)
temp <- c(-0.4,-0.2,0.2,0.4)
cols <- c("green", "red", "yellow", "purple")
dens4 <- rep(0, length(grid.y))
plot(sort(unique(y)), counts/sum(counts), type = 'h', xlim = c(0,45), 
     ylim = c(0,0.12), xlab = "number of earthquakes", ylab = "frequency",
     main = "mixture size 4")
ytick <- seq(0, 0.12, by = 0.02)
axis(side = 2, at = ytick, labels = ytick)
for (i in 1:4){
  dens <- p4[i]*dpois(grid.y, lambda4[i])
  lines(grid.y+temp[i], dens, type = "h", lty = 1, lwd = 2, col = cols[i])
  dens4 <- dens4+dens
}
lines(grid.y, dens4, type = "h", lty = 3, lwd = 2, col = "blue")
#
# 5 components.  
#
temp <- mix(y ~ 1, nstates = 5, family = poisson(), data = earthquake)
fit5 <- fit(temp)
params <- getpars(fit5)
p5 <- unname(params[1:5])
lambda5 <- unname(exp(params[6:10]))
aic[5] <- AIC(fit5)
bic[5] <- BIC(fit5)
temp <- c(-0.4, -0.2, 0, 0.2, 0.4)
cols <- c("green", "red", "yellow", "purple", "orange")
dens5 <- rep(0, length(grid.y))
plot(sort(unique(y)), counts/sum(counts), type = 'h', xlim = c(0,45), 
     ylim = c(0,0.12), xlab = "number of earthquakes", ylab = "frequency",
     main = "mixture size 5")
ytick <- seq(0, 0.12, by = 0.02)
axis(side = 2, at = ytick, labels = ytick)
for (i in 1:5){
  dens <- p5[i]*dpois(grid.y, lambda5[i])
  lines(grid.y+temp[i], dens, type = "h", lty = 1, lwd = 2, col = cols[i])
  dens5 <- dens5+dens
}
lines(grid.y, dens5, type = "h", lty = 3, lwd = 2, col = "blue")
#
# Compare the AIC and BIC
#
aic
bic
which.min(aic)
which.min(bic)
#
# AIC suggests a 3 component mixture and BIC a 2 component mixture.
#
# ------------------------------------------------------------------------------
#
# Look at the time series of the data
#
par(mfrow = c(1,1))
plot(x, y, type = "l", ylab = "earthquakes", 
     main = "Number of yearly earthquakes from 1900 to 2006")
#
# Check the acf.
#
acf(y, main = "ACF of yearly earthquakes from 1900 to 2006")
#
# Fit a hidden Markov model.
#
nstates <- 3
hmminit <- depmix(y ~ 1, nstates = nstates, family = poisson(), 
                  data = earthquake)
set.seed(1)
hmmfit <- fit(hmminit)
#
# Extract parameters
#
params <- getpars(hmmfit)
(pinit <- params[1:nstates])
(trans.matrix <- t(matrix(params[(nstates+1):(nstates+nstates^2)], 
                         nrow = nstates)))
(lambda <- unname(exp(params[(nstates+nstates^2+1):(2*nstates+nstates^2)])))
temp <- depmixS4::posterior(hmmfit)
state.probs <- temp[ , 2:(nstates+1)]
head(state.probs)
(states <- temp[ , 1])
#
# Plot the data and add the most likely lambda value (associated with the most
# probable state) for each year. 
#
plot(x, y, type = 'l', lwd = 2, xlab = "Year", ylab = "Number of earthquakes", 
     ylim = c(0,50), main = "State changes")
lines(x, lambda[states], type = 'l', lwd = 2, col = 'green')
#
# Compare the estimated lambda values with those for the mixture model.
#
rates <- rbind(lambda3, lambda)
row.names(rates) <- c("lambda mix", "lambda hmm")
rates
#
# Now estimate the steady state distribution of the MC.  
#
library(markovchain)
ergodicmc <- new("markovchain", transitionMatrix = trans.matrix)
p.hmm <- steadyStates(ergodicmc)
steady.state.probs <- rbind(p3, p.hmm)
row.names(steady.state.probs) <- c("mixture", "hmm")
steady.state.probs
#
# The states may be named differently in the two models but the results are 
# clearly similar.
#
# Show the weighted components and marginal predictive distribution for both the
# mixture model and the HMM.
# 
par(mfrow = c(1, 2))
temp <- c(-0.2, 0, 0.2)
cols <- c("green", "red", "yellow")
dens3 <- rep(0, length(grid.y))
plot(sort(unique(y)), counts/sum(counts), type = 'h', xlim = c(0,45),
     ylim = c(0,0.12), xlab = "number of earthquakes", ylab = "frequency", 
     main = "Marginal density: mixture")
ytick <- seq(0, 0.12, by = 0.02)
axis(side = 2, at = ytick, labels = ytick)
for (i in 1:3){
  dens <- p3[i]*dpois(grid.y, lambda3[i])
  lines(grid.y+temp[i], dens, type = "h", lty = 1, lwd = 2, col = cols[i])
  dens3 <- dens3+dens
}
lines(grid.y, dens3, type = "h", lty = 3, lwd = 2, col = "blue")
#
dens.hmm <- rep(0, length(grid.y))
plot(sort(unique(y)), counts/sum(counts), type = 'h', xlim = c(0,45),
     ylim = c(0,0.12), xlab = "number of earthquakes", ylab = "frequency", 
     main = "Marginal density: hmm")
ytick <- seq(0, 0.12, by = 0.02)
axis(side = 2, at = ytick, labels = ytick)
for (i in 1:3){
  dens <- p.hmm[i]*dpois(grid.y, lambda[i])
  lines(grid.y+temp[i], dens, type = "h", lty = 1, lwd = 2, col = cols[i])
  dens.hmm <- dens.hmm+dens
}
lines(grid.y, dens.hmm, type = "h", lty = 3, lwd = 2, col = "blue")
par(mfrow = c(1, 1))
#
# Model selection.
#
aic.hmm <- rep(NA, 5)
bic.hmm <- rep(NA, 5)
aic.hmm[3] <- AIC(hmmfit)
bic.hmm[3] <- BIC(hmmfit)
nstates <- 2
hmminit <- depmix(y ~ 1, nstates = nstates, family = poisson(), 
                  data = earthquake)
temp <- fit(hmminit)
aic.hmm[2] <- AIC(temp)
bic.hmm[2] <- BIC(temp)
nstates <- 4
hmminit <- depmix(y ~ 1, nstates = nstates, family = poisson(), 
                  data = earthquake)
temp <- fit(hmminit)
aic.hmm[4] <- AIC(temp)
bic.hmm[4] <- BIC(temp)
nstates <- 5
hmminit <- depmix(y ~ 1, nstates = nstates, family = poisson(), 
                  data = earthquake)
temp <- fit(hmminit)
aic.hmm[5] <- AIC(temp)
bic.hmm[5] <- BIC(temp)
#
# The results are similar to the mixture model.  BIC selects a model with 2 
# states and AIC selects the 3 state model.
#