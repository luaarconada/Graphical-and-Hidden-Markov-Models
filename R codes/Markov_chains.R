rm(list=ls())
if (!require("markovchain"))
  install.packages("markovchain")
library(markovchain)
#
# General properties of Markov chains
#
(trans_matrix <- matrix(c(0,0.5,0.25,0.25,0.5,0,0.5,0,0.5,0.5,0,0,0,0,0,1),
                        nrow = 4, byrow = TRUE))
sillymc <- new("markovchain", transitionMatrix = trans_matrix)
#
# Check probability distribution of state in 2 steps given we start in state 3. 
# Try changing to 20 or 200.
#
init <- c(0,0,1,0)
steps <- 2
fin <- init*sillymc^steps
fin
#
# Plot the transition matrix structure.
#
plot(sillymc)
#
# Now look at an ergodic chain.
#
trans_matrix <- matrix(c(0.5, 0.5, 0, 0.25, 0.25, 0.5, 0.2, 0.4, 0.4),nrow=3,byrow=TRUE)
ergodicmc <- new("markovchain", transitionMatrix = trans_matrix)
init <- c(0,1,0)
steps <- 2
init*ergodicmc^steps
steps <- 20
init*ergodicmc^steps
steps <- 200
init*ergodicmc^steps
steadyStates(ergodicmc)
#
# Data 
#
data(rain)
head(rain)
help(rain)
#
# Calculate the numbers of transitions from each state.
#
raindata <-rain$rain
transition_numbers <- createSequenceMatrix(raindata)
transition_numbers
#
# Fit the Markov chain model.
#
rainmc <- markovchainFit(data=raindata)
#
# Look at the estimated transition matrix
#
est_trans <- rainmc$estimate
est_trans
plot(est_trans)
#
# One step ahead forecasting.
#
init <- c(0,1,0)
onestep <- init%*%est_trans@transitionMatrix
onestep
#
# Estimate the steady state distribution.
#
steadyStates(est_trans)
#
# Bayesian prediction
#
bayes_fit <- markovchainFit(data = raindata, method = "laplace", laplacian = 1)
bayes_fit$estimate # This gives a posterior mean estimate.
#
# Bayesian estimation via Monte Carlo
#
n <- 1000
ss <- matrix(rep(NA,n*3),nrow = n)
dirichlet_matrix <- transition_numbers+1
tmean <- matrix(rep(0,9), nrow = 3)
for (i in 1:n){
  transi <- t(matrix(rgamma(9, dirichlet_matrix, 1), ncol = 3))
  transi <- transi/apply(transi,1,sum)
  tmean <- tmean+transi
  mc <- new("markovchain", transitionMatrix = transi) 
  ss[i,] <- steadyStates(mc)
}
#
# Mean transition matrix. (Compare with the analytical estimate calculated 
# previously).
#
tmean <- tmean/n
tmean
#
# Mean steady state distribution
#
apply(ss,2,mean)
#
# Histograms of the densities of the steady state probabilities
#
par(mfrow=c(1,3))
hist(ss[,1], probability = TRUE, xlab = expression(pi[1]), main = "")
lines(density(ss[,1]))
hist(ss[,2],probability = TRUE, xlab = expression(pi[2]), main= "")
lines(density(ss[,2]))
hist(ss[,3],probability = TRUE, xlab = expression(pi[3]), main = "")
lines(density(ss[,3]))
par(mfrow = c(1,1))
