library(markovchain)
trans_matrix <- matrix(c(0,0.5,0.25,0.25,0.5,0,0.5,0,0.5,0.5,0,0,0,0,0,1),
                       nrow=4,byrow=TRUE)
trans_matrix
sillymc <- new("markovchain",transitionMatrix=trans_matrix)
#
# Look at state probabilities 1 and 2 steps ahead and more
#
init <- c(0,0,1,0)
fin <- init*sillymc
fin
steps <- 2
fin <- init*sillymc^steps
fin
steps <- 20
fin <- init*sillymc^steps
fin
steps <- 200
fin <- init*sillymc^steps
fin
#
# Plot of the MC structure
#
plot(sillymc)
summary(sillymc)
#
# Check the types of states in the chain.
#
transientStates(sillymc)
recurrentStates(sillymc)
absorbingStates(sillymc)
#
# Now define an ergodic chain.
#
trans_matrix <- matrix(c(0.5,0.5,0,0.25,0.25,0.5,0.2,0.4,0.4),nrow=3,byrow=TRUE)
trans_matrix
ergodicmc <- new("markovchain",transitionMatrix=trans_matrix)
#
# Check the types of states
#
transientStates(ergodicmc)
recurrentStates(ergodicmc)
absorbingStates(ergodicmc)
#
# Check the period
#
period(ergodicmc)
#
# Calculate the state probabilities various steps ahead
#
init <- c(0,1,0)
init*ergodicmc
steps <- 2
init*ergodicmc^steps
steps <- 20
init*ergodicmc^steps
steps <- 200
init*ergodicmc^steps
#
# Calculate the steady state distribution.
#
steadyStates(ergodicmc)
#
# Real data example of inference for a MC.
#
data(rain)
head(rain)
help(rain)
raindata <- rain$rain
#
# Calculation matrix showing the numbers of state transitions
#
transition_numbers <- createSequenceMatrix(raindata)
transition_numbers
rainmc <- markovchainFit(data=raindata)
#
# Estimate the probability transition matrix
#
est_trans <- rainmc$estimate
est_trans
plot(est_trans)
#
# One step ahead forecast
#
init <- c(0,1,0)
onestep <- init%*%est_trans@transitionMatrix
onestep
#
# Steady state probability.
#
steadyStates(est_trans)
#
# Bayesian estimation.
#
bayes_fit <- markovchainFit(data=raindata,method="laplace",laplacian=1)
bayes_fit$estimate
n <- 1000
ss <- matrix(rep(NA,n*3),nrow=n)
dirichlet_matrix <- transition_numbers+1
tmean <- matrix(rep(0,9),nrow=3)
for (i in 1:n){
  transi <- t(matrix(rgamma(9,dirichlet_matrix,1),ncol=3))
  transi <- transi/apply(transi,1,sum)
  tmean <- tmean+transi
  mc <- new("markovchain",transitionMatrix = transi)
  ss[i,] <- steadyStates(mc)
}
#
# Calculate the posterior means of the transition matrix and steady state 
# probabilities.
#
tmean <- tmean/n
tmean
apply(ss,2,mean)
