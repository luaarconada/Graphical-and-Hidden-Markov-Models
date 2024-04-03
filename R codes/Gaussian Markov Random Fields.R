#
# Gaussian MRF example.
#
rm(list = ls())
#
# Install required packages.
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graph", version = "3.8")
BiocManager::install("Rgraphviz", version = "3.8")
if (!require("gRbase"))
  install.packages("gRbase")
if (!require("rbmn"))
  install.packages("rbmn")
if (!require("bnlearn"))
  install.packages("bnlearn")
if (!require("gRim"))
  install.packages("gRim")
if (!require("glasso"))
  install.packages("glasso")
library(graph)
library(bnlearn)
library(rbmn)
library(gRbase)
library(gRim)
library(glasso)
#
# Example of not being able to see conditional independences clearly from a covariance matrix
#
Q <- matrix(c(1,0,0.2,0,1,0.5,0.2,0.5,1),nrow=3)
Sigma <- solve(Q)
Sigma <- rbind(Sigma,c(rep(0,3)))
Sigma <- cbind(Sigma,c(rep(0,3),2.124))
rownames(Sigma) <- c("X1","X2","X3","X4")
colnames(Sigma) <- c("X1","X2","X3","X4")
round(Sigma,3)
Q <- solve(Sigma)
round(Q,3)
#
# Load up a body characteristics data set
#
help(boco)
data(boco)
head(boco)
#
# Estimate the mean.
#
apply(boco,2,mean)
#
# Estimate the covariance and precision matrices
#
covboco <- cov.wt(boco,method="ML")$cov
round(covboco,digits=2)
precboco <- solve(covboco)
round(precboco,digits=2)
#
# Calculate partial correlations
#
parcorrboco <- cov2pcor(covboco)
round(100*parcorrboco)
#
# Calculate the saturated model
#
satmod <- cmod(~.^.,data=boco)
#
# Set up an alternative model based on removing the terms with partial correlation <= 0.05
#
mod1 <- update(satmod,list(dedge= ~A:LF+A:TB+H:LF+H:LL+W:C+C:LL+C:LB+TF:TB+LF:LB+AF:LB+TL:TB+LL:LB))
plot(mod1)
prec1 <- mod1$fitinfo$K
cov1 <- solve(prec1)
parcorr1 <- cov2pcor(cov1)
round(prec1,digits=2)
round(cov1,digits=2)
round(100*parcorr1)
#
# Test the goodness of fit
#
satmod
mod1
pchisq(mod1$fitinfo$dev,mod1$fitinfo$dimension[4])
#
# Select an optimum model via stepwise regression using AIC ...
#
maic <- stepwise(satmod)
maic
plot(maic)
precaic <- maic$fitinfo$K
round(precaic,digits=2)
parcorraic <- cov2pcor(precaic)
round(100*parcorraic)
#
# Check this fits the data by testing the deviance.
# This fits.
#
pchisq(maic$fitinfo$dev,maic$fitinfo$dimension[4])
#
# ... using BIC
#
mbic <- stepwise(satmod,k=log(nrow(boco)))
mbic
lot(mbic)
precbic <- mbic$fitinfo$K
round(precbic,digits=2)
parcorrbic <- cov2pcor(precbic)
round(100*parcorrbic)
#
# This doesn't fit according to the deviance test.
#
pchisq(mbic$fitinfo$dev,mbic$fitinfo$dimension[4])
#
# Thresholding.
#
threshold <- 0.05
adjmat <- parcorrboco
adjmat <- abs(adjmat)
adjmat[adjmat<threshold] <- 0
diag(adjmat) <- 0
adjmat[adjmat>0] <- 1
g.thresh <- as(adjmat,"graphNEL")
mthresh <- cmod(g.thresh,data=boco)
plot(mthresh)
pchisq(mthresh$fitinfo$dev,mthresh$fitinfo$dimension[4])
#
# Glasso
#
corrboco <- cov2cor(covboco)
reslasso <- glasso(corrboco,rho=0.1)
AM <- reslasso$wi != 0
diag(AM) <- F
g.lasso <- as(AM,"graphNEL")
nodes(g.lasso) <- names(boco)
mlasso <- cmod(g.lasso,data=boco)
plot(mlasso)
pchisq(mlasso$fitinfo$dev,mlasso$fitinfo$dimension[4])