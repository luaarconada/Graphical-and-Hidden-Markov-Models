rm(list = ls())
#
# Install required packages.
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graph")
BiocManager::install("Rgraphviz")
if (!requireNamespace("rbmn", quietly = TRUE))
  install.packages("rbmn")
if (!requireNamespace("bnlearn", quietly = TRUE))
  install.packages("bnlearn")
library(graph)
library(bnlearn)
library(rbmn)
#
# Load up a body characteristics data set
#

help(boco)
data(boco)
head(boco)
#
# Learn the dependency structure of the Bayesian net.
#
iambex=iamb(boco)
iambex
plot(iambex)
hcex=hc(boco)
hcex
plot(hcex)
mmex <- mmhc(boco)
mmex
plot(mmex)
#
# Note that we could try to alter the implied dependence structure but here it is not obvious which way the dependencies go.
#
#
# Learn conditional probability structure.
#
gbn <- bn.fit(mmex, data = boco)
gbn
#
# See that what is being done is classical regression.
#
regr.res <- lm(LB ~ TB, data = boco)
regr.res
regr.res <- lm(AB ~ TB + LB, data = boco)
regr.res
#
# Look at the complete, joint normal distribution. (Uses routines from rbmn)
#
gbn.rbmn <- bnfit2nbn(gbn)
gema.rbmn <- nbn2gema(gbn.rbmn)
mn.rbmn <- gema2mn(gema.rbmn)
print8mn(mn.rbmn)
#
# Make a prediction.
#
cpquery(gbn,event=(C>100),evidence=list(A=53,W=95,H=178),method="lw")
cpquery(gbn,event=(C>100 & TF>15),evidence=list(A=53,W=95,H=178),method="lw")