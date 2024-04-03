rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graph")
BiocManager::install("Rgraphviz")
if (!require("bnclassify"))
  install.packages("bnclassify")
library(graph)
library(Rgraphviz)
library(bnclassify)
data(car)
head(car)
#
# Use a naive Bayes model to predict the acceptability of a car according to it's features.
# First set up and plot the model.
#
nbfit <- nb("class",car)
plot(nbfit)
#
# Use MLEs for parameter estimation.
#
nbfit <- lp(nbfit,car,smooth=0)
nbfit$.params
#
# Now do with simple Laplacian smoothing.
#
nbfit <- lp(nbfit,car,smooth=1)
nbfit$.params
#
# Now look at (in sample prediction).  We'll make a table to show the correct and misclassified cars.
#
carpred <- predict(nbfit,car)
successtab <- table(carpred,true=car$class)
successtab
#
# Set up a more complex model. First try using AIC.
#
ode_cl_aic <- tan_cl('class', car, score = 'aic')
plot(ode_cl_aic)
aicfit <- lp(ode_cl_aic,car,smooth=1)
carpred_aic <- predict(aicfit,car)
successtab_aic <- table(carpred_aic,true=car$class)
#
# Now try with BIC.  We get the Naive Bayes model
#
ode_cl_bic <- tan_cl('class', car, score = 'bic')
plot(ode_cl_bic)
#
# Now try with a semi-naive Bayes forward selection method
#
fssj <- fssj('class', car, k = 5, epsilon = 0)
plot(fssj)
fsfit <- lp(fssj,car,smooth=1)
carpred_fs <- predict(fsfit,car)
successtab_fs <- table(carpred_fs,true=car$class)
successtab_fs
