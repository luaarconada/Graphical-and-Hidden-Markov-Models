#
# Gaussian network example
#
rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graph")
BiocManager::install("Rgraphviz")
if (!requireNamespace("bnlearn", quietly = TRUE))
  install.packages("bnlearn")
library(graph)
library(bnlearn)
#
# The data are marks on various different exams.
#
help(marks)
data(marks)
#
# The following is a possible suggestion for the dependence structure.
#
marks.dag = model2network("[ALG][ANL|ALG][MECH|ALG:VECT][STAT|ALG:ANL][VECT|ALG]")
plot(marks.dag)
#
# We can fit the model using linear regression.
#
gbn <- bn.fit(marks.dag,marks)
gbn$MECH
#
# Check that what we are doing at each node really is regression.
#
regr.res <- lm(MECH ~ ALG+VECT, data = marks)
regr.res
#
# Look at the complete joint normal distribution.
#
gbn2mvnorm(gbn)
#
# Make a query given certain evidence.
#
cpquery(gbn, event = (ALG > 50), evidence = list(STAT = 50), 
              method = "lw")
cpquery(gbn,event = (ALG > 50  & STAT > 40), evidence = 
                list(MECH = c(40,60),STAT = 50), method = "lw")
#
# Learn the network structure.  (There are many possible algorithms)
# This is a hill climbing algorithm.
#
plot(hc(marks))
#
# Use a blacklist to prohibit odd dependencies or ...
#
(b <- data.frame(from = c("VECT","MECH"), to = c("ALG","ALG")))
plot(hc(marks,blacklist = b))
#
# ... use a whitelist to force certain dependencies.
#
(w <- data.frame(from = "ALG", to = "MECH"))
plot(hc(marks,whitelist = w))

