#
# Bayesian nets.
#
rm(list=ls())
#
# Packages.
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graph")
BiocManager::install("Rgraphviz")
if (!require("gRain"))
  install.packages("gRrain")
if (!require("bnlearn"))
  install.packages("bnlearn")
library(graph)
library(gRain)
library(bnlearn)
library(Rgraphviz)
#
# Data
#
data(asia)
dim(asia)
head(asia)
#
# Graphical representation
#
dag <- model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
graphviz.plot(dag)
# plot(dag) could also be used.
#
# Defining probabilities via expert judgements
#
yn <- c("yes", "no")
cptA <- matrix(c(0.01, 0.99), ncol=2, dimnames=list(NULL, yn))
cptS <- matrix(c(0.5, 0.5), ncol=2, dimnames=list(NULL, yn))
cptT <- matrix(c(0.05, 0.95, 0.01, 0.99), ncol=2, dimnames=list("T"=yn, "A"=yn))
cptL <- matrix(c(0.1, 0.9, 0.01, 0.99), ncol=2, dimnames=list("L"=yn, "S"=yn))
cptB <- matrix(c(0.6, 0.4, 0.3, 0.7), ncol=2, dimnames=list("B"=yn, "S"=yn))
cptX <- matrix(c(0.98, 0.02, 0.05, 0.95), ncol=2, dimnames=list("X"=yn, "E"=yn))
# cptE and cptD are 3-d matrices, which don't exist in R, so
# we need to build these manually as below.
cptE <- c(1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0)
dim(cptE) <- c(2, 2, 2)
dimnames(cptE) <- list("E"=yn, "L"=yn, "T"=yn)
cptD <- c(0.9, 0.1, 0.7, 0.3, 0.8, 0.2, 0.1, 0.9)
dim(cptD) <- c(2, 2, 2)
dimnames(cptD) <- list("D"=yn, "E"=yn, "B"=yn)
expertfit <- custom.fit(dag, dist=list(A=cptA, S=cptS, T=cptT, L=cptL,
                                       B=cptB, E=cptE, X=cptX, D=cptD))
expertfit
#
# Using maximum likelihood.
#
mlefit = bn.fit(dag,asia)
mlefit
mlefit$E
#
# Bayesian approach / Laplacian smoothing
#
bayesfit <- bn.fit(dag,asia,method="bayes")
bayesfit
bayesfit$E
#
# Prediction
#
v <- mlefit$T  # Try this with bayesfit or expertfit instead.
pTgivenA <- v$prob
v <- mlefit$A
pA <- v$prob
pT <- unname(pTgivenA[2,1]*pA[1]+pTgivenA[2,2]*pA[2])
pT
#
# Exact prediction: moralization.
#
mlefit_gr <- as.grain(mlefit)
mg <- moralize(mlefit_gr$dag)
plot(mg)
#
# Triangulation.
#
plot(triangulate(mg))
#
# Prediction
#
mlefit_gr <- compile(mlefit_gr)
mlefit_gr <- propagate(mlefit_gr)
querygrain(mlefit_gr,nodes="T",type="marginal")
#
# A more complicated prediction.
#
mlefit_gr_ev <- setFinding(mlefit_gr,nodes=c("D","X"),states=c("yes","yes"))
querygrain(mlefit_gr_ev,nodes="L",type="marginal")
#
# Approximate prediction.
#
cpquery(mlefit, (L=="yes"), (D=="yes" & X=="yes"),n=100)
cpquery(mlefit, (L=="yes"), (D=="yes" & X=="yes"),n=1000)
cpquery(mlefit, (L=="yes"), (D=="yes" & X=="yes"),n=10000)
cpquery(mlefit, (L=="yes"), (D=="yes" & X=="yes"),n=100000)
cpquery(mlefit, (L=="yes"), (D=="yes" & X=="yes"),n=1000000)
#
# Likelihood weighting.
# 
cpquery(mlefit, (L=="yes"), evidence=list(D="yes",X="yes"),method="lw")
cpquery(mlefit, (L=="yes"), evidence=list(D="yes",X="yes"),method="lw",n=10000)
cpquery(mlefit, (L=="yes"), evidence=list(D="yes",X="yes"),method="lw",n=100000)
cpquery(mlefit, (L=="yes"), evidence=list(D="yes",X="yes"),method="lw",n=1000000)
#
# Learning the independence structure.
#
# Chi-sq test.  Are D and X independent?
#
ci.test("D","X",data=asia)
#
# Are D and X independent given E?
#
ci.test("D","X","E",data=asia)
#
# Different structure learning algorithms.
#
help("structure-learning")
#
# Algorithms based on the Markov blanket.
#
mb_dag <- iamb(asia)
mb_dag
plot(mb_dag)
#
# Blacklists to ban silly dependencies ...
#
bl <- data.frame(from = "D",to = "B")
mb_dag_bl <- iamb(asia, blacklist = bl)
mb_dag_bl
plot(mb_dag_bl)
#
# ... and whitelists to force obvious dependencies ...
#
wl <- data.frame(from = c("S","E"), to = c("B","D"))
wl
mb_dag_wl <- iamb(asia,whitelist = wl)
mb_dag_wl
plot(mb_dag_wl)
#
# ... or both.
#
mb_dag_blwl <- iamb(asia, blacklist = bl, whitelist = wl)
plot(mb_dag_wl)
#
# Score based methods.
#
score_dag <- hc(asia)
score_dag
plot(score_dag)
#
# Hybrid algorithm.
#
hybrid_dag <- mmhc(asia)
hybrid_dag
plot(hybrid_dag)
#
# Sometimes these algorithms cannot give a full ordering ...
# 
plot(pc.stable(asia))
