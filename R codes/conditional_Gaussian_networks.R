rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graph")
BiocManager::install("Rgraphviz")
if (!require("rbmn"))
  install.packages("rbmn")
if (!require("bnlearn"))
  install.packages("bnlearn")
if (!require("DAAG"))
  install.packages("DAAG")
if (!require("ggplot2"))
  install.packages("ggplot2")
#
# Data
#
library(DAAG)
data(ais)
help(ais)
head(ais)
#
# Plot data:  look at hemoglobin levels vs sport.  There is a dependence on sport type ...
#
library(graph)
library(ggplot2)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
ggplot(ais, aes(x = sport, y = hg, fill = sport)) + geom_boxplot() + 
  scale_fill_manual(values = colorRampPalette(YlOrBr)(10))
#
# Also plot versus gender: again we see a dependence
#
ggplot(ais, aes(x = sex, y = hg, fill = sex)) + geom_boxplot() 
#
# Compare hemoglobin levels with hematocrit levels.  There is a clear, linear dependence.
#
plot(ais$hc,ais$hg,xlab="hc",ylab="hg")
#
# Also compare hemoglobin levels with lean body mass levels.  There is a slight, positive relation.
#
plot(ais$lbm,ais$hg,xlab="lbm",ylab="hg")
#
# Just look at a subset of the data.
#
ais.sub <- ais[c("hc", "hg", "sport", "lbm","sex")]
#
# Try to learn the dependence structure.
#
structure <- hc(ais.sub, score = "bic-cg")
plot(structure)
#
# It doesn't make sense that "sex"" depends on "sport".  We should change stop
# this using a blacklist.
#
b <- c("sport","sex")
structure <- hc(ais.sub, score = "bic-cg",blacklist=b)
plot(structure)
#
# Now let's fit the network.
#
bn <- bn.fit(structure,ais.sub)
bn
#
# Check that what is happening for continuous nodes is regression
#
bn$hg
lm(hg ~ hc + lbm,data=ais.sub)
bn$hc
lm(hc ~ sex,data=ais.sub) # summing intercept + sexm gives the value when sex = 1.
#
# Look at sport: we can see the conditioning on sex.
#
bn$sport
#
# Now make a query.
#
cpquery(bn,event =(sex == "f"),evidence = list(sport = "B_Ball"), method = "lw")
cpquery(bn, event = (lbm > 60 & sex == "f"), evidence = list(sport = "B_Ball"), 
        method = "lw")