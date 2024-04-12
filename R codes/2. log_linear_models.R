rm(list = ls())
#
# Install required packages.
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("graph")
BiocManager::install("Rgraphviz")
#
# When we install gRim, we also install the package gRbase.
#
if (!require("gRim"))
  install.packages("gRim")
library(graph)
library(gRim)
#
# A simple 2 way contingency table example.
#
freq <- matrix(c(4,10,11,10,8,14),nrow=3)
rownames(freq) <- c("Humanities","Natural Sciences","Social Sciences")
colnames(freq) <- c("Male","Female")
freq <- as.table(freq)
names(attributes(freq)$dimnames) <- c("Major","Gender")
addmargins(freq)
#
# Calculate MLEs for the cell probabilities without assuming independence.
#
prop.table(freq)
#
# Calculate MLEs for the cell probabilities when we assume independence.
#
X <- chisq.test(freq)
X$expected/sum(freq)
#
# Perform a chi-squared test for independence. 
#
freq # are the observed values,
X$expected # are the expected values,
X$statistic # is the test statistic,
X$p.value # is the p-value
#
# Calculation of the deviance.
#
dmod(~Major:Gender,data=freq)  #The saturated model
dmod(~Major+Gender,data=freq)  #The independent model (the plus makes it so we are treating both things as independent)
pchisq(2.38, df = 2, lower.tail = FALSE) # 0.3042213, pretty similar to the x$p.value 0.3120586 computed
# deviance = 2.38, pretty similar to the x$statistic 2.329129 computed
#
# Fit the Poisson regression model.
#
glm(Freq~Major+Gender,family=poisson,data=as.data.frame(freq))
# The residual deviance is our x$statistics from before or the deviance from dmod
#
# 3 way table example.
#
counts <- c(911,3,44,2,538,43,456,279)
dimn <- list(alcohol=c("Yes","No"),tobacco=c("Yes","No"),marijuana=c("Yes","No"))
drug_use <- as.table(array(counts,dim=c(2,2,2),dimnames=dimn)) #Contingency table format
as.data.frame(drug_use) #Original format
xtabs(Freq~.,data=as.data.frame(drug_use)) #Gets back to contingency table format again
#
# Dependence graph of a log-linear model.  Use either one of the following two commands.
#
m1 <- dmod(~alcohol*tobacco+tobacco*marijuana,
           data=drug_use) #2 way interactions
plot(m1)
summary(m1)
m1 <- dmod(list(c("alcohol","tobacco"),
                c("tobacco","marijuana")),data=drug_use) #The same
plot(m1)
summary(m1)
#
# Compare the dependence graphs of the full model and the model with only 2 way interactions.
#
m2 <- dmod(~alcohol*tobacco+tobacco*marijuana+marijuana*alcohol,data=drug_use) #2 way interactions 
plot(m2)
m3 <- dmod(~alcohol*tobacco*marijuana,data=drug_use) #Full model
plot(m3)
# We have the same graph but a different model, in the second one we have the term XYZ, while in the first we had
# the 3 pairs
#
# Perform maximum likelihood estimation for the 2 way interactions model.
#
c <- loglin(drug_use,margin=list(c(1,2),c(1,3),c(2,3)),fit=TRUE)
c$fit
#
# Check if different models are graphical.
#
isGraphical(~alcohol*tobacco+tobacco*marijuana) # m1
isGraphical(~alcohol*tobacco+tobacco*marijuana+marijuana*alcohol) # m2
isGraphical(~alcohol*tobacco*marijuana) # m3
#
# Check if the models are decomposable.
#
isDecomposable(~alcohol*tobacco+tobacco*marijuana)
isDecomposable(~alcohol*tobacco+tobacco*marijuana+marijuana*alcohol)
#
# Another way of checking for graphical and decomposable models
#
m1$modelinfo$properties
m2$modelinfo$properties
m3$modelinfo$properties
#
# Check RIP ordering.
#
ug1 <- ug(~alcohol*tobacco+tobacco*marijuana)
rip(ug1)
#
# Fit with ips.  Note that this model is not decomposable.
#
cc <- loglin(drug_use,list(c(1,2),c(1,3),c(2,3)),fit = TRUE)
cc$fit
#
# Calculate AIC and BIC for different models
#
dmod(~.^.,data=drug_use)
dmod(~alcohol*tobacco+tobacco*marijuana+marijuana*alcohol,data=drug_use)
dmod(~alcohol*tobacco+tobacco*marijuana,data=drug_use)
dmod(~alcohol*tobacco+marijuana*alcohol,data=drug_use)
dmod(~tobacco*marijuana+marijuana*alcohol,data=drug_use)
dmod(~alcohol+tobacco*marijuana,data=drug_use)
dmod(~tobacco+alcohol*marijuana,data=drug_use)
dmod(~marijuana+alcohol*tobacco,data=drug_use)
dmod(~alcohol+tobacco+marijuana,data=drug_use)
#
# Summarize the characteristics of the best fitting model and the best fitting graphical model
#
summary(dmod(~alcohol*tobacco+tobacco*marijuana,data=drug_use))
summary(dmod(~.^.,data=drug_use))
#
# Searching for the best decomposable model in a high dimensional table: using AIC
#
data(reinis)
msat <- dmod(~.^.,data=reinis)
mbaic <- backward(msat)
plot(mbaic)
mmain <- dmod(~.^1,data=reinis)
mfaic <- forward(mmain)
plot(mfaic)
#
# ... using BIC  
#
msat <- dmod(~.^.,data=reinis)
mbbic <- backward(msat,k=log(sum(reinis)))
plot(mbbic)
mmain <- dmod(~.^1,data=reinis)
mfbic <- forward(mmain,k=log(sum(reinis)))
plot(mfbic)
#
# We can also search without the restriction to decomposable models.  The fitted 
# decomposable (lhs) and unrestricted (rhs) models do not always have the same graphs.
#
par(mfrow=c(1,2))
mbaic_unrestricted <- backward(msat,type="unrestricted")
plot(mbaic)
plot(mbaic_unrestricted)
mfaic_unrestricted <- forward(mmain,type="unrestricted")
plot(mfaic)
plot(mfaic_unrestricted)
mbbic_unrestricted <- backward(msat,type="unrestricted",k=log(sum(reinis)))
plot(mbbic)
plot(mbbic_unrestricted)
mfbic_unrestricted <- forward(mmain,type="unrestricted",k=log(sum(reinis)))
plot(mfbic)
plot(mfbic_unrestricted)
#
# Note that the fitted model in the unrestricted case is not necessarily decomposable.
#
mbaic$modelinfo$properties
mbaic_unrestricted$modelinfo$properties
mfaic$modelinfo$properties
mfaic_unrestricted$modelinfo$properties
mbbic$modelinfo$properties
mbbic_unrestricted$modelinfo$properties
mfbic$modelinfo$properties
mfbic_unrestricted$modelinfo$properties
#
# Goodness of fit checks: none of the fitted models are rejected by a chisq test
#
pchisq(mbaic$fitinfo$dev,mbaic$fitinfo$dimension[4])
pchisq(mfaic$fitinfo$dev,mfaic$fitinfo$dimension[4])
pchisq(mbbic$fitinfo$dev,mbbic$fitinfo$dimension[4])
pchisq(mfbic$fitinfo$dev,mfbic$fitinfo$dimension[4]) 
pchisq(mbaic_unrestricted$fitinfo$dev,mbaic$fitinfo$dimension[4])
pchisq(mfaic_unrestricted$fitinfo$dev,mfaic$fitinfo$dimension[4])
pchisq(mbbic_unrestricted$fitinfo$dev,mbbic$fitinfo$dimension[4])
pchisq(mfbic_unrestricted$fitinfo$dev,mfbic$fitinfo$dimension[4])
#
# Which of all the models should be selected according to aic?
#
aic <- rep(NA,8)
aic[1] <- mbaic$fitinfo$aic
aic[2] <- mbaic_unrestricted$fitinfo$aic
aic[3] <- mfaic$fitinfo$aic
aic[4] <- mfaic_unrestricted$fitinfo$aic
aic[5] <- mbbic$fitinfo$aic
aic[6] <- mbbic_unrestricted$fitinfo$aic
aic[7] <- mfbic$fitinfo$aic
aic[8] <- mfbic_unrestricted$fitinfo$aic
aic
#
# The unrestricted model generated by backward and forward search are the same.  
#
# Now select according to BIC.
#
bic <- rep(NA,8)
bic[1] <- mbaic$fitinfo$bic
bic[2] <- mbaic_unrestricted$fitinfo$bic
bic[3] <- mfaic$fitinfo$bic
bic[4] <- mfaic_unrestricted$fitinfo$bic
bic[5] <- mbbic$fitinfo$bic
bic[6] <- mbbic_unrestricted$fitinfo$bic
bic[7] <- mfbic$fitinfo$bic
bic[8] <- mfbic_unrestricted$fitinfo$bic
bic
#
# Plot the graphs of the two optimal (non-decomposable) models.
#
plot(mbaic_unrestricted)
plot(mbbic_unrestricted)

