rm(list=ls())
if(!require(bsts)){install.packages("bsts")}
library(bsts)     # load the bsts package
data(iclaims)     # bring the initial.claims data into scope
#
# Plot data
#
ts.plot(initial.claims$iclaimsNSA,lwd=2,ylab='claims')
#
# Introduce different model components
#
ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
#
# Fit model
#
model1 <- bsts(initial.claims$iclaimsNSA,
               state.specification = ss,
               niter = 1000)
plot(model1)
plot(model1, "components")  
#
# Out of sample prediction (for 2013 along with the previous 3 years)
# 
pred1 <- predict(model1, horizon = 12)
plot(pred1, plot.original = 156)