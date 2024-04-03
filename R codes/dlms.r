#
# Bayesian analysis of the Lake Huron data via a first order polynomial DLM.
#
# Load up data
#
rm(list=ls())
if(!require(dlm)){install.packages("dlm")}
library(dlm)
plot.ts(LakeHuron)
#
# MODEL 0: Simple linear regression vs time
#
t <- c(1875:1972)
model0 <- lm(LakeHuron ~ t)
f <- model0$fitted.values
lines(t,f,lwd=2,col='red')
r <- model0$residuals
plot(t,r,col='red')
lines(t,rep(0,length(t)))
acf(r,lwd=2,col='red',main='ACF of residuals')
# 
# MODEL 1: Kalman filter with filter estimates of variances.
#
mod1 <- dlmModPoly ( order=1 )
# Estimate the filtered values of the state vector 
filt1 <- dlmFilter ( LakeHuron, mod1 )
names(filt1)
filt1$y
filt1$mod
filt1$m
#
# Plot time series and fit.
#
plot ( LakeHuron, type="p", xlab="Year", ylab="Depth",
       main="Lake Huron" )
lines ( 1875:1972, filt1$m[-1] )
#
# Obtain confidence intervals 
#
var <- dlmSvd2var(filt1$U.C, filt1$D.C)
sd <- sqrt(unlist(var))
lines ( 1875:1972, filt1$m[-1] + 2*sd[-1], lty=3, col="red" )
lines ( 1875:1972, filt1$m[-1] - 2*sd[-1], lty=3, col="red" )
#
# Estimate  predicted values of the state vectors
#
filt1$a
dlmSvd2var ( filt1$U.R, filt1$D.R )
points ( 1875:1972, filt1$a, lty=2, col="green" ) # forecast
points ( 1875:1972, filt1$f, lty=2, col="magenta" ) # forecast
#
# Model 2: using fixed values of variances. Do V and W matter?
#
mod2 <- dlmModPoly ( order=1, dV=10, dW=1 )
filt2 <- dlmFilter ( LakeHuron, mod2 )
mod3 <- dlmModPoly ( order=1, dV=1, dW=10 )
filt3 <- dlmFilter ( LakeHuron, mod3 )
plot ( LakeHuron, type="p", xlab="Year", ylab="Depth",
       main="Lake Huron" )
lines ( 1875:1972, filt1$m[-1], col="red" )
lines ( 1875:1972, filt2$m[-1], col="blue" )
lines ( 1875:1972, filt3$m[-1], col="green" )
#
# Model 3: Try using discount factors!
#
dlmFilterDF <- function (y, mod, simplify = FALSE, DF) 
{
  ## storage.mode(y) <- "double"
  mod1 <- mod
  yAttr <- attributes(y)
  ytsp <- tsp(y)
  y <- as.matrix(y)
  timeNames <- dimnames(y)[[1]]
  stateNames <- names(mod$m0)
  m <- rbind(mod$m0, matrix(0, nr = nrow(y), nc = length(mod$m0)))
  a <- matrix(0, nr = nrow(y), nc = length(mod$m0))
  f <- matrix(0, nr = nrow(y), nc = ncol(y))
  U.C <- vector(1 + nrow(y), mode = "list")
  D.C <- matrix(0, 1 + nrow(y), length(mod$m0))
  U.R <- vector(nrow(y), mode = "list")
  D.R <- matrix(0, nrow(y), length(mod$m0))
  U.W <- vector(nrow(y), mode = "list")
  D.W <- matrix(0, nrow(y), length(mod$m0))
  Wliste <- vector(nrow(y), mode = "list")
  P <- vector(nrow(y), mode = "list")
  tmp <- La.svd(mod$V, nu = 0)
  Uv <- t(tmp$vt)
  Dv <- sqrt(tmp$d)
  Dv.inv <- 1/Dv
  Dv.inv[abs(Dv.inv) == Inf] <- 0
  sqrtVinv <- Dv.inv * t(Uv)
  sqrtV <- Dv * Uv
  tmp <- La.svd(mod$C0, nu = 0)
  U.C[[1]] <- t(tmp$vt)
  D.C[1, ] <- sqrt(tmp$d)
  for (i in seq(length = nrow(y))) {
    tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
    a[i, ] <- mod$GG %*% m[i, ]
    P[[i]] <- mod$GG %*% crossprod(D.C[i,] * t(U.C[[i]])) %*% t(mod$GG)
    Wliste[[i]] <- P[[i]]* ((1-DF)/DF)
    svdW <- La.svd( Wliste[[i]] , nu = 0)
    sqrtW <- sqrt(svdW$d) * svdW$vt
    U.W[[i]] <- t(svdW$vt)
    D.W[i, ] <- sqrt(svdW$d)
    tmp <- La.svd(rbind(D.C[i, ] * t(mod$GG %*% U.C[[i]]), 
                        sqrtW), nu = 0)
    U.R[[i]] <- t(tmp$vt)
    D.R[i, ] <- tmp$d
    f[i, ] <- mod$FF %*% a[i, ]
    D.Rinv <- 1/D.R[i, ]
    D.Rinv[abs(D.Rinv) == Inf] <- 0
    tmp <- La.svd(rbind(sqrtVinv %*% mod$FF %*% U.R[[i]], 
                        diag(x = D.Rinv, nrow = length(D.Rinv))), nu = 0)
    U.C[[i + 1]] <- U.R[[i]] %*% t(tmp$vt)
    foo <- 1/tmp$d
    foo[abs(foo) == Inf] <- 0
    D.C[i + 1, ] <- foo
    m[i + 1, ] <- a[i, ] + crossprod(D.C[i + 1, ] * t(U.C[[i
                                                           + 1]])) %*% tF.Vinv %*% as.matrix(y[i, ] - f[i,])
  }        
  m <- drop(m)
  a <- drop(a)
  f <- drop(f)
  attributes(f) <- yAttr
  ans <- list(m = m, U.C = U.C, D.C = D.C, a = a, U.R = U.R, 
              D.R = D.R, f = f, U.W=U.W, D.W=D.W)
  ans$m <- drop(ans$m)
  ans$a <- drop(ans$a)
  ans$f <- drop(ans$f)
  attributes(ans$f) <- yAttr
  if (!is.null(ytsp)) {
    tsp(ans$a) <- ytsp
    tsp(ans$m) <- c(ytsp[1] - 1/ytsp[3], ytsp[2:3])
    class(ans$a) <- class(ans$m) <- if (length(mod$m0) > 
                                        1) 
      c("mts", "ts")
    else "ts"
  }
  if (!(is.null(timeNames) && is.null(stateNames))) {
    dimnames(ans$a) <- list(timeNames, stateNames)
    dimnames(ans$m) <- list(if (is.null(timeNames)) NULL else c("", 
                                                                timeNames), stateNames)
  }
  if (simplify) 
    return(c(mod = list(mod1), ans))
  else {
    attributes(y) <- yAttr
    ans <- c(y = list(y), mod = list(mod1), ans)
    class(ans) <- "dlmFiltered"
    return(ans)
  }
}
mod3 <- dlmModPoly (order=1, dV=1)
modFilt <- dlmFilterDF(LakeHuron, mod3, DF=0.9)
plot.ts(LakeHuron)
lines(modFilt$m,col="red")
modFilt <- dlmFilterDF(LakeHuron, mod3, DF=0.8)
lines(modFilt$m,col="blue")
modFilt <- dlmFilterDF(LakeHuron, mod3, DF=0.95)
lines(modFilt$m,col="green")
#
# Model 4: local linear regression model
#
mod4 <- dlmModPoly(dV = 1, dW = c(0, 0.5))
modFilt <- dlmFilter(LakeHuron, mod4)
plot.ts(LakeHuron)
lines(modFilt$f,col="blue")
#
# Model 4: Try running the first order polynomial model with unknown variances using FFBS. This is slow!
#
gibbsout <- dlmGibbsDIG(LakeHuron,mod=dlmModPoly(1),shape.y=0.1,rate.y=0.1,shape.theta=0.1,rate.theta=0.1,n.sample=5000,thin=10)
burn <- 500
attach(gibbsout)
#
# Check convergence of running means
#
ts.plot(ergMean(dV[-burn]),ylab=expression(bar(V)),xlab='iterations')
ts.plot(ergMean(dW[-burn]),ylab=expression(bar(V)),xlab='iterations')
#
# Check the ACFs
#
acf(dV[-burn],lwd=3,col='red')
acf(dW[-burn],lwd=3,col='red')
#
# Plot posterior distributions of V and W
#
hist(dV[-burn],probability=TRUE,nclass=20,xlab='V',ylab='f',main='')
lines(density(dV[-burn]),lwd=3,col='red')
hist(dW[-burn],probability=TRUE,nclass=20,xlab='W',ylab='f',main='')
lines(density(dW[-burn]),lwd=3,col='red')
#
# Model 5: Simple bootstrap filter model for out of sample prediction.
#
thetagen <- function(y,m0,C0,V,W,N){
  theta <- rnorm(N,m0,sqrt(C0))
  w <- rep(1/N,N)
  n <- length(y)
  ypred <- matrix(rep(0,3*length(y)),nrow=length(y))
  for (i in 1:n){
    theta <- rnorm(N,theta,sqrt(W))
    err <- rnorm(1,0,sqrt(V))
    ypred[i,] <- quantile(theta+err,probs=c(0.05,0.5,0.95))
    w <- dnorm(y[i],theta,sqrt(V))
    w1 <- w/sum(w)
    theta = sample(theta,size=N,replace=TRUE,prob=w1)
  }
  return(list("theta"=theta,"ypred"=ypred,"w"=w))
}
cc <- thetagen(LakeHuron,600,100,0.2,0.5,10000)
thetan <- cc$theta
ypred <- cc$ypred
hist(thetan,probability=TRUE,xlab=expression(theta[n]),ylab='f',nclass=20,main='')
lines(density(thetan),col='blue',lwd=3)
#
# Compare with true posterior density for theta[n].
#
m0 <- 600
C0 <- 100
V <- 0.2
W <- 0.5
n <- length(LakeHuron)
C <- C0
m <- m0
for (i in 1:n){
  R <- C+W
  Q <- R+V
  A <- R/Q
  e <- LakeHuron[i]-m
  C <- R-(A^2)*Q
  m <- m+A*e
}
lines(seq(min(thetan),max(thetan),length.out=1000),dnorm(seq(min(thetan),max(thetan),length.out=1000),m,sqrt(C)),col='red',lwd=3)
#
# Look also at out of sample predictives.
#
ts.plot(LakeHuron,col='red',lwd=3)
lines(c(1875:1972),ypred[,1],col='blue',lwd=3,lty=3)
lines(c(1875:1972),ypred[,2],col='blue',lwd=3)
lines(c(1875:1972),ypred[,3],col='blue',lwd=3,lty=3)