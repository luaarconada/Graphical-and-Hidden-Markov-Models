#
# Importance sampling and sampling importance resampling (SIR).
#
# Available functions:
#
# imp_sampler_example: uses a sample from a uniform distribution to estimate the 
#   mean of a beta distribution.
# SIR_example: uses sample importance resampling from a uniform distribution to 
#   generate an approximate sample from a beta distribution.
# importance_sampler_failure_example: tries to estimate the mean and variance 
#   of a heavy tailed distribution (t) by sampling from a light tailed 
#   distribution (normal).
# SIR_failure_example: tries to estimate a normal density by sampling from a 
#   shifted t distribution.  If the shift is large, the number of unique data is 
#   very small
# bayesian_SIR_example: estimates the posterior distribution for a model with a
#   non-conjugate prior.
#
imp_sampler_example <- function(n = 1000, pars = c(5, 2), plot = FALSE){
  #
  # Uses a uniform distribution to estimate the mean of a beta distribution. 
  #
  if (is.vector(n) == FALSE || length(n) > 1 || class(n) != "numeric" || 
      is.na(n) == TRUE || is.infinite(n) == TRUE || n != round(n) || n <= 0){
    stop("Error: n must be a positive integer.")
  }
  if (is.vector(pars) == FALSE || length(pars) > 2 || class(pars) != 
      "numeric" || any(pars <= 0) == TRUE || any(is.na(pars) == TRUE) == TRUE || 
      any(is.infinite(pars) == TRUE)){
    stop("Error: pars must be a vector containing 2 non-negative numbers.")
  }
  if (is.vector(plot) == FALSE || length(plot) > 1 || class(plot) != 
      "logical"){
    stop("Error: plot must be either TRUE or FALSE")
  }
  #
  # Generate a sample from a uniform distribution and calculate the 
  # unnormalized weights, that is the true density divided by the uniform 
  # density at each point.
  #
  y <- runif(n)
  weights <- dbeta(y,pars[1], pars[2])
  #
  # If necessary, plot a histogram of the sampled data along with the sampling 
  # density and the true density.  This shows that the sample data don't look 
  # like the true density and motivates the use of SIR.
  #
  if (plot == TRUE){
    ygrid <- c(0:1000)/1000
    ftrue <- dbeta(ygrid, pars[1], pars[2])
    cc <- hist(y, breaks = ceiling(sqrt(n)), plot = FALSE)
    ylim <- c(0,max(1, ftrue, cc$density))
    hist(y, breaks = ceiling(sqrt(n)), freq = FALSE, xlim = c(0,1), ylim = ylim, 
        xlab = "y", ylab = "f", main = "")
    lines(ygrid, ftrue, lwd = 2, col = "blue")
    lines(c(0,1), c(1,1), lwd = 2, col = "green")
    legendplace <- ifelse(pars[1] > pars[2], "topleft", "topright")
    legend(legendplace,c("Sampling density", "True density"), col = c("green", 
                                                              "blue"), lwd = 2)
  }
  #
  # Calculate the normalized weights and the effective sample size.
  #
  normalized_weights <- weights/sum(weights)
  ess <- 1/sum(normalized_weights^2)
  return(list("True mean" = pars[1]/sum(pars), "Estimated mean" = 
        sum(normalized_weights*y), "Effective sample size" = ess))
}
#
SIR_example <- function(n = 1000, pars = c(5,2)){
  #
  # Estimates the density of a Beta distribution using sampling importance 
  # resampling.
  #
  if (is.vector(n) == FALSE || length(n) > 1 || class(n) != "numeric" || 
      is.na(n) == TRUE || is.infinite(n) == TRUE || n != round(n) || n <= 0){
    stop("Error: n must be a positive integer.")
  }
  if (is.vector(pars) == FALSE || length(pars) > 2 || class(pars) != 
      "numeric" || any(pars <= 0) == TRUE || any(is.na(pars) == TRUE) == TRUE || 
      any(is.infinite(pars) == TRUE)){
    stop("Error: pars must be a vector containing 2 non-negative numbers.")
  }
  #
  # Generate a sample from a uniform distribution and calculate the 
  # unnormalized weights, that is the true density divided by the uniform 
  # density at each point. Then resample with probabilities proportional to the
  # weights.
  #
  y <- runif(n) 
  weights <- dbeta(y, pars[1], pars[2]) 
  y <- sample(y, n, TRUE, weights) 
  #
  # Plot a histogram of the resampled 
  # values and the estimated and true densities.
  #
  gridy <- c(0:1000)/1000
  truefy <- dbeta(gridy, pars[1], pars[2])
  cc <- hist(y, breaks = ceiling(sqrt(n)), plot = FALSE)
  dd <- density(y, from = 0, to = 1)
  ylim = c(0,max(cc$density, truefy, dd$y))
  hist(y, breaks = ceiling(sqrt(n)), freq = FALSE, ylim = ylim, xlab = 
         "y", ylab = "f", main = "")
  lines(dd$x, dd$y, lwd = 2, col = "red")
  lines(gridy, truefy, lwd = 2, col = "blue")
  legendplace <- ifelse(pars[1] > pars[2], "topleft", "topright")
  legend(legendplace,c("Estimated density", "True density"), col = c("red", 
                                                              "blue"), lwd = 2)
}
#
importance_sampler_failure_example <- function(n = 1000, df = 3){
  #
  # Tries to estimate the mean and variance of a t distribution using 
  # values sampled from a normal. You will see that the estimate of the mean is 
  # not too bad, but the variance is severely underestimated.  For an importance
  # sampler, it is necessary that the distribution used for sampling has longer 
  # tails than the true density.
  #
  if (is.vector(n) == FALSE || length(n) > 1 || class(n) != "numeric" || 
      is.na(n) == TRUE || is.infinite(n) == TRUE || n != round(n) || n <= 0){
    stop("Error: n must be a positive integer.")
  }
  if (is.vector(df) == FALSE || length(df) > 1 || class(df) != "numeric" || 
      is.na(df) == TRUE || is.infinite(df) == TRUE || df <= 0 || df <= 2){
    stop("Error: df must be a number greater than 2.")
  }
  #
  # Generate a sample from the closest normal distribution to the desired t 
  # distribution and calculate the associated, unnormalized weights as the ratio 
  # of the normal to the t.
  #
  sigma <- sqrt(df/(df-2)) 
  y <- rnorm(n, 0, sigma) 
  weights <- dt(y, df)/dnorm(y, 0, sigma)
  #
  # Return the true and estimated means and variances and effective sample size.
  #
  normalized_weights <- weights/sum(weights)
  est_mean <- sum(normalized_weights*y)
  est_var <- sum(normalized_weights*y^2)-est_mean^2
  return(list("True mean and variance" = c(0,df/(df-2)), 
              "Estimated mean and variance" = c(est_mean,est_var), 
              "Effective sample size " = 1/sum(normalized_weights^2)))
}
#
SIR_failure <- function(n = 1000, shift = 5, df = 4){
  #
  # Tries to use SIR to estimate the density of a normal(0,1) distribution using 
  # values sampled from a shifted t. The final sample is very small!
  #
  if (is.vector(n) == FALSE || length(n) > 1 || class(n) != "numeric" || 
      is.na(n) == TRUE || is.infinite(n) == TRUE || n != round(n) || n <= 0){
    stop("Error: n must be a positive integer.")
  }
  if (is.vector(df) == FALSE || length(df) > 1 || class(df) != "numeric" || 
      is.na(df) == TRUE || is.infinite(df) == TRUE || df <= 0 || df <= 2){
    stop("Error: df must be a number greater than 2.")
  }
  if (is.vector(shift) == FALSE || length(shift) > 1 || class(shift) != 
      "numeric" ||is.na(shift) == TRUE || is.infinite(shift) == TRUE || 
      abs(shift) > 5){
    stop("Error: shift must be a real number with absolute value <= 5.")
  }
  #
  # Generate an initial sample from a normal distribution with the same variance
  # as the t distribution, calculate the associated weights and draw a histogram
  # of these.  The weight distributin will be very skewed.
  #
  inity <- shift+rt(n, df) # Sample from the closest normal distribution.
  weights <- dnorm(inity)/dt(inity-shift, df) # Calculate the weights.
  hist(weights/sum(weights), xlab = "normalised weights", ylab = "f", main = 
         "Histogram of normalized weights")
  #
  # Resample to try and get an approximate sample from the true distribution and
  # draw a histogram of the sampled values along with the true and estimated 
  # densities.
  #
  y <- sample(inity, n, TRUE, weights)
  cc <- hist(y, breaks = ceiling(sqrt(n)), plot = FALSE)
  dd <- density(y)
  gridy <- -4+c(0:1000)*8/1000
  truefy <- dnorm(gridy)
  fsample <- dt(gridy-shift,df)
  ylim = c(0,max(cc$density, truefy, dd$y, fsample))
  hist(y, breaks = ceiling(sqrt(n)), freq = FALSE,  xlim = c(-4,4), 
       ylim = ylim, xlab = "y", ylab = "f", main = "")
  lines(dd$x, dd$y, lwd = 2, col = "red")
  lines(gridy, truefy, lwd = 2, col = "blue")
  lines(gridy, fsample, lwd = 2, col = "green")
  legendplace <- ifelse(shift < 0, "topleft", "topright")
  legend(legendplace,c("Estimated density", "True density", "Sampling density"), 
         col = c("red", "blue", "green"), lwd = 2)
  nuniquey <- length(unique(y))
  #
  # Returns the effective sample size of the importance sampler and the number 
  # of unique values in the final sample.
  #
  normalized_weights <- weights/sum(weights)
  return(list("effective sample size" = 1/sum(normalized_weights^2), 
                  "number of unique y" = nuniquey))
}
#
bayesian_SIR_example <- function(n = 1000, size = 12, x = 9, mu = 0.5, 
                                 sigma = 1){
  #
  # Produces a sample from the posterior distribution of p given x, where 
  # X | size, p ~ Binomial(size, p) and p ~ Truncated N(mu, sigma).
  # Requires: truncnorm.
  #
  if (is.vector(n) == FALSE || length(n) > 1 || class(n) != "numeric" || 
      is.na(n) == TRUE || is.infinite(n) == TRUE || n != round(n) || n <= 0){
    stop("Error: n must be a positive integer.")
  }
  if (is.vector(size) == FALSE || length(size) > 1 || class(size) != 
      "numeric" || is.na(size) == TRUE || is.infinite(size) == TRUE || 
      size != round(size) || size <= 0){
    stop("Error: size must be a positive integer.")
  }
  if (is.vector(x) == FALSE || length(x) > 1 || class(x) != "numeric" || 
      is.na(x) == TRUE || is.infinite(x) == TRUE || x != round(x) || x < 0 || 
      x > size){
    stop("Error: x must be a non-negative integer <= size.")
  }
  if (is.vector(mu) == FALSE || length(mu) > 1 || class(mu) != "numeric" ||
      is.na(mu) == TRUE || is.infinite(mu) == TRUE){
    stop("Error: mu must be a real number.")
  }
  if (is.vector(sigma) == FALSE || length(sigma) > 1 || class(sigma) != 
      "numeric" || is.na(sigma) == TRUE || is.infinite(sigma) == TRUE || 
      sigma <= 0){
    stop("Error: sigma must be a positive, real number.")
  }
  if (!require("truncnorm") == TRUE){
    install.packages("truncnorm")
  }
  #
  # Generate a sample from a uniform distribution and calculate the unnormalized
  # weights.
  #
  p <- runif(n)
  weights <- dbinom(x,size,p)*truncnorm::dtruncnorm(p, 0, 1, mu, sigma)
  #
  # Estimate the integrating constant. 
  #
  est_ic <- mean(weights)
  #
  # Calculate the normalized weights and effective sample size and estimate the
  # posterior mean of p.
  #
  normalized_weights <- weights/(est_ic*n)
  ess <- 1/sum(normalized_weights^2)
  est_mean_p <- sum(normalized_weights*p)
  #
  # Calculate the integrating constant and posterior mean via numerical 
  # integration.
  #
  fic <- function(p){
    f <- truncnorm::dtruncnorm(p,0,1,mu,sigma)*dbinom(x,size,p)
  }
  ni_ic <- integrate(fic,0,1)$value
  fmeanp <- function(p){
    f <- p*fic(p)/ni_ic
  }
  ni_mean_p <- integrate(fmeanp,0,1)$value
  #
  # Resample and plot a histogram of the sampled values, the estimated density 
  # and the density estimated using numerical integration.
  #
  p <- sample(p, n, TRUE, weights)
  gridp <- c(0:1000)/1000
  f_ni <- fic(gridp)/ni_ic
  f_est <- density(p)
  cc <- hist(p, breaks = ceiling(sqrt(n)), plot = FALSE)
  ylim <- c(0,max(c(f_ni,f_est$y,cc$density)))
  hist(p, breaks = ceiling(sqrt(n)), freq = FALSE, xlim = c(0,1), ylim = ylim,
                           xlab = "p", ylab = "f(p|x)", main = "")
  lines(gridp, f_ni, lwd = 2, col = "blue")
  lines(f_est$x, f_est$y, lwd = 2, col = "red")
  legendplace <- ifelse(est_mean_p > 0.5, "topleft", "topright")
  legend(legendplace,c("SIR density", "NI density"), 
         col = c("red", "blue"), lwd = 2)
  #
  # Return the true and estimated integrating constant and mean.
  #
  return(list("True integrating constant and mean" = 
            c(ni_ic,ni_mean_p), "Estimated integrating constant and mean" = 
              c(est_ic, est_mean_p)))
}
