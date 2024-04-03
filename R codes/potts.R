#
# Potts model and similar.
#
rm(list = ls())
if (!require("mrf2d") == TRUE){
  install.packages("mrf2d")
}
library(mrf2d)
#
# Sample from a Potts model.
#
th <- expand_array(-1, family = "onepar", mrfi(1), C = 1)
z_sample <- rmrf2d(init_Z = c(200,200), mrfi = mrfi(1), theta = th)
cplot(z_sample)
#
# Sample with border constraints.
#
border <- matrix(FALSE, nrow = 100, ncol = 100)
border[1, ] <- border[100, ] <- border[, 1] <- border[, 100] <- TRUE
initial <- matrix(sample(0:1, 100*100, replace = TRUE), nrow = 100, ncol = 100)
initial[border] <- 0
z_border <- rmrf2d(initial, mrfi = mrfi(1), theta = th, fixed_region = border)
cplot(z_border)
#
# Inference on a real data set.
#
data("bold5000", package = "mrf2d")
cplot(bold5000)
Rnn <- mrfi(1)
theta_nn <- expand_array(-1, family = "onepar", C = 3, mrfi = Rnn)
fit_brain <- fit_ghm(bold5000, Rnn, theta_nn, equal_vars = TRUE)
fit_brain_ind <- fit_ghm(bold5000, Rnn, theta_nn*0, equal_vars = TRUE)
summary(fit_brain)
summary(fit_brain_ind)
cplot(fit_brain$Z_pred, legend = TRUE)
cplot(fit_brain_ind$Z_pred, legend = TRUE)