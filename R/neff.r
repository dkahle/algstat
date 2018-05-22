# Effective Sample Size
# 
# Calculates the effective sample size of one MCMC chain. This function is very similar to Stan function rstan:::ess_rfun. 
#
# 
# the parameter chain should be the output of the metropolis function
#
#

neff <- function(chain) {
  PRs <- computeUProbsCpp(chain$steps)
  chain_length <- length(PRs)
  auto_cov <- acf(PRs, lag.max = chain_length - 1, plot = FALSE, type = "covariance")$acf[,,1]
  mean_var <- auto_cov[1] * chain_length / (chain_length - 1)
  rho_hat_sum <- 0
  for (t in 2:chain_length) {
    rho_hat <- 1 - (mean_var - auto_cov[t]) / auto_cov[1]
    if (is.nan(rho_hat)) rho_hat <- 0
    if (rho_hat < 0) break
    rho_hat_sum <- rho_hat_sum + rho_hat
  }
  ess <- chain_length
  if (rho_hat_sum > 0) ess <- ess / (1 + 2 * rho_hat_sum)
ess
}