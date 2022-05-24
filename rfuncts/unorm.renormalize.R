unorm.renormalize <- function(mcmc_samps, L, K){
  mu <- mcmc_samps[1:L]
  sigma <- mcmc_samps[(L+1):(2*L)]
  omega <- mcmc_samps[(2*L + 1):(3*L)]
  dens_clusters <- mcmc_samps[(3*L + 1): (4*L)]
  grouped_mixtures <- list()
  grouped_weights <- c()
  for (k in 1:K){
    omega_k <- omega[dens_clusters == k]/sum(omega[dens_clusters == k])
    mu_k <- mu[dens_clusters == k]
    sigma_k <- sigma[dens_clusters == k]
    n_k <- length(omega_k)
    grouped_pars <- cbind(omega_k, mu_k, sigma_k)
    grouped_mixtures[[k]] <- grouped_pars
    grouped_weights[k] <- sum(omega[dens_clusters == k])
  }
  return(list(mixtures = grouped_mixtures, weights = grouped_weights))
}