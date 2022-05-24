# creates distance matrix for L univariate normal distributions
unorm.D <- function(mu, sigma){
  # mu = length L vector of component means
  # sigma = length L vector of component sds
  pars <- cbind(mu, sigma)
  L <- length(mu)
  D <- matrix(0, nrow = L, ncol = L)
  for (l in 1:(L-1)){
    for (j in l:L) {
      D[l,j] = unorm.hellinger(mu_1 = pars[l, 1], mu_2 = pars[j, 1], sigma_1 = pars[l,2], sigma_2 = pars[j,2])
    }
  }
  D = D + t(D)
  return(D)
}