# clusters univariate normal distribution
unorm.cluster <- function(pars, L, K) {
  # pars = matrix where each row has 2*L entries, the first L being the means mu and the second being standard devs sigma
  # L = number of overfitted clusters
  # K = number of clusters
  
  mu <- pars[1:L]
  sigma <- pars[-(1:L)]
  
  # first, create distance matrix
  D <- unorm.D(mu = mu, sigma = sigma)
  D <- as.dist(D) # converting to distance matrix for hclust
  
  # single linkage clustering
  cl <- hclust(D, method = "single")
  labs <- cutree(cl, k = K) # cluster labels

  return(labs)
}