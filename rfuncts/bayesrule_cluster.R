# for clustering the data, for use with lapply
bayesrule_cluster <- function(new_mixtures, y, K) {
  mixtures <- new_mixtures$mixtures
  L_mix <- unlist(lapply(mixtures, nrow))
  mixtures <- do.call(rbind, mixtures)
  labels <-  rep(1:K, L_mix)
  clusters <- sapply(y, cluster_dnpmix, mixtures = mixtures, L_mix = L_mix, K = K, labels = labels)
  return(clusters)
}

cluster_dnpmix <- function(x, mixtures, L_mix, K, labels){
  # evalues newly weighted mixture
  dens <- rep(0,K)
  for (k in 1:K){
    mixtures_k <- mixtures[ labels == k,]
    if (!is.matrix(mixtures_k)) {
      mixtures_k <- matrix(mixtures_k, nrow = 1)
    }
    dens_comps <- apply(mixtures_k, 1, function(z) z[1]  * dnorm(x, mean = z[2], sd = z[3]))
    dens[k] <- sum(dens_comps)
  }
  return(which.max(dens))
}