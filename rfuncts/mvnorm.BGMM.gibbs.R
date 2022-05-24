gmm.gibbs <- function(R, B, y, alpha, kappa, w, nu, M, K, stops = 1) {
  # Gibbs sampler from McLachlan and Peel
  require(dirmult)
  require(mvtnorm)
  require(LaplacesDemon)
  
  # K = number of clusters
  # alpha = Dirichlet parameter
  # kappa = scaling parameter for GMM prior
  # w = matrix of mean parameters for mu's
  # nu = first parameter in Wishart distribution
  # M = second parameter in Wishart distribution
  
  # prelims
  n <- nrow(y)
  p <- ncol(y)
  # storage
  pi <- matrix(NA, nrow = R, ncol = K)
  mu <- list() # list of matricies for each cluster
  for (i in 1:R) {
    mu[[i]] <- matrix(0, nrow = K, ncol = p)
  }
  Sigma <- list() # list of list, each inner list containing Sigma_k samples
  for (r in 1:R) {
    Sigma[[r]] <- list() 
  }
  c <- matrix(NA, nrow = R, ncol = n) # cluster labels
  cluster.means <- matrix(NA, nrow = K, ncol = p)
  # initialization
  Z <- t(rmultinom(n, 1, rep(1/K, K)))
  c[1, ] <- apply(Z, 1, function(x) which.max(x == 1))
  M_inv <- list()
  tau <- c()
  for (k in 1:K){
    M_inv[[k]] <- solve(M[[k]])
    Sigma[[1]][[k]] <- diag(p)
  }
  # configuring stops
  stops <- (1:(R/stops)) * stops
  # sampling
  print("Sampling")
  for (r in 2:R) {
    # compute cluster sizes
    cluster.sizes <- as.vector(table(c[r-1, ]))
    # sample pi
    alpha_star <- alpha + cluster.sizes
    pi[r, ] <- rdirichlet(1, alpha_star)
    # sample normal parameters
    for (k in 1:K){
      # mu
      cluster <- matrix(y[c[r-1, ] == k, ])
      cluster.means[k, ] <- colMeans(cluster)
      w_star = (cluster.sizes[k] * cluster.means[k, ] + kappa[k] * w[k, ])/(cluster.sizes[k] + kappa[k])
      mu[[r]][k, ] <- rmvnorm(1, w_star, (1/(cluster.sizes[k] + kappa[k])) * Sigma[[r-1]][[k]])
      # Sigma
      V <- tcrossprod(cluster[1, ] - cluster.means[k, ], cluster[1, ] - cluster.means[k, ])
      for (j in 2:cluster.sizes[k]){
        V <- V + tcrossprod(cluster[j, ] - cluster.means[k, ], cluster[j, ] - cluster.means[k, ])
      }
      V = (1/cluster.sizes[k]) * V
      M_star = M_inv[[k]] + cluster.sizes[k] * V + ( (cluster.sizes[k] * nu[k])/(cluster.sizes[k] + nu[k])) * 
        tcrossprod(cluster.means[k, ] - w[k, ], cluster.means[k, ] - w[k, ])
      M_star = round(M_star, 6)
      Sigma[[r]][[k]] <- rinvwishart(cluster.sizes[k] + nu[k], M_star)
    }
    # sample multinomial vectors and cluster labels
    for (i in 1:n) {
      for (k in 1:K) {
        tau[k] <- pi[r, k] * dmvnorm(x = y[i, ], mean = mu[[r]][k, ], sigma = Sigma[[r]][[k]])
      }
      tau <- tau/sum(tau)
      Z[i, ] <- rmultinom(1,1,tau)
    }
    c[r, ] <- apply(Z, 1, function(x) which(x == 1))
    
    # print stops
    if (r %in% stops){
      print(r)
    }
  }
  
  # discard Burn-in
  pi <- pi[(B+1):R, ]
  mu <- mu[(B+1):R]
  Sigma <- Sigma[(B+1):R]
  c <- c[(B+1):R, ]
  return(list(pi = pi, mu = mu, Sigma = Sigma, c = c))
}