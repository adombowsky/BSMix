# hellinger distance for two univariate normal distributions
unorm.hellinger <- function(mu_1, mu_2, sigma_1, sigma_2){
  # mu_1, mu_2 = means of 2 normals
  # sigma_1, sigma_2 = standard deviation of two normals
  h_sq = 1 - sqrt( (2*sigma_1 * sigma_2)/(sigma_1^2 + sigma_2^2)) * exp(-(mu_1 - mu_2)^2/(4 * (sigma_1^2 + sigma_2^2)))
  h = sqrt(h_sq)
  return(h)
}