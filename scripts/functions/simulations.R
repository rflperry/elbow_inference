#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Functions to simulate data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

simulate_global_null <- function (n, p, sigma){
  return(array(rnorm(n*p, sd=sigma), dim=c(n, p)))
}

simulate_matrix <- function (n, p, sigma, rank, m, degree=1, eigen=TRUE, offset=0){
  UV <- array(rnorm(n*p, sd=sigma), dim=c(n, p))
  duv <- svd(UV)
  U <- duv$u
  V <- duv$v
  vals <- c(rev(seq(1:(rank))), rep(0, p-rank))^degree * m * (n*p)^(1/4) + offset
  if (eigen) {
    vals <- sqrt(vals)
  }
  mean_mat <- U %*% diag(vals) %*% t(V)
  results <- {}
  results$mean_mat <- mean_mat
  results$obsv_mat <- mean_mat + array(rnorm(n*p, sd=sigma), dim=c(n, p))
  return(results)
}

simulate_Choi <- function (n, p, sigma, rank, m, degree=1, eigen=FALSE){
  UV <- array(rnorm(n*p, sd=sigma), dim=c(n, p))
  duv <- svd(UV)
  U <- duv$u
  V <- duv$v
  vals <- c(rev(seq(1:(rank))), rep(0, p-rank))^degree * m * (n*p)^(1/4) * sigma
  if (eigen) {
    vals <- sqrt(vals)
  }
  mean_mat <- U %*% diag(vals) %*% t(V)
  results <- {}
  results$mean_mat <- mean_mat
  results$obsv_mat <- mean_mat + array(rnorm(n*p, sd=sigma), dim=c(n, p))
  return(results)
}