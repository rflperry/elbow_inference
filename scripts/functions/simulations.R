#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Functions to simulate data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

simulate_global_null <- function (n, p, sigma){
  return(array(rnorm(n*p, sd=sigma), dim=c(n, p)))
}

simulate_matrix <- function (n, p, sigma, rank, m, degree=1, eigen=TRUE, offset=0, thin_c=NA){
  UV <- array(rnorm(n*p, sd=sigma), dim=c(n, p))
  duv <- svd(UV)
  U <- duv$u
  V <- duv$v
  # vals <- c(rev(seq(1:(rank))), rep(0, p-rank))^degree * m * (n*p)^(1/4) + offset
  if ( rank < min(n, p)) {
    vals <- c(seq(0, 5, length.out=rank+1), rep(0, p-rank-1))^degree * m * (50*10)^(1/4) + offset
  } else {
    vals <- (seq(0, 5, length.out=rank+1)[2:(rank+1)])^degree * m * (50*10)^(1/4) + offset
  }
  
  if (eigen) {
    vals <- sqrt(vals)
  }
  mean_mat <- U %*% diag(vals) %*% t(V)
  results <- {}
  results$mean_mat <- mean_mat
  H <- array(rnorm(n*p, sd=sigma), dim=c(n, p))
  
  if(!is.na(thin_c)) {
    results$obsv_mat <- mean_mat + H * sqrt(1 + c^2)
    results$obsv_mat2 <- mean_mat - H * sqrt(1 + 1/c^2)
  } else {
    results$obsv_mat <- mean_mat + H
  }
  
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