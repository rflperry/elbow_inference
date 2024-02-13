source("./scripts/functions/hypothesis_tests.R")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	MLE computation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

proportional_density <- function(sval, sval_other, n, p, sigma, theta, num_stab=0, log=TRUE) {
  log_val <- (-sval**2 / 2 + sval * theta) / sigma**2 - num_stab
  log_val <- log_val + ((n-p) * log(sval) + sum(log(abs(sval_other**2 - sval**2))))
  if (log) {
    return ( log_val )
  } else {
    return (exp (log_val))
  }
}

compute_neg_loglik <- function(vals, k, n, p, sigma, bounds, eigen=FALSE, theta=0) {
  # vals are in descending order
  if (k == 0 || k == (p+1)) { return(NA) }
  
  if ( !is.list(bounds) ) {
    bounds <- list(bounds)
  }
  
  if (eigen) {
    vals <- sqrt(vals)
    bounds <- lapply(bounds, sqrt)
  }
  # bounds due to sval_ordering
  sval_ub <- c(Inf, vals, 0)[k]
  sval_lb <- c(Inf, vals, 0)[k+2]
  
  val_obs <- vals[k]
  num_stab <- (- val_obs**2  / 2 + val_obs * theta) / sigma**2 
  
  numerator <- 0
  denominator <- 0
  
  num_error <- 0
  denom_error <- 0
  # Given a list of bounds for s_k, and [k-1, k+1] sval bounds:
  # 1a) if bounds are outside of sval bounds, skip
  # 1b) else, truncate bounds by sval bounds
  # 2) compute denominator density in bounds
  for (i in seq(1, length(bounds))) {
    # Step 1a
    ub <- bounds[[i]][2]
    lb <- bounds[[i]][1]
    if (ub < sval_lb || lb > sval_ub) { next }
    
    # Step 1b
    ub <- min(ub, sval_ub)
    lb <- max(lb, sval_lb)
    
    # Step 2
    denom_integral <- integrate(
      integrand, lb, ub, vals[-k], n, p, sigma,
      eigen=eigen, theta=theta, num_stab=num_stab
    )
    
    denominator <- denominator + denom_integral$value # * exp(num_stab) # * exp( num_stab )
    denom_error <- denom_error + denom_integral$abs.error # * exp( num_stab )
    
    # if (length(bounds) > 1) {
    #   denominator <- denominator + exp(denom_integral$value + num_stab) # * exp( num_stab )
    #   denom_error <- denom_error + denom_integral$abs.error # * exp( num_stab )
    # } else {
    #   denominator <- denom_integral$value + num_stab # * exp( num_stab )
    #   denom_error <- denom_integral$abs.error # * exp( num_stab )
    # }
  }
  # if (length(bounds) > 1) {
  #   denominator <- log(denominator)
  # }
  
  # Step 3
  numerator <- proportional_density(vals[k], vals[-k], n, p, sigma, theta=theta, num_stab = num_stab, log=TRUE)
  
  neg_loglik <- - ( numerator - log(denominator) )
  
  return( neg_loglik )
}

get_mle <- function(vals, k, n, p, sigma, bounds=c(0, Inf), eigen=eigen) {

  fn <- function(theta) {
    nll <- compute_neg_loglik(
      vals, k, n, p, sigma, bounds=bounds, eigen=eigen, theta=theta
    )
    return( nll )
  }
  
  opt <- NA
  try({
    opt <- optim(0, fn, method = "L-BFGS-B")$par
  }, silent = TRUE)
  if (is.na(opt)) {
    warning("L-BFGS-B optimizer failed, using robust Nelder and Mead optimizer")
    opt <- optim(0, fn)$par
  }

  #, method = "L-BFGS-B", lower = 0, upper = Inf)
  # opt <- optimize(fn, interval = c(-10, 10))
  return(opt)
}


get_median <- function(
  test, vals, k, n, p, sigma, bounds=c(0, Inf), eigen=FALSE, root_lb=-20, root_ub=20, num_stab=0
) {
  alpha <- 0.5  # median
  
  edit_bounds <- function(lbub) {
    ub <- min(lbub[2], c(Inf, vals, 0)[k])
    lb <- max(lbub[1], c(Inf, vals, 0)[k+2])
    
    if (eigen) {
      lb <- sqrt(lb)
      ub <- sqrt(ub)
    }
    return(c(lb, ub))
  }
  
  if ( !is.list(bounds) ) {
    bounds <- edit_bounds(bounds)
  }
  else {
    bounds <- lapply(bounds, edit_bounds)
  }
  
  if (eigen) {
    vals <- sqrt(vals)
  }
  
  int_ub <- 0
  try({
    for (i in seq(1, root_ub)) {
      test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=i)
      int_ub <- i
    }
  }, silent=TRUE)
  
  int_lb <- 0
  try({
    for (i in seq(-1, root_lb)) {
      test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=i)
      int_lb <- i
    }
  }, silent=TRUE)
  
  if (test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=int_lb) > (1 - alpha)) {
    # If integration LB is above 1-alpha, we cannot construct CI
    ci_lb <- NA
  } else if (test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=int_lb) > alpha) {
    # If we can't integrate LB low enough, set to negative infinity
    ci_lb <- -Inf
  } else if (test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=int_ub) < alpha) {
    # Otherwise, we need to chech that the upper bound makes sense
    ci_lb <- NA
  } else {
    # We can integrate lower bound
    ci_lb <- uniroot(
      function(theta) {
        test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=theta) - alpha
      },
      interval = c(int_lb, int_ub)
    )$root
  }
  
  return(ci_lb)
}

get_pve <- function(theta, val_k, val_sum) {
  pve <- theta / ( val_sum - val_k + theta)
  return( pve )
}


# fn <- function(theta) {
#   return(
#     compute_neg_loglik(
#       vals, k, n, p, sigma, bounds=bounds, eigen=eigen, theta=theta
#     )
#   )
# }
# 
# thetas <- seq(-1, 5, 0.1)
# neg_liks <- lapply(thetas, fn)
# plot(thetas, neg_liks)
