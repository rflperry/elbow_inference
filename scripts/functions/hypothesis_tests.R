#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Integrands and p-value computation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

# Integration integrand, proportional to likelihood.
# No reason to variable transform here since can be done earlier with the bounds.
# Note that integration requires log=FALSE. Set to TRUE only if 
integrand <- function(val_obsv, val_other, n, p, sigma, eigen=FALSE, theta=0, num_stab=0) {
  sapply(
    val_obsv, 
    function(val) {
      sval <- val
      sval_other <- val_other
      log_val <- (-sval**2 / 2 + sval * theta) / sigma**2 - num_stab
      log_val <- log_val + ((n-p) * log(sval) + sum(log(abs(sval_other**2 - sval**2))))
      return (exp (log_val))
    }
  )
}

compute_pvalue_integral <- function(vals, k, n, p, sigma, bounds, eigen=FALSE, theta=0) {
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
  # 3a) if upperbound < s_k, numerator = 0
  # 3b) else if lower bound < s_k, \numerator = int(s_k, upperbound)
  # 3c) else numerator = denominator
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
    
    denominator <- denominator + denom_integral$value
    denom_error <- denom_error + denom_integral$abs.error
    
    # Step 3
    if (ub < val_obs) {
      numerator <- numerator + 0
    } else if (lb < val_obs) {
      num_integral <- integrate(
        integrand, val_obs, ub, vals[-k], n, p, sigma,
        eigen=eigen, theta=theta, num_stab=num_stab
      )
      numerator <- numerator + num_integral$value
      num_error <- num_error + num_integral$abs.error
    } else {
      numerator <- numerator + denom_integral$value
      num_error <- num_error + denom_integral$abs.error
    }
  }
  
  # try({
  if ( ((numerator - num_error) / (denominator + denom_error) ) > 1) {
    warning(cat("Invalid numerator/denominator:", numerator, '/', denominator))
  }
  # })

  if (denominator == 0) {
    pvalue = 1
  } else {
    pvalue <- numerator / denominator
  }
  if (pvalue > 1) {
    # Previous check ensures that this is only due to numerical error
    pvalue <- 1
  }
# 
#   # pvalue = 1 / (\int_c^r / \int_r^{r-1} + 1)
#   num_stab <- val_obs * theta / sigma**2
#   numerator <- integrate(integrand, lb, val_obs, vals[-k], n, p, sigma, eigen=eigen, theta=theta, num_stab=num_stab)$value
#   denominator <- integrate(integrand, val_obs, ub, vals[-k], n, p, sigma, eigen=eigen, theta=theta, num_stab=num_stab)$value
#   pvalue <- 1 / (numerator/ denominator + 1)

  return(pvalue)
}

# Choi pvalue that the kth eigenvalue is greater than zero
# bounds parameter unused and solely for compatibility
choi_test_pvalue <- function(vals, k, n, p, sigma, bounds=NA, eigen=FALSE, theta=0) {
  pvalue <- compute_pvalue_integral(
    vals, k, n, p, sigma, bounds=c(0, Inf), eigen=eigen, theta=theta
    )
  # if(!is.na(pvalue) && (pvalue < 0 || pvalue > 1+1e-5)){
  #   stop(paste0("Choi pvalue invalid ", pvalue))
  # }
  return(pvalue)
}

# SI pvalue that the rth (smallest eigenvalue greater than the cutoff)
# is greater than zero
si_test_pvalue <- function(vals, k, n, p, sigma, bounds=c(0, Inf), eigen=FALSE, theta=0) {
  pvalue <- compute_pvalue_integral(
    vals, k, n, p, sigma, bounds=bounds, eigen=eigen, theta=theta
    )
  # if(!is.na(pvalue) && (pvalue < 0 || pvalue > 1+1e-5)){
  #   stop(paste0("SI pvalue invalid ", pvalue))
  # }
  return(pvalue)
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	P-value wrappers with alpha levels, and other tests
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

choi_test <- function(svals, k, sigma=1, alpha=0.05) {
  pvalue <- choi_test_pvalue(svals, k, sigma=1)
  return(as.integer(pvalue < alpha))
}

test_pvalues <- function(svals, r, sigma=1) {
  bounds <- c(Inf, svals, 0)
  s_obs <- svals[r]
  numerator <- integrate(integrand, s_obs, bounds[r], svals[-r], sigma)$value
  choi_denominator <- integrate(integrand, bounds[r+2], bounds[r], svals[-r], sigma)$value
  si_denominator <- integrate(integrand, cutoff, bounds[r], svals[-r], sigma)$value
  choi_pvalue <- numerator / choi_denominator
  si_pvalue <- numerator / si_denominator
  return(c(choi_pvalue, si_pvalue, si_denominator / choi_denominator))
}

si_test <- function(svals, cutoff, sigma=1, alpha=0.05) {
  pvalue <- si_test_pvalue(svals, cutoff, sigma=1)
  return(as.integer(pvalue < alpha))
}

# Our sequential tests using the Choi pvalue with a multiple correction test
choi_seq_test <- function(svals, cutoff, sigma=1, alpha=0.05) {
  r <- which.min(c(svals >= cutoff, 0)) - 1
  if (r == 0) {
    return(NA)
  }
  return(choi_test(svals, r, sigma, alpha=alpha/(r+1)))
}

# muirhead test under null rank <= k-1, i.e. s_k \neq 0 (reject implies rank >= k)
muirhead_test_pvalue <- function(evals, k, n, p) {
  if (k == p) {
    return(1)
  }
  q <- p - k + 1
  lq <- sum((evals)[k:p]) / q
  vk <- prod(evals[k:p]) / lq^q # 
  test_stat <- -(n - k + 1 - (2*q^2 + q + 2) / (6*q) + sum(lq^2 / (evals[1:(k-1)] - lq)^2))*log(vk)
  return( pchisq(test_stat, df=(q+2)*(q-1)/2, lower.tail=FALSE) )
}

# pseudo rank test under null rank <= k-1, i.e. s_k \neq 0 (reject implies rank >= k)
pseudorank_test_pvalue <- function(evals, k, sn, p) {
  mup <- (sqrt(n - 1/2) + sqrt(p - k - 1/2))^2
  sigmap <- (sqrt(n - 1/2) + sqrt(p - k - 1/2)) * (1/sqrt(n - 1/2) + 1/sqrt(p - k - 1/2))^(1/3)
  test_stat <- (evals[k] - mup) / sigmap
  return( ptw(test_stat, lower.tail=FALSE) )
}