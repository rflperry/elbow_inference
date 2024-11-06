source("./scripts/functions/hypothesis_tests.R")

conf_interval_solver <- function(
    test, vals, k, n, p, sigma, bounds=c(0, Inf), alpha=0.05, eigen=FALSE, root_lb=-20, root_ub=20, num_stab=0
  ) {
  alpha <- alpha / 2  # Since two sided CI
  
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
  
  int_ub <- NA
  try({
    for (i in seq(1, root_ub)) {
      test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=i)
      int_ub <- i
    }
  }, silent=TRUE)
  
  int_lb <- NA
  try({
    for (i in seq(-1, root_lb)) {
      test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=i)
      int_lb <- i
    }
  }, silent=TRUE)
  
  if ( is.na(int_ub) | is.na(int_lb) ) {
    warning('No bounds are valid')
  }

  if (test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=int_lb) > (1 - alpha)) {
    # If integration LB is above 1-alpha, we cannot construct CI
    ci_lb <- NA
    ci_ub <- NA
  } else if (test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=int_lb) > alpha) {
    # If we can't integrate LB low enough, set to negative infinity
    ci_lb <- -Inf
  } else if (test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=int_ub) < alpha) {
    # Otherwise, we need to chech that the upper bound makes sense
    ci_lb <- NA
    ci_ub <- NA
  } else {
    # We can integrate lower bound
    ci_lb <- uniroot(
      function(theta) {
        test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=theta) - alpha
      },
      interval = c(int_lb, int_ub)
    )$root
  }
  
  if (!is.na(ci_lb)) {
    # If we have a valid lower bound, need to express the upper bound
    if (test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=int_ub) < (1 - alpha)) {
      ci_ub <- Inf
    } else {
      ci_ub <- uniroot(
        function(theta) {
          test(vals, k, n, p, sigma, bounds=bounds, eigen=FALSE, theta=theta) - (1 - alpha)
        },
        interval = c(int_lb, int_ub)
      )$root
    }
  }
  
  return(c(ci_lb, ci_ub))
}

get_pve_ci <- function(ci_lower, ci_upper) { #, val_k, val_sum) {
  
  if (is.na(ci_lower) | is.na(ci_upper) ) {
    return(c(NA, NA))
  }
  
  if ((ci_lower <= 0) & (ci_upper >=0 )) {
    sq_ci_lower = 0
  } else {
    sq_ci_lower <- min( ci_lower^2, ci_upper^2 )
  }
  sq_ci_upper <- max( ci_lower^2, ci_upper^2 )
  # sq_ci_lower <- if_else(
  #   sign(ci_lower) == sign(ci_upper),
  #   pmin( ci_lower^2, ci_upper^2 ), 0)
  # sq_ci_upper <- pmax( ci_lower^2, ci_upper^2 )
  # estimator_ci_lower <- sq_ci_lower / ( val_sum - val_k + sq_ci_lower)
  # if (sq_ci_upper == Inf | sq_ci_lower == Inf) {
  #   estimator_ci_upper <- 1
  # } else {
  #   estimator_ci_upper <- sq_ci_upper / ( val_sum - val_k + sq_ci_upper)
  # }
  
  
  return(c( sq_ci_lower, sq_ci_upper ))
}

get_nc_chsq_ci <- function(frob_hat, df, alpha) {
  ci_lb <- uniroot(
    function(frob) {
      pchisq(frob_hat, df = df, ncp = frob) - (1 - alpha / 2)
    },
    interval = c(0, n*p*100)
  )$root


  ci_ub <- uniroot(
    function(frob) {
      pchisq(frob_hat, df = df, ncp = frob) - ( alpha / 2)
    },
    interval = c(0, df*100)
  )$root

  return(c(ci_lb, ci_ub))
}
