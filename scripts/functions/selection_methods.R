#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Define selection rule bounds
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#



zg_lik <- function(vals, q) {
  p <- length(vals)
  vals1 <- vals[1:q]
  vals2 <- vals[(q+1):p]
  
  # mean ests
  mu1 <- mean(vals1)
  mu2 <- mean(vals2)
  
  # pooled var
  var1 <- var(vals1)
  var2 <- var(vals2)
  pooled_var <- ( (q - 1) * var1 + (p - q - 1) * var2 ) / (p - 2)
  
  log_lik1 <- -q/2 * log(2 * pi * pooled_var) - (1/(2 * pooled_var)) * sum((vals1 - mu1)^2)
  log_lik2 <- -(p-q)/2 * log(2 * pi * pooled_var) - (1/(2 * pooled_var)) * sum((vals2 - mu2)^2)
  
  return(log_lik1 + log_lik2)
}


select_r_zg <- function(vals, return_liks=FALSE) {
  # computes logliks, will include NAs but that's okay
  logliks <- lapply(seq(1, length(vals)), zg_lik, vals=vals)
  
  if (return_liks) {
    return(c(which.max(logliks), logliks))
  } else {
    return(which.max(logliks))
  }
}


# Function to calculate the menger curvature, which we wish to minimize with the elbow
menger_curvature <- function(x1, x2, x3, y1, y2, y3) {
  tri_area <- abs( (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) ) / 2
  
  # Calculate the menger curvature
  dist_prod <- sqrt( (x1 - x2)^2 + (y1 - y2)^2 ) *
    sqrt( (x1 - x3)^2 + (y1 - y3)^2 ) *
    sqrt( (x2 - x3)^2 + (y2 - y3)^2 )
  
  return(4*tri_area / dist_prod)
}

select_r_menger <- function(vals) {
  p <- length(vals)
  mcs <- c()
  for ( i in seq(2, p-1) ) {
    mcs <- c(mcs, menger_curvature(1, i, p, vals[1], vals[i], vals[p]))
  }
  r <- which.max(mcs) # for the selected index, but cancels with desired target
  
  return(c(r, mcs))
}


select_r_variance <- function(vals, threshold) {
  r <- which.max( (cumsum(vals) / sum(vals) - threshold) >= 0)
  return(r)
}

select_r_cutoff <- function(vals, cutoff) {
  r <- which.min(c(vals >= cutoff, 0)) - 1
  return(c(r, cutoff, Inf))
}

compute_lineq_bounds <- function(A, b, x, r) {
  # Given the system Ax <= b, solves for the interval of values which x[r] can
  # take within this solution set, holding all other values fixed
  #
  # Parameters
  # ----------
  # A : matrix, dim=c(p, n)
  # b : vector, dim=c(p)
  # x : vector, dim=c(n)
  # r : integer, in [1, n]
  #
  # Returns
  # -------
  # c(lb, ub) : upper and lower bounds which x[r] can take satisfying Ax <= b

  ub_indices <- which(A[, r] > 0)
  lb_indices <- which(A[, r] < 0)
  
  if (length(ub_indices) == 0) {
    ub <- Inf
  } else {
    ub <- min( (b[ub_indices] - A[ub_indices, -r] %*% x[-r]) / A[ub_indices, r])
  }
  
  if (length(lb_indices) == 0) {
    lb <- 0
  } else {
    lb <- max( (b[lb_indices] - A[lb_indices, -r] %*% x[-r]) / A[lb_indices, r])
  }
  
  # If solution set doesn't exist
  if (lb > ub) {
    lb <- x[r]
    ub <- x[r]
  }

  return(c(lb, ub))
}

compute_elbow_bounds_broken <- function(vals, k, inequality=TRUE) {
  # Computes bounds on the kth value such that
  # the maximized index is always >= k
  
  if (k == 1) {
    return(c(vals[2], Inf))
  }
  
  # elbow matrix
  p <- length(vals)
  A <- matrix(0, p-2, p)
  for (i in 1:nrow(A)) {
    # objective function
    A[i, i:(i+2)] <- c(1, -2, 1)
  }
  
  kappas <- as.vector(A %*% vals)

  if (!inequality) {
    # We want bounds on s_k such that elbow_r is the maximum, i.e. 
    # Ax <= A_r x <=> (A - A_r)x <= 0 for index k
    r <- which.max(kappas)
    Ar <- A
    Ar[,(r):(r+2)] <- t(t(Ar[,(r):(r+2)]) - c(1, -2, 1))
    bounds <- compute_lineq_bounds(Ar, rep(0, p-2), vals, k)
    return(bounds)
  }
  
  # Two sets of inequalities, either of which must hold
  # Note 1: from the start we work with only the first k-1 rows of A
  # Note 2: if we are studying s_k, then k <= r. If k < r, then elbow_r is independent of s_k
  Ak <- A[1:(k-1),]
  if (is.vector(Ak)) {
    Ak <- t(Ak)
  }
  # Computes the bounds of the kth value s_k can take such that
  # the maximum [k+1:p] elbow (dependent of s_k) is larger than
  # any [1:k-1] elbow (some dependent on s_k)
  # i.e. A_{1:k-1} x <= max_{h \in \{k+1:p\} (elbow_h)
  # Note only relevant if k < p-2
  if (k < p-2) {
    bounds1 <- compute_lineq_bounds(Ak, rep(max(kappas[(k+1):(p-2)]), k-1), vals, k)
  }
  
  # Computes the bounds of the kth value s_k can take such that
  # the kth elbow is larger than any elbow [1:k-1]
  # i.e. A_{1:k-1} x <= A_k x <=>  (A_{1:k-1} - A_k )x <= 0
  Ak[,(k):(k+2)] <- t(t(Ak[,(k):(k+2)]) - c(1, -2, 1))
  bounds2 <- compute_lineq_bounds(Ak, rep(0, k-1), vals, k)
  # stopifnot(bounds2[1] <= bounds2[2])
  
  if (k == p-2) {
    return(bounds2)
  } else {
    stopifnot(bounds1[1] <= bounds1[2])
  }

  if ( bounds2[2] < bounds1[1] || bounds1[2] < bounds2[1]) {
    # Case 1: disjoint intervals
    return( list(bounds1, bounds2) )
  } else {
    # Case 2: intersecting intervals
    # Since either (not both) must hold, min lower bound and max upper bound suffice.
    return(c( min(bounds1[1], bounds2[1]), max(bounds1[2], bounds2[2]) ))
  }
}

compute_elbow_bounds <- function(vals, k, verbose=FALSE) {
  # Computes bounds on the kth value such that
  # the maximized index is always >= k
  
  # Compute derivatives
  p <- length(vals)
  A <- matrix(0, p-2, p)
  for (i in 1:nrow(A)) {
    # objective function
    A[i, i:(i+2)] <- c(1, -2, 1)
  }
  
  kappas <- as.vector(A %*% vals)
  kappas <- c(-Inf, kappas, -Inf) # pad for indexing
  
  # Compute constants
  if (k >= 4) {
    c1 <- max(kappas[2:(k-2)])
  } else {
    c1 <- -Inf
  }
  if (k <= p-3) {
    c2 <- max(kappas[(k+2):(p-1)])
  } else {
    c2 <- -Inf
  }
  
  # First two cases, k=1 or k=2
  if (k == 1) {
    if (verbose) {
      return(list(c(0, Inf), 1))
    }
    return(c(0, Inf))
  } else if (k == 2) {
    interval <- c(
      min(
        1/3 * (vals[k-1] + 3*vals[k+1] - vals[k+2]),
        1/2 * (vals[k-1] + vals[k+1] - c2)
      )
      , Inf
    )
    if (verbose) {
      return(list(interval, 1))
    }
    return(interval)
  } 
  
  # Compute A and B
  a <- max(
    1/3 * (vals[k-1] + 3*vals[k+1] - vals[k+2]),
    c1 + 2*vals[k+1] - vals[k+2]
  )
  b <- 1/2 * (vals[k-1] + vals[k+1] - c2)
  c <- 2*vals[k-1] - vals[k-2] + c2
  
  # Conditions 1,2
  if ( !(c1 <= c2 & c2 > -Inf) ) {
    # If not condition 1 => A
    code <- 1
    bounds <- c(a, Inf)
  } else if ( !(vals[k-2] - 2*vals[k-1] <= -2*vals[k+1] + vals[k+2]) ) {
    # If not condition 2 => B
    code <- 2
    bounds <- c(b, c)
  } else if (c < b) {
    # stop("B interval ordering incorrect.")
    bounds <- c(a, Inf)
    code <- -1
  } else{
    if (a < c) {
      # If overlapping, union
      code <- 3
      bounds <- c(min(a, b), Inf)
    } else {
      # If disjoint
      code <- 4
      bounds <- list( c(b, c), c(a, Inf) )
    }
  }
  if (verbose) {
    return(list(bounds, code))
  }
  return(bounds)
}


compute_constraint_regions <- function(vals, i, func, tol=1000) {
  # Takes in a list of values and an index i
  # Takes in a function maping (vals, i) -> index k
  # Given a tolerance, finds all windows where k >= i
  
  if (i == 1) {
    # Unchanged by selection event
    return(list(c(vals[i+1], Inf)))
  }
  
  regions <- c()
  inv_tol <- (vals[i-1] - vals[i+1]) / tol

  is_valid = FALSE

  for ( check_val in seq(vals[i+1], vals[i-1], inv_tol) ) {
    k <- func(c(vals[1:(i-1)], check_val, vals[(i+1):length(vals)]))
    if (!is_valid & (i <= k )) {
      # If wasn't in a valid region but k >= i, append and make valid
      regions <- c(regions, check_val)
      is_valid <- TRUE
    } else if (is_valid & (i > k)) {
      # If was in a valid region but now k < i, append prior and make invalid
      regions <- c(regions, check_val - inv_tol)
      is_valid <- FALSE
    }
  }

  if (is_valid) {
    regions <- c(regions, vals[i-1])
  }
  region_mat <- matrix(regions, ncol=2, byrow=TRUE)
  region_list <- lapply(seq_len(nrow(region_mat)), function(i) region_mat[i,])

  return(region_list)

}


select_r_elbow <- function (vals) {
  # Inputs eigenvalues and selects that which is followed by the discrete derivative maximizer.
  # Computes bounds between which the selected eigenvalue can vary.
  # Returns: selection index, lower bounds, upper bounds, list of discrete derivatives.

  # elbow
  p <- length(vals)
  A <- matrix(0, p-2, p)
  for (i in 1:nrow(A)) {
    # objective function
    A[i, i:(i+2)] <- c(1, -2, 1)
  }
  kappas <- A %*% vals
  r <- which.max(kappas) # kappas on length p-2, so r corresponds to eval index r+1

  # compute bounds, inputting the matrix corresponding to the linear constraint As < b
  # for the selected s_r
  A[,(r):(r+2)] <- t(t(A[,(r):(r+2)]) - c(1, -2, 1))
  bounds <- compute_lineq_bounds(A, rep(0, p-2), vals, r)

  return(c(r, bounds[1], bounds[2], kappas))
}
