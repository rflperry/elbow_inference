source("./scripts/functions/hypothesis_tests.R")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Randomized hypothesis test and selection tools
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

get_randomized_svals <- function(X, tau, sigma=1) {
  # Parameters
  # ----------
  # X : data matrix
  # tau : randomization noise fraction, relative to variance sigma
  # sigma : data variance (default=1)
  #
  # Returns
  # -------
  # duv : SVD of noisy data matrix
  #
  n <- dim(X)[1]
  p <- dim(X)[2]
  noise <- array(rnorm(n*p, sd=tau*sigma), dim=c(n, p))
  return(svd(X + noise))
}

randomized_integrand <- function(s, svals, rand_svals, Lk, Lr, r, n, p, sigma, tau, lb, ub, num_const=0, theta=0) {
  f_integrand <- function(s_rand, s, const=0) {
    # Integrating over \tilde{s}
    # from c to next singular value
    f_int_helper <- function(s_rand) {
      exp(
        -(1/tau**2) / 2 * s_rand**2 +
          (1/tau**2) * s_rand * (s * Lr[r] + sum(svals[-r] * Lr[-r]) ) +
          (n - p) * log(s_rand) + sum(log(abs(rand_svals[-r]**2 - s_rand**2))) -
          const # sign(const) * log(const)
      ) # * s_rand ** (n - p) * prod(abs(rand_svals[-r]**2 - s_rand**2))
    }
    
    # num_const <- 0 # f_int_helper(rand_svals[r])
    return(
      sapply(
        s_rand,
        function(t) {
          f_int_helper(t)
        }
      )
    )
  }
  
  g_int_helper <- function(s, const=0) {
    exp(
      -(1/sigma**2 + 1/tau**2) / 2 * s**2 +
        s * ( (theta / sigma**2) + sum(rand_svals[-r] * Lk[-r]) ) +
        (n - p) * log(s) + sum(log(abs(svals[-r]**2 - s**2))) -
        const
    )
  }
  num_const <- 0 # g_int_helper(svals[r])

  return(
    sapply(
      s,
      function(t) {
        f_term <- integrate(f_integrand, lb, ub, t, num_const)$value
        g_term <- g_int_helper(t, num_const)
        return(f_term * g_term)
      }
    )
  )
}


randomized_pvalue <- function(duv, rand_duv, r, n, p, sigma, tau, lb=0, ub=Inf, eigen=FALSE, theta=0) {
  tau <- sigma * tau
  svals <- duv$d
   rand_svals <- rand_duv$d

   Lk <- diag(t(rand_duv$v) %*% duv$v * t(rand_duv$u) %*% duv$u)
   Lr <- (t(rand_duv$v[,r]) %*% duv$v * t(rand_duv$u[,r]) %*% duv$u)[1,]

   if (r == 0 || r == (p+1)) { return(NA) }
   ub <- min(ub, c(Inf, rand_svals, 0)[r])
   lb <- max(lb, c(Inf, rand_svals, 0)[r+2])
   
   if (eigen) {
     svals <- sqrt(svals)
     rand_svals <- sqrt(rand_svals)
     ub <- sqrt(ub)
     lb <- sqrt(lb)
   }
   
   # numerator <- integrate(
   #   randomized_integrand, svals[r], ub,
   #   svals, rand_svals, Lk, Lr, r, n, p, sigma, tau, lb=lb)$value
   # denominator <- integrate(
   #   randomized_integrand, lb, ub,
   #   svals, rand_svals, Lk, Lr, r, n, p, sigma, tau, lb=lb)$value

   # pvalue = 1 / (\int_c^r / \int_r^{r-1} + 1)
   numerator <- integrate(
     randomized_integrand, c(Inf, svals, 0)[r+2], svals[r],
     svals, rand_svals, Lk, Lr, r, n, p, sigma, tau, lb=lb, ub=ub)$value
   denominator <- integrate(
     randomized_integrand, svals[r], c(Inf, svals, 0)[r],
     svals, rand_svals, Lk, Lr, r, n, p, sigma, tau, lb=lb, ub=ub)$value
   pvalue <- 1 / (numerator/ denominator + 1)
   
   if(!is.na(pvalue) && (pvalue < 0 || pvalue > 1)){
     stop(paste0("Randomized pvalue invalid ", pvalue))
   }
   
   return(pvalue)
}

# 
# n <- 10
# p <- 5
# X1 <- array(rnorm(n*p, sd=sigma), dim=c(n, p))
# X2 <- array(rnorm(n*p, sd=sigma*0.2), dim=c(n, p))
# duv <- svd(X1)
# rand_duv <- svd(X2 + X1)
# r <- 3
# 
# randomized_pvalue(duv, rand_duv, r, n, p, 1, 0.2)
# 
# s_rand <- rand_svals[r]
# s <- svals[r]
# exp(
#   -(1/tau**2) / 2 * s_rand**2 +
#     (1/tau**2) * s_rand * (s * Lr[r] + sum(svals[-r] * Lr[-r]) ) +
#   (n - p) * log(s_rand) + sum(log(abs(rand_svals[-r]**2 - s_rand**2)))
# )

#* s_rand ** (n - p) * prod(abs(rand_svals[-r]**2 - s_rand**2))

# 
# i <- 1
# for (i in seq(1, 5)) {
#   j <- 2
#   # print(t(rand_duv$v[,i]) %*% duv$v[,j] %*% t(rand_duv$u[,i]) %*% duv$u[,j])
#   # print(sum(diag(duv$v[,j] %*% t(duv$u[,j]) %*% rand_duv$u[,i] %*% t(rand_duv$v[,i]))))
#   print(t(rand_duv$v[,i]) %*% duv$v[,j] * t(duv$u[,j]) %*% rand_duv$u[,i])
# }
# 
# t(t(diag(t(rand_duv$v) %*% duv$v))) %*% diag(t(duv$u) %*% rand_duv$u)
# 
# duv$v[,j] %*% t(duv$u[,j]) %*% rand_duv$u[,i] %*% t(rand_duv$v[,i])
