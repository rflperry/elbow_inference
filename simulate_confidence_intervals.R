library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(scales)
library(latex2exp)
library(argparse)
library(RMTstat)

source("./scripts/functions/confidence_intervals.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")
source("./scripts/functions/estimation.R")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("-n", type = "integer", default = 50)
parser$add_argument("-p", type = "integer", default = 10)
parser$add_argument("-m", type = "integer", default = 1)
parser$add_argument("--rank", type = "integer", default = 5)
parser$add_argument("--reps", type = "integer", default = 1000)
parser$add_argument("--alpha", type = "double", default = 0.1)
parser$add_argument("--method", choices = c("zg", "elbow"), default = "zg")
parser$add_argument("--sigmas",
  type = "double",
  nargs = "+",
  default = c(0.1, 0.2)
)
parser$add_argument("--c", type = "double", default = sqrt(1))
parser$add_argument("--choi", type = "logical", default = FALSE)
parser$add_argument("--mle", type = "logical", default = TRUE)
parser$add_argument("--var_est", action='store_true', default = TRUE)
parser$add_argument("--signal_alpha_frac", type = "double", default = 0.75)
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Run Simulation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

eigen <- TRUE

# sigmas <- c(0.1, 0.15, 0.2, 0.3, 0.4, 0.7, 1)
selection_rule <- "zg" # options include "zg", "elbow"
if (selection_rule == "zg") {
  # sigmas <- c(0.1, 0.2) # c(0.2, 0.4, 0.7, 1)
  degree <- 0.4
} else if (selection_rule == "elbow") {
  # sigmas <- c(0.1, 0.2, 0.3, 1, 2)
  degree <- 1
}

# Output args to global variables
list2env(as.list(args), envir = .GlobalEnv)

fname <- "sim_conf_ints"

for (arg_name in names(args)) {
  if (!arg_name %in% c("sigmas")) {
    arg_value <- args[[arg_name]]
    fname <- paste0(fname, "_", arg_name, "=", if (is.vector(arg_value)) paste(arg_value, collapse = "-") else arg_value)
  }
  # Format the output: prepend the variable name
  # cat(arg_name, ":", if (is.vector(arg_value)) paste(arg_value, collapse = " ") else arg_value, "\n")
}
print(fname)

results <- c()

for (rep in 1:reps) {
  if (rep %% 100 == 0) {
    print(paste0("Rep: ", rep))
  }
  for (sigma in sigmas) {
    try({
      # Simulate data
      sim <- simulate_matrix(n, p, sigma, rank, m, degree = degree, eigen = eigen, thin_c = c)

      # Estimate sigma if specified
      if (var_est) {
        # #  Transform data instead
        # I <- diag(1, n) # Identity matrix
        # J <- matrix(1 / n, n, n) # Matrix with all elements equal to 1/n
        # C <- I - J # Centering matrix
        # duv <- svd(C)

        # H <- duv$u %*% sqrt(diag(duv$d))
        # sim$obsv_mat <- t(H) %*% as.matrix(sim$obsv_mat)
        # SCALAR <- max(sim$obsv_mat) / 2 * sqrt(n)
        # sim$obsv_mat <- sim$obsv_mat / SCALAR

        # sim$mean_mat <- t(H) %*% as.matrix(sim$mean_mat) / SCALAR

        duv <- svd(sim$obsv_mat)
        vals <- duv$d
        sigma_hat <- sqrt(median(sqrt(vals))^2 / (max(n, p) * qmp(0.5, svr = max(n, p) / min(n, p))))

      } else{
        # duv <- svd(sim$obsv_mat)
        # vals <- duv$d
        sigma_hat <- sigma
      }

      sigma1 <- sqrt(sigma_hat^2 * (1 + c^2))
      sigma2 <- sqrt(sigma_hat^2 * (1 + 1 / c^2))

      # Thin
      W <- array(rnorm(n*p, sd=sigma_hat), dim=c(n, p))
      sim$obsv_mat1 <- sim$obsv_mat + W * sqrt(1 + c^2)
      sim$obsv_mat2 <- sim$obsv_mat - W * sqrt(1 + 1/c^2)

      # Estimate svals on thin 1
      duv <- svd(sim$obsv_mat1)
      vals <- duv$d
      if (eigen) {
        vals <- duv$d^2
      }

      X2_frob_norm <- sqrt(sum(sim$obsv_mat2^2))
      frob2_hat <- X2_frob_norm^2 - n * p * sigma2^2
      frob_norm <- sqrt(sum(sim$mean_mat^2))

      frob_ci <- sqrt(sigma2^2 * get_nc_chsq_ci(X2_frob_norm^2 / sigma2^2, n * p, alpha * (1 - signal_alpha_frac)))
      frob_mle <- sqrt(sigma2^2 * get_mle_nc_chisq(X2_frob_norm^2 / sigma2^2, n * p))

      # Perform selection
      if (selection_rule == "zg") {
        r <- select_r_zg(vals)
      } else if (selection_rule == "elbow") {
        r <- select_r_elbow(vals)[1]
      }

      for (k in seq(1, r)) {
        # Base quantities to save
        tr_mean <- (t(duv$u[, k]) %*% sim$mean_mat %*% duv$v[, k])[1]

        # Compute intervals
        base_results <- c(rep, n, p, sigma, c, m, rank, r, k, tr_mean, frob_norm, vals[k], sum(vals), X2_frob_norm, frob_ci, frob_mle, frob2_hat, sigma_hat)

        if (choi) {
          # Choi p-value
          ci <- conf_interval_solver(
            choi_test_pvalue, vals, k, n, p,
            sigma = sigma1, alpha = alpha * signal_alpha_frac, eigen = eigen
          )
          # ci[1] <- ci[1] / frob_ci[2]
          # ci[2] <- ci[2] / frob_ci[1]

          choi_mle <- NaN
          if (mle) {
            try({
              choi_mle <- get_mle_signal(
                vals, k, n, p, sigma1,
                eigen = eigen
              )
            })
          }

          # choi_median <- get_median_signal(
          #   choi_test_pvalue, vals, k, n, p,
          #   sigma = sigma1, eigen = eigen
          # )

          results <- c(results, "Choi", base_results, ci, choi_mle)
        }

        # Selective inference p-value (R >= k)
        # Acquire bounds
        if (selection_rule == "zg") {
          bounds <- compute_constraint_regions(vals, k, select_r_zg, tol = 1000)
        } else if (selection_rule == "elbow") {
          bounds <- compute_elbow_bounds(vals, k)
        }

        si_mle <- NaN
        if (mle) {
          try({
            si_mle <- get_mle_signal(
              vals, k, n, p, sigma1,
              bounds = bounds, eigen = eigen
            )
          })
        }
        # si_median <- get_median_signal(
        #   si_test_pvalue, vals, k, n, p,
        #   sigma = sigma1, bounds = bounds, eigen = eigen
        # )

        ci <- conf_interval_solver(
          si_test_pvalue, vals, k, n, p,
          sigma = sigma1, alpha = alpha * signal_alpha_frac, bounds = bounds, eigen = eigen
        )
        results <- c(results, "Selective inference", base_results, ci, si_mle)
      }
    })
  }
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Process data and save
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

header <- c(
  "method", "rep", "n", "p", "sigma", "thin_c", "signal_strength",
  "rank", "selection_r", "tested_k", "signal", "frob_norm", "val_k", "val_sum",
  "X2_frob_norm", "frob_ci_lower", "frob_ci_upper", "frob_mle", "frob2_hat", "sigma_hat", "ci_lower", "ci_upper", "mle"
)

results_df <- data.frame(
  data = t(array(results, dim = c(length(header), length(results) / length(header))))
)
colnames(results_df) <- header
results_df[, 2:length(header)] <- sapply(results_df[, 2:length(header)], as.numeric)
head(results_df)
print(nrow(results_df[complete.cases(results_df), ]) / nrow(results_df))

save(results_df, file = paste0("data/", fname, ".RData"))
