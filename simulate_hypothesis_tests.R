library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis) 
library(stringr)
library(scales)
library(latex2exp)
library(argparse)

source("./scripts/functions/hypothesis_tests.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")
source("./scripts/functions/estimation.R")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("-n", type = "integer", default=50)
parser$add_argument("-p", type = "integer", default=10)
parser$add_argument("--rank", type = "integer", default=5)
parser$add_argument("-m", type = "integer", default=1)
parser$add_argument("--reps", type = "integer", default=1000)
parser$add_argument("--alpha", type = "double", default=0.1)
parser$add_argument("--method", choices = c("zg", "elbow"), default="zg")
parser$add_argument("--sigmas", 
                    type = "double", 
                    nargs = "+",
                    default = c(0.1, 0.2, 0.4, 0.5, 0.7, 1))
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Run Simulation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

eigen = TRUE

# sigmas <- c(0.1, 0.15, 0.2, 0.3, 0.4, 0.7, 1)
selection_rule = "zg" # options include "zg", "elbow"
if (selection_rule == "zg") {
  # sigmas <- c(0.1, 0.2) # c(0.1, 0.2, 0.4, 0.5, 0.7, 1)
  degree <- 0.4
} else if (selection_rule == "elbow") {
  # sigmas <- c(0.1, 0.2, 0.3, 1, 2)
  degree <- 1
}

# Output args to global variables
list2env(as.list(args), envir = .GlobalEnv)

fname <- "sim_hypo_tests"

for (arg_name in names(args)) {
  if (!arg_name %in% c("sigmas")) {
    arg_value <- args[[arg_name]]
    fname <- paste0(fname, "_", arg_name, "=", if (is.vector(arg_value)) paste(arg_value, collapse = "-") else arg_value)
  }
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
      sim <- simulate_matrix(n, p, sigma, rank, m, degree=degree, eigen=eigen)
      
      # Extract quantities
      duv <- svd(sim$obsv_mat)
      
      if (eigen) {
        vals <- duv$d^2
      } else {
        vals <- duv$d
      }

      frob_norm <- sqrt(sum(sim$mean_mat^2))

      noise_op_norm <- (svd(sim$obsv_mat - sim$mean_mat)$d)[1]
      
      # Perform selection
      if (selection_rule == "zg") {
        r <- select_r_zg(vals)
      } else if (selection_rule == "elbow") {
        r <- select_r_elbow(vals)[1]
      }
      
      for (k in seq(1, r)) {
        # Base quantities to save
        tr_mean <- (t(duv$u[,k]) %*% sim$mean_mat %*% duv$v[,k])[1]
        
        # Compute intervals
        base_results <- c(rep, n, p, sigma, m, rank, r, k, tr_mean, frob_norm, vals[k], sum(vals), sqrt(vals[1]), noise_op_norm)
        
        # Choi p-value
        choi_pvalue <- choi_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma, eigen=eigen)
        results <- c(results, 'Choi', base_results, choi_pvalue)

        # Acquire bounds for selective p-value
        if (selection_rule == "zg") {
          bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)
        } else if (selection_rule == "elbow") {
          bounds <- compute_elbow_bounds(vals, k)
        }

        si_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma, bounds=bounds, eigen=eigen)
        results <- c(results, 'Selective inference', base_results, si_pvalue)
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
  "method", "rep", "n", "p", "sigma", "signal_strength",
  'rank', 'selection_r', 'tested_k', 'signal', "frob_norm", 'val_k', 'val_sum', 'mean_op_norm', 'noise_op_norm', "p_value")

results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
head(results_df)
print(nrow(results_df[complete.cases(results_df),]) / nrow(results_df))

save(results_df, file=paste0("data/", fname, ".RData"))