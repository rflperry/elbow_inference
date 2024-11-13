library(dplyr)
library(tidyr)

source("./scripts/functions/confidence_intervals.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Run simulation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

n = 50 # 50
p = 10 # 10
rank = 5
m = 1
eigen = TRUE
reps <- 10000
alpha <- 0.1

# sigmas_list <- c(0.1, 0.15, 0.2, 0.3, 0.4, 0.7, 1)
selection_rule = "zg" #  "elbow" # options include "zg", "elbow"
if (selection_rule == "zg") {
  sigmas_list <- c(0.2)
  degree <- 0.4
} else if (selection_rule == "elbow") {
  sigmas_list <- c(0.2)
  degree <- 1
}

c <- sqrt(1)
signal_alpha_frac <- 0.75

args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  selection_rule <- args[1]
  reps <- as.numeric(args[2])
  c <- sqrt(as.numeric(args[3]))
  signal_alpha_frac <- as.numeric(args[4])
}

results <- c()

for (rep in 1:reps) {
  for (sigma in sigmas_list) {
    try({
      set.seed(rep)
      # Simulate data
      sim <- simulate_matrix(n, p, sigma, rank, m, degree=degree, eigen=eigen, thin_c=c)
      
      # Extract quantities
      duv <- svd(sim$obsv_mat)
      
      if (eigen) {
        vals <- duv$d^2
      } else {
        vals <- duv$d
      }
      
      sigma1 <- sqrt( sigma^2 * (1 + c^2) )
      sigma2 <- sqrt( sigma^2 * (1 + 1/c^2) )
      
      frob_hat <- sqrt(sum(sim$obsv_mat2^2)) # - n * p * sigma2^2)
      frob2_consistent <- frob_hat^2 - n * p * sigma2^2
      frob_norm <- sqrt(sum(sim$mean_mat^2))
      
      # Deducer library
      # frob_ci <- sqrt( chi.noncentral.conf(frob_hat^2 / sigma2^2, df=n*p, conf=1-alpha/4) )

      frob_ci <- sqrt( sigma2^2 * get_nc_chsq_ci(frob_hat^2 / sigma2^2, n*p, alpha * ( 1- signal_alpha_frac)) )
      
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
        base_results <- c(rep, n, p, sigma, c, m, rank, r, k, tr_mean, frob_norm, vals[k], sum(vals), frob_hat, frob2_consistent, frob_ci)
        
        # Choi p-value
        ci <- conf_interval_solver(
          choi_test_pvalue, vals, k, n, p, sigma=sigma1, alpha=alpha * signal_alpha_frac, eigen=eigen
        )
        # ci[1] <- ci[1] / frob_ci[2]
        # ci[2] <- ci[2] / frob_ci[1]
        
        # Triangle inequality method
        # ci[1] <- ci[1] / ( frob_hat + sigma2 * eps )
        # ci[2] <- ci[2] / ( frob_hat - sigma2 * eps )
  
        results <- c(results, 'Choi', base_results, ci)
        
        # Choi p-value, w/ Bonferroni correction
        # ci <- conf_interval_solver(
        #   choi_test_pvalue, vals, k, n, p, sigma=sigma, alpha=alpha / (p-2), eigen=eigen
        # )
        # results <- c(results, 'Choi (Bf.)', base_results, ci)
        
        # Selective inference p-value (R >= k)
        # Acquire bounds
        if (selection_rule == "zg") {
          bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)
        } else if (selection_rule == "elbow") {
          bounds <- compute_elbow_bounds(vals, k)
        }
        
        # Choi p-value
        ci <- conf_interval_solver(
          si_test_pvalue, vals, k, n, p, sigma=sigma1, alpha=alpha * signal_alpha_frac, bounds=bounds, eigen=eigen
        )
        # ci[1] <- ci[1] / frob_ci[2]
        # ci[2] <- ci[2] / frob_ci[1]
        
        # Triangle inequality method
        # ci[1] <- ci[1] / ( frob_hat + sigma2 * eps )
        # ci[2] <- ci[2] / ( frob_hat - sigma2 * eps )

        results <- c(results, 'Selective inference', base_results, ci)
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
  'rank', 'selection_r', 'tested_k', 'signal', "frob_norm", 'val_k', 'val_sum', "X2_frob_norm", "frob2_hat",
  "frob_ci_lower", "frob_ci_upper", 'ci_lower', 'ci_upper')
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
head(results_df)
print(nrow(results_df[complete.cases(results_df),]) / nrow(results_df))

save(results_df, file=paste0("data/", selection_rule, "-c2=", c^2, "-alpha_frac=", signal_alpha_frac,  "-coverage_vs_k.RData"))

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

library(ggplot2)
library(viridis) 
library(stringr)
library(scales)
library(latex2exp)

# selection_rule <- "elbow"

theme_update(text = element_text(size=10, family="Times"))
load(paste0("data/", selection_rule, "-c2=", c^2, "-alpha_frac=", signal_alpha_frac, "-coverage_vs_k.RData"))

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    # pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[1, ],
    # pve_ci_upper = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[2, ],
    # pve = signal^2 / ( val_sum - val_k + signal^2),
    pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
    pve_ci_upper = mapply(get_pve_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2,
    pve = (signal / frob_norm )^2,
  )
head(results_df)

plot_df <- results_df %>%
  mutate(
    method = case_when(
      method == "Selective inference" ~ "Selective",
      method == "Choi" ~ "Unselective",
    ),
    precision = 1/sigma,
    pve_covered = as.numeric((pve_ci_lower <= pve) & (pve <= pve_ci_upper)),
    frob_covered = as.numeric((frob_ci_lower <= frob_norm) & (frob_norm <= frob_ci_upper)),
    signal_covered = as.numeric((signal <= pve) & (pve <= pve_ci_upper)),
    pve_ci_width = pve_ci_upper - pve_ci_lower,
  ) %>% 
  group_by(
    method, sigma, precision, tested_k
  ) %>%
  summarise(
    mean_pve_ci_lower = median(pve_ci_lower, na.rm=TRUE),
    mean_pve = median(pve, na.rm=TRUE),
    mean_pve_ci_upper = median(pve_ci_upper, na.rm=TRUE),
    coverage = mean(pve_covered, na.rm=TRUE),
    frob_coverage = mean(frob_covered, na.rm=TRUE),
    mean_pve_width = median(pve_ci_width, na.rm=TRUE),
    se = sqrt(coverage * (1 - coverage) / n())
  )
head(plot_df)

# Scatter plot of coverage
g <- ggplot(plot_df, aes(x = tested_k, y = coverage, color=method, lty=method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=coverage-1.96*se, ymax=coverage+1.96*se), width=.1) +
  geom_hline(yintercept=0.9, linetype="dashed",color = "gray", size=0.5) +
  labs(
    x = TeX(r"( Index $k \leq r(s\{X\})$ of target $PVE_k(X)$ )"),
    y = "Coverage",
    col = "",
    linetype = "",
  ) +
  scale_color_manual(
    values = c("Unselective" = hue_pal()(3)[1], "Selective" = hue_pal()(3)[3]),
    labels = c("Non-selective", "Selective")
    ) +
  scale_linetype_manual(
    values = c("Unselective" = "dashed", "Selective" = "solid"),
    labels = c("Non-selective", "Selective")
    ) +
  # scale_shape_manual(values = c("Unselective" = 16, "Selective" = 16)) +
  theme_bw() +
  theme(
    # legend.direction = "horizontal",
    # legend.position = "bottom",
    # legend.position = c()
    # legend.box = "vertical",
    legend.spacing.y = unit(0.05, 'in'),
    legend.margin=margin(t=0.05, unit='in'),
    legend.box.spacing=unit(0.1, "in"),
    # legend.box.just = "left"
    legend.just = c("right", "top"),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),
    legend.key.size=unit(2,"lines")
  )
show(g)
ggsave(paste0('./figures/Figure3-', selection_rule, "-c2=", c^2, "-alpha_frac=", signal_alpha_frac, '-coverage_vs_k.png'), width = 5.5, height = 2, unit = "in")
