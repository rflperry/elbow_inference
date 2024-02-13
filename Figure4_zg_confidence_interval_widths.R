library(dplyr)
library(tidyr)

source("./scripts/functions/confidence_intervals.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")
source("./scripts/functions/estimation.R")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Run Simulation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

n = 50
p = 10
rank = 5
m = 1
eigen = TRUE
reps <- 1000
alpha <- 0.1

# sigmas_list <- c(0.1, 0.15, 0.2, 0.3, 0.4, 0.7, 1)
selection_rule = "elbow" # options include "zg", "elbow"
if (selection_rule == "zg") {
  sigmas_list <- c(0.1, 0.2, 0.4, 0.5, 0.7, 1)
  degree <- 0.4
} else if (selection_rule == "elbow") {
  sigmas_list <- c(0.1, 0.3, 1, 2)
  degree <- 1
}

results <- c()
rep <- 1

for (rep in rep:reps) {
  for (sigma in sigmas_list) {
    try({
      set.seed(rep)
      # Simulate data
      sim <- simulate_matrix(n, p, sigma, rank, m, degree=degree, eigen=eigen)
      
      # Extract quantities
      duv <- svd(sim$obsv_mat)
      
      if (eigen) {
        vals <- duv$d^2
      } else {
        vals <- duv$d
      }
      
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
        base_results <- c(rep, n, p, sigma, m, rank, r, k, tr_mean, vals[k], sum(vals))
        
        # Choi p-value
        ci <- conf_interval_solver(
          choi_test_pvalue, vals, k, n, p, sigma=sigma, alpha=alpha, eigen=eigen
        )
        
        choi_mle <- NaN
        try({
          choi_mle <- get_mle(
            vals, k, n, p, sigma, eigen=eigen
          )
        })
        
        choi_median <- get_median(
          choi_test_pvalue, vals, k, n, p, sigma=sigma, eigen=eigen
        )
        results <- c(results, 'Choi', base_results, ci, choi_mle, choi_median, mean(ci))
        
        # Choi p-value, w/ Bonferroni correction
        ci <- conf_interval_solver(
          choi_test_pvalue, vals, k, n, p, sigma=sigma, alpha=alpha / (p-2), eigen=eigen
        )
        results <- c(results, 'Choi (Bf.)', base_results, ci, choi_mle, choi_median, mean(ci))
        
        # Selective inference p-value (R >= k)
        # Acquire bounds
        if (selection_rule == "zg") {
          bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)
        } else if (selection_rule == "elbow") {
          bounds <- compute_elbow_bounds(vals, k)
        }
        
        si_mle <- NaN
        try({
          si_mle <- get_mle(
            vals, k, n, p, sigma, bounds=bounds, eigen=eigen
          )
        })
        si_median <- get_median(
          si_test_pvalue, vals, k, n, p, sigma=sigma, bounds=bounds, eigen=eigen
        )
        
        ci <- conf_interval_solver(
          si_test_pvalue, vals, k, n, p, sigma=sigma, alpha=alpha, bounds=bounds, eigen=eigen
        )
        results <- c(results, 'Selective inference', base_results, ci, si_mle, si_median, mean(ci))
        
        # Selective inference p-value (R = r)
        # bounds <- compute_elbow_bounds(vals, k, inequality=FALSE)
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
  'rank', 'selection_r', 'tested_k', 'true_pivot', 'val_k', 'val_sum', 'ci_lower', 'ci_upper', 'mle', 'median', 'midpoint')
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
head(results_df)
print(nrow(results_df[complete.cases(results_df),]) / nrow(results_df))

save(results_df, file=paste0("data/", selection_rule, "-ci_vs_sigma.RData"))

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

theme_update(text = element_text(size=10, family="Times"))
color_map <- c("Choi et al. (2017)" = hue_pal()(3)[1], "Choi et al. (2017) [Bf]" = hue_pal()(3)[2], "Selective Inference" = hue_pal()(3)[3])
scale_map <- c("Choi et al. (2017)" = 16, "Choi et al. (2017) [Bf]" = 17, "Selective Inference" = 15)
linetype_map <- c( "Selective Inference" = "solid", "Choi et al. (2017)" = "dashed")

# load(paste0("data/", selection_rule, "-ci_vs_sigma.RData"))

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[1, ],
    pve_ci_upper = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[2, ],
    pve_mle = mle^2 / ( val_sum - val_k + mle^2),
    pve_median = median^2 / ( val_sum - val_k + median^2),
    pve = true_pivot^2 / ( val_sum - val_k + true_pivot^2),
    pve_naive = val_k / val_sum
  )
head(results_df)

plot_df <- results_df %>%
  # drop_na() %>%
  mutate(
    method = case_when(
      method == "Selective inference" ~ "Selective", # "Selective Inference",
      method == "Choi" ~ "Unselective", # "Choi et al. (2017)",
      method == "Choi (bf.)" ~ "Choi et al. (2017) [Bf.]",
    ),
    precision = 1/sigma,
    pve_covered = as.numeric((pve_ci_lower <= pve) & (pve <= pve_ci_upper)),
    pve_ci_width = pve_ci_upper - pve_ci_lower,
  ) %>% subset(
    method != "Choi et al. (2017) [Bf.]"
  ) %>%
  subset(selection_r == 5) %>%
  group_by(
    method, sigma, precision, tested_k
  ) %>%
  summarise(
    mean_pve_ci_lower = median(pve_ci_lower, na.rm=TRUE),
    mean_pve = median(pve, na.rm=TRUE),
    mean_pve_naive = median(pve_naive, na.rm=TRUE),
    mean_pve_ci_upper = median(pve_ci_upper, na.rm=TRUE),
    coverage = mean(pve_covered, na.rm=TRUE),
    mean_pve_width = median(pve_ci_width, na.rm=TRUE),
    # pve_mle_error = mean((pve_mle - pve)),
    # pve_median_error = mean((pve_median - pve)),
    mean_pve_mle = median(pve_mle, na.rm=TRUE),
    mean_pve_median = median(pve_median, na.rm=TRUE),
    # # midpoint_error = mean((midpoint - true_pivot)^2)
    # mle_mse = mean((pve_mle -pve)^2),
    # median_mse = mean((pve_median - pve)^2),
    # # midpoint_ = mean((midpoint - true_pivot)^2)
  )
head(plot_df)

# Scatter plot with different colors and shapes
g <- ggplot(plot_df %>% subset(sigma %in% c(0.3, 1)), aes(x = tested_k, y = mean_pve)) +
  geom_point(aes(
    x = tested_k, y = mean_pve_naive,
    color = "Naive",
    shape = "Naive")) +
  geom_point(aes(color = "Estimand", shape = "Estimand")) +
  geom_point(aes(
    y=mean_pve_mle, color=method, group=method, shape=method
  ),  position = position_dodge(0.6)) +
  geom_errorbar(aes(
    y=mean_pve_mle, ymin=mean_pve_ci_lower, ymax=mean_pve_ci_upper, color=method, group=method, lty=method,
  ), width=0.5, position = position_dodge(0.6)) +
  labs(
    x = "Index k",
    y = "PVE",
    col = "",
    shape = "",
    lty = "",
  ) +
  scale_color_manual(
    values = c("Estimand" = hue_pal()(3)[2], "Unselective" = hue_pal()(3)[1], "Selective" = hue_pal()(3)[3], "Naive" = "black"),
    labels = c(
      unname(TeX(r"( True $PVE^{\Theta}_k(X)$ )")),
      unname(TeX(r"( Non-selective $PVE^{\Theta}_k(X)$ )")),
      unname(TeX(r"( Selective $PVE^{\Theta}_k(X)$ )")),
      unname(TeX(r"( Naive $\widehat{PVE}_k(X)$ )")))
    ) +
  scale_shape_manual(
    values = c("Estimand" = 17, "Unselective" = 17, "Selective" = 17, "Naive" = 16),
    labels = c(
      unname(TeX(r"( True $PVE^{\Theta}_k(X)$ )")),
      unname(TeX(r"( Non-selective $PVE^{\Theta}_k(X)$ )")),
      unname(TeX(r"( Selective $PVE^{\Theta}_k(X)$ )")),
      unname(TeX(r"( Naive $\widehat{PVE}_k(X)$ )")))
    ) +
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
    legend.text = element_text(size = 8)
  ) +
  facet_wrap(~ sigma)
show(g)
ggsave(paste0('./figures_FINAL/Figure4-', selection_rule, '-coverage_screeplot.png'), width = 5.5, height = 2, unit = "in")

