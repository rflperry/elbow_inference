library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis) 
library(stringr)
library(scales)
library(latex2exp)
source("./scripts/functions/confidence_intervals.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/estimation.R")

theme_update(text = element_text(size=5, family="Times"))
color_map <- c("Choi et al. (2017)" = hue_pal()(3)[1], "Choi et al. (2017) [Bf]" = hue_pal()(3)[2], "Selective Inference" = hue_pal()(3)[3])
scale_map <- c("Choi et al. (2017)" = 16, "Choi et al. (2017) [Bf]" = 17, "Selective Inference" = 15)

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Toy Example
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

####### scree plot w/ CIs


n <- 50
p <- 10
eigen <- TRUE

# vals <- c(28, 23, 20, 15, 11, 6.3, 5.5, 3.5, 2.9, 2.1)
# vals[6:10] <- vals[6:10] * 0
# vals[1:5] <- vals[1:5] 

# i <- 1
# 
# i <- i + 1

# i <- 17 for sigma <- 0.2 looks nice
# i <- 19 for sigma <- 0.1, 2x miscoverave vs 1x
# i <- 45, sigma # 0.2
i <- 45
sigma <- 0.17
set.seed(i)
vals <- c(5, 4, 3, 2, 1, 0, 0, 0, 0, 0)
vals[1:5] <- (vals[1:5] + 1)^1.2

UV <- array(rnorm(n*p, sd=1), dim=c(n, p))
duv <- svd(UV)
U <- duv$u
V <- duv$v
mean_mat <- U %*% diag(sqrt(vals)) %*% t(V)
results <- {}
results$mean_mat <- mean_mat
results$obsv_mat <- mean_mat + array(rnorm(n*p, sd=sigma), dim=c(n, p))

# Extract quantities
duv <- svd(results$obsv_mat)
vals <- duv$d^2

# Perform selection
r <- select_r_zg(vals)

# Compute population terms
toy_results <- c()
alpha <- 0.1
for (k in seq(p)) {
  # Theta
  theta <- (t(duv$u[,k]) %*% mean_mat %*% duv$v[,k])[1]

  # PVE
  pve <- theta^2 / ( sum(vals) - vals[k] + theta^2)
  
  base_results <- c(k, theta, pve, vals[k] / sum(vals))
  
  # Confidence intervals on theta, PVE
  if (k <= r) {
    ci <- conf_interval_solver(
      choi_test_pvalue, vals, k, n, p, sigma=sigma, alpha=alpha, eigen=TRUE
    )
    pve_ci <- get_signal_squared_ci(ci[1], ci[2], vals[k], sum(vals))
    mle <- NaN
    try({
      mle <- get_mle(
        vals, k, n, p, sigma, eigen=eigen
      )
    })
    pve_mle <- mle^2 / ( sum(vals) - vals[k] + mle^2)
    
    toy_results <- c(toy_results, base_results, "Unselective", ci, pve_ci, mle, pve_mle)
    
    # Selective inference p-value (R >= k)
    bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)
    
    ci <- conf_interval_solver(
      si_test_pvalue, vals, k, n, p, sigma=sigma, alpha=alpha, bounds=bounds, eigen=TRUE
    )
    pve_ci <- get_signal_squared_ci(ci[1], ci[2], vals[k], sum(vals))
    
    mle <- NaN
    try({
      mle <- get_mle(
        vals, k, n, p, sigma, bounds=bounds, eigen=eigen
      )
    })
    pve_mle <- mle^2 / ( sum(vals) - vals[k] + mle^2)

    toy_results <- c(toy_results, base_results, "Selective", ci, pve_ci, mle, pve_mle)
  } else {
    toy_results <- c(toy_results, c(base_results, NA, NA, NA, NA, NA, NA, NA))
  }
}

header <- c(
  "Index", "theta", "PVE", "PVE_naive", "Method", 'ci_lower', 'ci_upper', 'pve_ci_lower', 'pve_ci_upper', 'mle', 'pve_mle')
results_df <- data.frame(
  data=t(array(toy_results, dim = c(length(header), length(toy_results)/length(header)))))
colnames(results_df) <- header
results_df[,1:4] <- sapply(results_df[,1:4], as.numeric)
results_df[,6:length(header)] <- sapply(results_df[,6:length(header)], as.numeric)
results_df$Index <- factor(results_df$Index)
head(results_df)

# Scatter plot with different colors and shapes
g <- ggplot(results_df, aes(x = Index, y = PVE)) +
  geom_point(aes(color = "Estimand", shape = "Estimand")) +
  geom_point(aes(
    x = Index, y = PVE_naive,
    color = "Naive",
    shape = "Naive")) +
  geom_point(data=results_df %>% drop_na(), aes(
    x=Index, y=pve_mle, color=Method, group=Method, shape = Method
  ),  position = position_dodge(0.6)) +
  geom_errorbar(data=results_df %>% drop_na(), aes(
    x=Index, y=pve_mle, ymin=pve_ci_lower, ymax=pve_ci_upper, color=Method, group=Method, lty=Method,
      ), width=0.5, position = position_dodge(0.6)) +
  labs(
       x = "Index K",
       y = "PVE",
       col = "",
       shape = "",
       lty = "",
       ) +
  scale_color_manual(values = c("Estimand" = hue_pal()(3)[2], "Naive" = "Black", "Unselective" = hue_pal()(3)[1], "Selective" = hue_pal()(3)[3])) +
  scale_shape_manual(values = c("Estimand" = 17, "Naive" = 16, "Unselective" = 17, "Selective" = 17)) +
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
  )
show(g)
ggsave(paste0('./figures_FINAL/Figure1_screeplot_zg.png'), width = 3.7, height = 2, unit = "in")
