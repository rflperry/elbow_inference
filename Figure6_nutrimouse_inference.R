library('CCA')
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis) 
library(stringr)
library(scales)
library(RMTstat)
library(latex2exp)

source("./scripts/functions/hypothesis_tests.R")
source("./scripts/functions/confidence_intervals.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/estimation.R")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	PreliMinaries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

eigen=TRUE
alpha <- 0.1


#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Load Data and transform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

data(nutrimouse)

raw_data <- nutrimouse$gene[,1:20]

n <- nrow(raw_data)
p <- ncol(raw_data)

# Transform data

I <- diag(1, n)  # Identity matrix
J <- matrix(1/n, n, n)  # Matrix with all elements equal to 1/n
C <- I - J  # Centering matrix

duv <- svd(C)

H <- duv$u %*% sqrt(diag(duv$d))

print(H %*% t(H))

data <- t(H) %*% as.matrix(raw_data)

# data <- scale(data, center = TRUE, scale = FALSE)

SCALAR <- max(data) / 2 * sqrt(n)

data <- data / SCALAR

duv <- svd(data)
vals <- duv$d^2
sigma <- sqrt(median(sqrt(vals))^2 / (max(n, p) * qmp(0.5, svr=max(n, p)/min(n, p))))

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Data thin data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

set.seed(6)
noise <- array(rnorm(n*p, sd=sigma), dim=c(n, p))
c <- sqrt(1)
data1 <- data + c * noise
data2 <- data - noise / c

signal_alpha_frac <- 0.9

sigma1 <- sigma * sqrt(1 + c^2)
sigma2 <- sigma * sqrt(1 + 1/c^2)

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	SVD
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

# SVD
duv <- svd(data1)
vals <- duv$d^2

r <- select_r_zg(vals / sum(vals))
print(r)

# Plot evals

g <- ggplot(
  data=data.frame(x=seq(1,length(vals)), y=vals / sum(vals)),
  aes(x=x, y=y)) +
  geom_point() +
  theme(text=element_text(size=10)) +
  scale_x_continuous(breaks= pretty_breaks()) +
  labs(
    x = TeX(r"( Index $k$ )"),
    y = TeX(r"( $\widehat{PVE}_k(x)$ )"),
  ) +
  theme_bw()+
  scale_y_continuous(limits = c(0, NA))
print(g)
ggsave(paste0('figures/Figure6-nutrimouse_screeplot.png'), width = 2, height = 2, unit = "in")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Test
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

frob_hat <- sqrt(sum(data2^2))
frob2_consistent <- frob_hat^2 - n * p * sigma2^2
frob_ci <- sqrt( sigma2^2 * get_nc_chsq_ci(frob_hat^2 / sigma2^2, n*p, alpha * ( 1- signal_alpha_frac)) )

# Selective inference p-value (R >= k)
results <- c()

for (k in seq(1, r)) {
  bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)

  si_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma1, bounds=bounds, eigen=eigen)
  
  ci <- conf_interval_solver(
    si_test_pvalue, vals, k, n, p, sigma=sigma1, alpha=alpha * signal_alpha_frac, bounds=bounds, eigen=eigen
  )
  
  si_mle <- NaN
  try({
    si_mle <- get_mle_signal(
      vals, k, n, p, sigma1, bounds=bounds, eigen=eigen
    )
  })
  
  results <- c(results, c('Selective', k, si_pvalue, ci), sum(vals), vals[k], si_mle, frob2_consistent, frob_ci)
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	CI Plots
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
header <- c(
  "method", "k", "pvalue", 'ci_lower', 'ci_upper', "val_sum", "val_k", "mle", "frob2_hat", "frob_ci_lower", "frob_ci_upper")
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
head(results_df)

theme_update(text = element_text(size=10, family="Times"))

plot_df <- results_df %>% mutate(
  pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
  pve_ci_upper = pmin(mapply(get_pve_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2, 1),
  pve_hat = mle^2 / frob2_hat,
  pve_naive = val_k / val_sum
)
head(plot_df)

g <- ggplot(
  plot_df,
  aes(x=k, y=pve_hat, color=method, group=method)) +
  # geom_line(size=0.5) +
  geom_point(size=2) +
  geom_point(aes(y=pve_naive), size=2, color="black", group="Naive") +
  geom_errorbar(aes(ymin=pve_ci_lower, ymax=pve_ci_upper), width=0.5) +
  labs(
    x = TeX(r"( Index $k$ )"),
    y = TeX(r"( PVE )"),
    col = "",
    linetype = "",
  ) +
  # labs(color='') +
  # theme_bw() +
  # theme(
  #   legend.position="none") +
  # labs(col='Method') +
  scale_color_manual(
    values = c("Selective" = hue_pal()(3)[3], "Naive" = "black"),
    labels = c(unname(TeX(r"( Selective $PVE_k(x)$ )")), unname(TeX(r"( Naive $\widehat{PVE}_k(x)$ )")))
    ) +
  # scale_linetype_manual(values = c("Selective" = "solid", unname(TeX(r"( $\widehat{PVE}_k(x)$ )")) = "blank")) +
  # scale_shape_manual(values = scale_map) +
  scale_x_continuous(breaks= c(1, 2, 3), labels=c("1", "2", "3")) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(0, NA))
  # theme(legend.position="none")
plot(g)
ggsave(paste0('figures/Figure6-nutrimouse_confidence_intervals.png'), width = 4, height = 2, unit = "in")


#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Iterate process over noise levels
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#


for (nl in c(0, 1, 2, 3, 4, 5, 6, 7, 8)) {
  set.seed(1)
  noise_level <- nl * sigma
  sigma2 <- sqrt(noise_level^2 + sigma^2)
  duv <- svd(data + array(rnorm(n*p, sd=noise_level), dim=c(n, p)))
  vals <- duv$d^2
  r <- select_r_zg(vals / sum(vals))
  print(r)
  df <- data.frame(x=seq(1,length(vals)), y=vals / sum(vals), r=r, noise_level=noise_level, label="Raw data")
  
  if (nl > 0) {
    title <- paste0("Raw data + ", nl, "x noise")
  } else {
    title <- "Raw data"
  }
  g <- ggplot(
    data=data.frame(x=seq(1,length(vals)), y=vals / sum(vals)),
    aes(x=x, y=y)) +
    geom_point() +
    scale_x_continuous(breaks= pretty_breaks()) +
    labs(
      x = TeX(r"( Index $k$ )"),
      y = TeX(r"( $\widehat{PVE}_k(x)$ )"),
      title = title
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10),
      legend.title = element_text(size = 10), 
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(limits = c(0, NA))
  ggsave(paste0('figures_Final/Figure6-nutrimouse_screeplot_', nl, '.png'), width = 2, height = 1.5, unit = "in")
  
  results <- c()
  
  for (k in seq(1, r)) {
    bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)
    
    si_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma2, bounds=bounds, eigen=eigen)
    
    ci <- conf_interval_solver(
      si_test_pvalue, vals, k, n, p, sigma=sigma2, alpha=alpha * signal_alpha_frac, bounds=bounds, eigen=eigen
    )
    
    si_mle <- NaN
    try({
      si_mle <- get_mle_signal(
        vals, k, n, p, sigma2, bounds=bounds, eigen=eigen
      )
    })
    
    results <- c(results, c('Selective', k, si_pvalue, ci), sum(vals), vals[k], si_mle)
  }
  
  header <- c(
    "method", "k", "pvalue", 'ci_lower', 'ci_upper', "val_sum", "val_k", "mle")
  results_df <- data.frame(
    data=t(array(results, dim = c(length(header), length(results)/length(header)))))
  colnames(results_df) <- header
  results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
  head(results_df)
  
  theme_update(text = element_text(size=10, family="Times"))
  
  plot_df <- results_df %>% mutate(
    pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[1, ],
    pve_ci_upper = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[2, ],
    pve_mle = mle^2 / ( val_sum - val_k + mle^2),
    pve_naive = val_k / (val_sum)
  )
  head(plot_df)
  
  g <- ggplot(
    plot_df,
    aes(x=k, y=pve_mle, color=method, group=method)) +
    # geom_line(size=0.5) +
    geom_point(size=2) +
    geom_point(aes(y=pve_naive), size=2, color="black", group="Naive") +
    geom_errorbar(aes(ymin=pve_ci_lower, ymax=pve_ci_upper), width=0.5) +
    labs(
      x = TeX(r"( Index $k$ )"),
      y = TeX(r"( PVE )"),
      col = "",
      linetype = ""# ,
      # title = title
    ) +
    scale_color_manual(
      values = c("Selective" = hue_pal()(3)[3], "Naive" = "black"),
      labels = c(unname(TeX(r"( Selective $PVE_k(x)$ )")), unname(TeX(r"( Naive $\widehat{PVE}_k(x)$ )")))
    ) +
    scale_x_continuous(breaks=pretty_breaks()) + # breaks= c(1, 2, 3), labels=c("1", "2", "3")) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10),
      legend.title = element_text(size = 10), 
      legend.text = element_text(size = 10)
    ) +
    scale_y_continuous(limits = c(0, NA))
  ggsave(paste0('figures_FINAL/Figure6-nutrimouse_confidence_intervals_', nl, '.png'), width = 4, height = 1.5, unit = "in")
}

# #
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 	Facet grid
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #

results <- c()

set.seed(6)
dt_noise <- array(rnorm(n*p, sd=sigma), dim=c(n, p))

set.seed(1)
noise <- array(rnorm(n*p, sd=sigma), dim=c(n, p))

c <- sqrt(1)
signal_alpha_frac <- 0.9

for (nl in c(0, 1, 2, 3, 4, 5, 6, 7, 8)) {
  nl <- nl / 2
  noise_level <- nl * sigma
  
  data_noised <- data + nl * noise
  
  data_noised1 <- data_noised + sqrt(1 + nl^2) * dt_noise * c
  data_noised2 <- data_noised + sqrt(1 + nl^2) * dt_noise / c
  
  sigma1 <- sigma * sqrt(1 + nl^2) * sqrt(1 + c^2)
  sigma2 <- sigma * sqrt(1 + nl^2) * sqrt(1 + 1/c^2)
  
  duv <- svd(data_noised1)
  vals <- duv$d^2
  r <- select_r_zg(vals / sum(vals))
  print(r)
  df <- data.frame(x=seq(1,length(vals)), y=vals / sum(vals), r=r, noise_level=noise_level, label="Raw data")
  
  frob_hat <- sqrt(sum(data_noised2^2))
  frob2_consistent <- frob_hat^2 - n * p * sigma2^2
  frob_ci <- sqrt( sigma2^2 * get_nc_chsq_ci(frob_hat^2 / sigma2^2, n*p, alpha * ( 1- signal_alpha_frac)) )

  for (k in seq(1, r)) {
    bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)
    
    naive_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma1, eigen=eigen)
    si_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma1, bounds=bounds, eigen=eigen)
    
    naive_ci <- conf_interval_solver(
      si_test_pvalue, vals, k, n, p, sigma=sigma1, alpha=alpha * signal_alpha_frac, eigen=eigen
    )
    
    ci <- conf_interval_solver(
      si_test_pvalue, vals, k, n, p, sigma=sigma1, alpha=alpha * signal_alpha_frac, bounds=bounds, eigen=eigen
    )
    
    si_mle <- NaN
    try({
      si_mle <- get_mle_signal(
        vals, k, n, p, sigma1, bounds=bounds, eigen=eigen
      )
    })
    
    naive_mle <- NaN
    try({
      naive_mle <- get_mle_signal(
        vals, k, n, p, sigma1, eigen=eigen
      )
    })
    
    lb <- vals[k+1]
    ub <- c(Inf, vals)[k]
    
    results <- c(results, c('Selective', k, si_pvalue, ci), sum(vals), vals[k], si_mle, bounds[[1]], nl, frob2_consistent, frob_ci)
    results <- c(results, c('Non-selective', k, naive_pvalue, naive_ci), sum(vals), vals[k], naive_mle, lb, ub, nl, frob2_consistent, frob_ci)
  }
}

header <- c(
  "method", "k", "pvalue", 'ci_lower', 'ci_upper', "val_sum", "val_k", "mle", "lb", "ub", "noise_level", "frob2_hat", "frob_ci_lower", "frob_ci_upper")
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
head(results_df)

plot_df <- results_df %>% mutate(
  pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
  pve_ci_upper = pmin(mapply(get_pve_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2, 1),
  pve_mle = mle^2 / frob2_hat,
  pve_naive = val_k / val_sum,
  int_width = ub - lb,
  pve_ci_width = pve_ci_upper - pve_ci_lower,
  ci_width = ci_upper - ci_lower,
  facet_label = paste0("c = ", noise_level)
)
tail(plot_df)

g <- ggplot(
  plot_df %>% subset(noise_level < 3),
  aes(x=k, y=pve_mle, color=method, group=method)) +
  geom_point(size=1, position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin=pve_ci_lower, ymax=pve_ci_upper), width=0.5, position = position_dodge(0.6)) +
  labs(
    x = TeX(r"( Index $k$ )"),
    y = TeX(r"( PVE )"),
    col = "",
    linetype = "",
  ) +
  scale_color_manual(
    values = c("Selective" = hue_pal()(3)[3], "Non-selective" = hue_pal()(3)[2]),
    labels = c(unname(TeX(r"( Selective $PVE_k(x)$ )")), unname(TeX(r"( Non-selective $PVE_k(x)$ )")))
  ) +
  scale_x_continuous(breaks=pretty_breaks()) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 10),
    legend.position="bottom"
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ facet_label)
plot(g)
ggsave(paste0('figures/Figure6-nutrimouse_confidence_intervals_appendix_pve.png'), width = 6, height = 3, unit = "in")

g <- ggplot(
  plot_df  %>% subset(noise_level < 3),
  aes(x=k, y=mle, color=method, group=method)) +
  geom_point(size=1, position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=0.5, position = position_dodge(0.6)) +
  labs(
    x = TeX(r"( Index $k$ )"),
    y = TeX(r"( $u_k\{x\}^T \Theta  v_k\{x\}$ )"),
    col = "",
    linetype = "",
  ) +
  scale_color_manual(
    values = c("Selective" = hue_pal()(3)[3], "Non-selective" = hue_pal()(3)[2]),
    labels = c(unname(TeX(r"( Selective $u_k\{x\}^T \Theta  v_k\{x\}$ )")), unname(TeX(r"( Non-selective $u_k\{x\}^T \Theta  v_k\{x\}$ )")))
  ) +
  scale_x_continuous(breaks=pretty_breaks()) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 10),
    legend.position="bottom"
  ) +
  facet_wrap(~ facet_label)
plot(g)
ggsave(paste0('figures/Figure6-nutrimouse_confidence_intervals_appendix_theta.png'), width = 6, height = 4, unit = "in")

g <- ggplot(
  plot_df, # %>% subset(noise_level == 8),# %in% c(0, 2, 4, 8)),
  aes(x=val_k - lb, y=ci_lower, shape=method, color=method)) +
  # geom_line(size=0.5) +
  geom_point(size=2) +
  # geom_point(aes(y=pve_naive), size=2, color="black", group="Naive") +
  # geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper, color=method, group=method), width=0.5, position = position_dodge(0.6)) +
  # scale_color_manual(
  #   values = c("Selective" = hue_pal()(3)[3], "Non-selective" = hue_pal()(3)[2], "Naive" = "black"),
  #   labels = c(unname(TeX(r"( Selective $PVE_k(x)$ )")), unname(TeX(r"( Non-selective $PVE_k(x)$ )")), unname(TeX(r"( Naive $\widehat{PVE}_k(x)$ )")))
  # ) +
  scale_x_continuous(breaks=pretty_breaks()) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 10)
  ) # +
  # scale_y_continuous(limits = c(0, NA)) # +
  # facet_wrap(~ noise_level)
plot(g)


# #
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 	Iterate a lot
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #

results <- c()

for (rep in seq(1, 20)) {
  set.seed(rep)
  nl <- 8
  noise_level <- nl * sigma
  sigma2 <- sqrt(noise_level^2 + sigma^2)
  duv <- svd(data + array(rnorm(n*p, sd=noise_level), dim=c(n, p)))
  vals <- duv$d^2
  r <- select_r_zg(vals / sum(vals))
  print(r)
  df <- data.frame(x=seq(1,length(vals)), y=vals / sum(vals), r=r, noise_level=noise_level, label="Raw data")
  
  for (k in seq(1, r)) {
    bounds <- compute_constraint_regions(vals, k, select_r_zg, tol=1000)
    
    naive_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma2, eigen=eigen)
    si_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma2, bounds=bounds, eigen=eigen)
    
    naive_ci <- conf_interval_solver(
      si_test_pvalue, vals, k, n, p, sigma=sigma2, alpha=alpha, eigen=eigen
    )
    
    ci <- conf_interval_solver(
      si_test_pvalue, vals, k, n, p, sigma=sigma2, alpha=alpha, bounds=bounds, eigen=eigen
    )
    
    si_mle <- NaN
    try({
      si_mle <- get_mle(
        vals, k, n, p, sigma2, bounds=bounds, eigen=eigen
      )
    })
    
    lb <- vals[k+1]
    ub <- c(Inf, vals)[k]
    
    results <- c(results, c('Selective', k, si_pvalue, ci), sum(vals), vals[k], si_mle, bounds[[1]], nl)
    results <- c(results, c('Non-selective', k, naive_pvalue, naive_ci), sum(vals), vals[k], si_mle, lb, ub, nl)
  }
}

header <- c(
  "method", "k", "pvalue", 'ci_lower', 'ci_upper', "val_sum", "val_k", "mle", "lb", "ub", "noise_level")
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
head(results_df)

plot_df <- results_df %>% mutate(
  pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[1, ],
  pve_ci_upper = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[2, ],
  pve_mle = mle^2 / ( val_sum - val_k + mle^2),
  pve_naive = val_k / (val_sum),
  int_width = ub - lb,
  pve_ci_width = pve_ci_upper - pve_ci_lower,
  ci_width = ci_upper - ci_lower
)

g <- ggplot(
  plot_df, # %>% subset(noise_level == 8),# %in% c(0, 2, 4, 8)),
  aes(x=val_k - lb, y=ci_lower, shape=method, color=method)) +
  # geom_line(size=0.5) +
  geom_point(size=2) +
  # geom_line() +
  # geom_point(aes(y=pve_naive), size=2, color="black", group="Naive") +
  # geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper, color=method, group=method), width=0.5, position = position_dodge(0.6)) +
  # scale_color_manual(
  #   values = c("Selective" = hue_pal()(3)[3], "Non-selective" = hue_pal()(3)[2], "Naive" = "black"),
  #   labels = c(unname(TeX(r"( Selective $PVE_k(x)$ )")), unname(TeX(r"( Non-selective $PVE_k(x)$ )")), unname(TeX(r"( Naive $\widehat{PVE}_k(x)$ )")))
  # ) +
  scale_x_continuous(breaks=pretty_breaks()) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 10)
  ) +
  facet_wrap(~ (k))
# scale_y_continuous(limits = c(0, NA)) # +
# facet_wrap(~ noise_level)
plot(g)
