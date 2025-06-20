library(CCA)
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
# 	Preliminaries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

eigen <- TRUE
alpha <- 0.1

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Load Data and transform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

data(nutrimouse)

raw_data <- nutrimouse$gene[, 1:20]

n <- nrow(raw_data)
p <- ncol(raw_data)

# Transform data

I <- diag(1, n) # Identity matrix
J <- matrix(1 / n, n, n) # Matrix with all elements equal to 1/n
C <- I - J # Centering matrix

duv <- svd(C)

H <- duv$u %*% sqrt(diag(duv$d))

print(H %*% t(H))

data <- t(H) %*% as.matrix(raw_data)

# data <- scale(data, center = TRUE, scale = FALSE)

SCALAR <- max(data) / 2 * sqrt(n)

data <- data / SCALAR

duv <- svd(data)
vals <- duv$d^2
sigma <- sqrt(median(sqrt(vals))^2 / (max(n, p) * qmp(0.5, svr = max(n, p) / min(n, p))))

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Data thin data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

set.seed(6)
noise <- array(rnorm(n * p, sd = sigma), dim = c(n, p))
# saveRDS(noise, file = "data/nutrimouse_noise_array.rds")
# # Uncomment out the line below if set.seed is not working for reproducibility purposes
# noise <- readRDS("data/nutrimouse_noise_array.rds")

c <- sqrt(1)
data1 <- data + c * noise
data2 <- data - noise / c

signal_alpha_frac <- 0.75

sigma1 <- sigma * sqrt(1 + c^2)
sigma2 <- sigma * sqrt(1 + 1 / c^2)

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
  data = data.frame(x = seq(1, length(vals)), y = vals / sum(vals)),
  aes(x = x, y = y)
) +
  geom_point() +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    x = TeX(r"( Index $k$ )"),
    y = TeX(r"( $\widehat{PVE}_k(X^{(1)})$ )"),
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA))
print(g)
ggsave(paste0("figures/Figure7-nutrimouse_screeplot.png"), width = 2, height = 2, unit = "in")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Test
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

frob_hat <- sqrt(sum(data2^2))
frob2_consistent <- frob_hat^2 - n * p * sigma2^2
frob_ci <- sqrt(sigma2^2 * get_nc_chsq_ci(frob_hat^2 / sigma2^2, n * p, alpha * (1 - signal_alpha_frac)))

# Selective inference p-value (R >= k)
results <- c()

for (k in seq(1, r)) {
  bounds <- compute_constraint_regions(vals, k, select_r_zg, tol = 1000)

  si_pvalue <- si_test_pvalue(vals, k = k, n = n, p = p, sigma = sigma1, bounds = bounds, eigen = eigen)

  ci <- conf_interval_solver(
    si_test_pvalue, vals, k, n, p,
    sigma = sigma1, alpha = alpha * signal_alpha_frac, bounds = bounds, eigen = eigen
  )

  si_mle <- NaN
  try({
    si_mle <- get_mle_signal(
      vals, k, n, p, sigma1,
      bounds = bounds, eigen = eigen
    )
  })

  results <- c(results, c("Selective", k, si_pvalue, ci), sum(vals), vals[k], si_mle, frob2_consistent, frob_ci)
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	CI Plots
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
header <- c(
  "method", "k", "pvalue", "ci_lower", "ci_upper", "val_sum", "val_k", "mle", "frob2_hat", "frob_ci_lower", "frob_ci_upper"
)
results_df <- data.frame(
  data = t(array(results, dim = c(length(header), length(results) / length(header))))
)
colnames(results_df) <- header
results_df[, 2:length(header)] <- sapply(results_df[, 2:length(header)], as.numeric)
head(results_df)

theme_update(text = element_text(size = 10, family = "Times"))
labels <- c(unname(TeX(r"( Selective $PVE_k$ )")), unname(TeX(r"( Naive $\widehat{PVE}_k$ )")))

plot_df <- results_df %>% mutate(
  pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
  pve_ci_upper = pmin(mapply(get_signal_squared_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2, 1),
  pve_hat = mle^2 / frob2_hat,
  pve_naive = val_k / val_sum
)
head(plot_df)

g <- ggplot(
  plot_df,
  aes(x = k, y = pve_hat, color = method, group = method)
) +
  geom_point(size = 2, group = "Selective", shape = 17) +
  geom_errorbar(aes(ymin = pve_ci_lower, ymax = pve_ci_upper), width = 0.5) +
  geom_point(aes(y = pve_naive), size = 2, shape = 16, color = "black", group = "Naive") +
  labs(
    x = TeX(r"( Index $k$ )"),
    y = TeX(r"( PVE )"),
    col = ""
  ) +
  scale_color_manual(
    values = c("Selective" = hue_pal()(3)[3], "Naive" = "black"),
    labels = labels,
  ) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("1", "2", "3")) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(0, NA))
plot(g)
ggsave(paste0("figures/Figure7-nutrimouse_confidence_intervals.png"), width = 4, height = 2, unit = "in")