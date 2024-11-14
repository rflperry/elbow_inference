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

# See Figure4_zg_confidence_interval_width.R

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

selection_rule <- "zg"
if (selection_rule == "zg") {
  plot_sigmas <- c(0.1, 0.2)
} else if (selection_rule == "elbow") {
  plot_sigmas <- c(0.3, 2)
}

theme_update(text = element_text(size=10, family="Times"))
color_map <- c("Choi et al. (2017)" = hue_pal()(3)[1], "Choi et al. (2017) [Bf]" = hue_pal()(3)[2], "Selective Inference" = hue_pal()(3)[3])
scale_map <- c("Choi et al. (2017)" = 16, "Choi et al. (2017) [Bf]" = 17, "Selective Inference" = 15)
linetype_map <- c( "Selective Inference" = "solid", "Choi et al. (2017)" = "dashed")

load(paste0("data/", selection_rule, "-c2=1-alpha_frac=0.5-ci_vs_sigma.RData"))

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
    pve_ci_upper = pmin(mapply(get_signal_squared_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2, 1),
    pve = (signal / frob_norm )^2,
    pve_mle = mle^2 / frob_mle^2, # (mle / frob_mle)^2,
    pve_hat = mle^2 / frob2_hat,
    # pve_median = median^2 / ( val_sum - val_k + median^2),
    pve_naive = val_k / val_sum,
    pve_ci_width = pve_ci_upper - pve_ci_lower,
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
  ) %>% subset(
    method != "Choi et al. (2017) [Bf.]"
  ) %>%
  subset(selection_r == 5) %>%
  subset(sigma %in% plot_sigmas) %>% 
  pivot_wider(
    id_cols = c(tested_k, rep, sigma, precision),
    names_from = method,
    values_from = pve_ci_width
  ) %>%
  group_by(
    sigma, precision, tested_k
  ) %>%
  summarise(
    mean_pve_ci_diff = mean(Selective / Unselective, na.rm=TRUE)
  )
head(plot_df)

# Scatter plot with different colors and shapes
g <- ggplot(
  plot_df ,
  aes(x = tested_k, y = mean_pve_ci_diff, col = as.factor(sigma)),
  ) +
  geom_point() +
  labs(
    x = "Index k",
    y = "Width ratio",
    col = unname(TeX(r"( $\sigma$ )"))
  ) +
  theme_bw()#  +
  # facet_wrap(~ sigma)
show(g)
ggsave(paste0('./figures/Figure-', selection_rule, '-width_ratio.png'), width = 5.5, height = 1.5, unit = "in")

