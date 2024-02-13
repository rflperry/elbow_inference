library(dplyr)
library(tidyr)

source("./scripts/functions/hypothesis_tests.R")
source("./scripts/functions/randomized_test.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")
source("./scripts/functions/confidence_intervals.R")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Run Simulation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

# See Figure4_zg_confidence_interval_width.R

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Process data and save
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#


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

selection_rule = "elbow"

theme_update(text = element_text(size=10, family="Times"))

load(paste0("data/", selection_rule, "-ci_vs_sigma.RData"))

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    method = case_when(
      method == "Selective inference" ~ "Selective", # "Selective Inference",
      method == "Choi" ~ "Unselective", # "Choi et al. (2017)",
      method == "Choi (bf.)" ~ "Choi et al. (2017) [Bf.]",
    ),
    sigma <- sigma^2,
    pve_ci_lower = mapply(get_pve_ci, ci_lower, ci_upper, val_k, val_sum)[1, ],
    precision = 1/sigma,
    reject_null = as.numeric((pve_ci_lower > 0)),
  ) %>% subset(
    method != "Choi et al. (2017) [Bf.]"
  )
head(results_df)

plot_df <- results_df %>%
  # drop_na() %>%
  subset(selection_r == 5) %>%
  group_by(
    method, sigma, precision, tested_k
  ) %>%
  summarise(
    power = mean(reject_null, na.rm=TRUE)
  )
head(plot_df)

######### plot detection probability
g_dp <- ggplot(
  results_df %>% group_by(precision) %>% summarise(detection_probability = mean(selection_r == rank)),
  aes(x=precision, y=detection_probability)) +
  geom_line() +
  geom_point(size=1) +
  xlab('1 / Variance') +
  ylab('Detection probability') +
  scale_x_log10() +
  theme(text=element_text(size=10)) +
  theme_bw()
plot(g_dp)
ggsave(paste0('figures_FINAL/Figure5-', selection_rule, '-detection_probability_vs_sigma.png'), width = 2, height = 2, unit = "in")

# Scatter plot with different colors and shapes
g <- ggplot(plot_df %>% subset(method == 'Selective'), aes(x = precision, y = power, color = as.factor(tested_k))) +
  geom_point() +
  geom_line() +
  labs(
    x = "1 / Variance",
    y = "Conditional power",
    col = "Tested index k"# TeX(r"( Tested index $k \leq r(s\{X\})$ of target $PVE^{\theta}_k(X)$ )"),
  ) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  # scale_color_manual(values = c("Estimand" = hue_pal()(3)[2], "Unselective" = hue_pal()(3)[1], "Selective" = hue_pal()(3)[3])) +
  # scale_shape_manual(values = c("Estimand" = 17, "Unselective" = 17, "Selective" = 17)) +
  theme_bw() +
  scale_x_log10() +
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
ggsave(paste0('./figures_FINAL/Figure5-', selection_rule, '-power_vs_sigma.png'), width = 4, height = 2, unit = "in")

