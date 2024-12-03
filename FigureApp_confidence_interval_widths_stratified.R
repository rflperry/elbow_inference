library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(scales)
library(latex2exp)
library(argparse)

source("./scripts/functions/confidence_intervals.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")
source("./scripts/functions/estimation.R")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------

args <- list()
args$input_file <- "data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=1000_signal_alpha_frac=0.75.RData"

parser <- ArgumentParser()
parser$add_argument("input_file", nargs = 1, help = "File to be displayed")
parser$add_argument("--sigmas",
  type = "double",
  nargs = "+",
  default = c(0.1, 0.2)
)
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

fname <- paste0("figures/FigureApp_", sub("\\.RData$", "_stratified.png", basename(args$input_file)))

# # selection_rule <- "elbow"
# if (selection_rule == "zg") {
#   plot_sigmas <- c(0.1, 0.2) # c(0.2, 0.4, 0.7, 1)
# } else if (selection_rule == "elbow") {
#   plot_sigmas <- sigmas # c(0.3, 1)
# }

theme_update(text = element_text(size = 10, family = "Times"))

load(args$input_file)

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
    pve_ci_upper = pmin(mapply(get_signal_squared_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2, 1),
    pve = (signal / frob_norm)^2,
    pve_mle = mle^2 / frob_mle^2,
    pve_hat = pmin(mle^2 / frob2_hat, 1),
    pve_naive = val_k / val_sum
  )
head(results_df)

summary_func <- median

plot_df <- results_df %>%
  # drop_na() %>%
  mutate(
    method = case_when(
      method == "Selective inference" ~ "Selective", # "Selective Inference",
      method == "Choi" ~ "Unselective", # "Choi et al. (2017)",
      method == "Choi (bf.)" ~ "Choi et al. (2017) [Bf.]",
    ),
    precision = 1 / sigma^2,
    pve_covered = as.numeric((pve_ci_lower <= pve) & (pve <= pve_ci_upper)),
    pve_ci_width = pve_ci_upper - pve_ci_lower,
  ) %>%
  subset(
    method == "Selective" &
      sigma %in% args$sigmas
  ) %>%
  group_by(
    sigma, method, tested_k, selection_r
  ) %>%
  summarise(
    mean_pve_ci_lower = summary_func(pve_ci_lower, na.rm = TRUE),
    mean_pve = summary_func(pve, na.rm = TRUE),
    mean_pve_naive = summary_func(pve_naive, na.rm = TRUE),
    mean_pve_ci_upper = summary_func(pve_ci_upper, na.rm = TRUE),
    coverage = mean(pve_covered, na.rm = TRUE),
    mean_pve_width = summary_func(pve_ci_width, na.rm = TRUE),
    mean_pve_mle = summary_func(pve_mle, na.rm = TRUE),
    mean_pve_hat = summary_func(pve_hat, na.rm = TRUE),
    count = n()
  )
head(plot_df)

# Scatter plot with different colors and shapes
g <- ggplot(plot_df, aes(x = tested_k)) +
  geom_point(aes(
    y = mean_pve_hat, color = "Selective", shape = "Selective", group = "1"
  )) +
  geom_point(aes(
    y = mean_pve, color = "Estimand", shape = "Estimand", group = "1"
  )) +
  geom_point(aes(
    y = mean_pve_naive, color = "Naive", shape = "Naive", group = "1"
  )) +
  geom_errorbar(aes(
    y = mean_pve_hat, ymin = mean_pve_ci_lower, ymax = mean_pve_ci_upper,
    color = "Selective", lty = "Selective", group = "2"
  ), width = 0.5) +
  labs(
    x = "Index k",
    y = "Median PVE quantities",
    col = "",
    shape = "",
    lty = "",
  ) +
  scale_color_manual(
    values = c("Estimand" = hue_pal()(3)[2], "Selective" = hue_pal()(3)[3], "Naive" = "black"),
    labels = c(
      TeX(r"( $PVE_k\{X^{(1)}\}$ )"),
      TeX(r"( $PVE_{k, MLE}\{X^{(1)}\}$ )"),
      TeX(r"( $\widehat{PVE}_k\{X^{(1)}\}$ )")
    )) +
  scale_linetype_manual(
    values = c(
      # "Estimand" = "blank",
      # "Naive" = "blank",
      "Selective" = "solid"
    ),
    labels = c(
      # "",
      # "",
      TeX(r"( CI for $PVE_k\{X^{(1)}\}$ )")
    )
  ) +
  scale_shape_manual(
    values = c("Estimand" = 4, "Selective" = 17, "Naive" = 16),
    labels = c(
      TeX(r"( $PVE_k\{X^{(1)}\}$ )"),
      TeX(r"( $PVE_{k, MLE}\{X^{(1)}\}$ )"),
      TeX(r"( $\widehat{PVE}_k\{X^{(1)}\}$ )")
    )
  ) +
  # guides(
  #   linetype = guide_legend(
  #     # nrow = 1,
  #     # byrow = TRUE,
  #     override.aes = list(
  #       shape = c(NA, NA, 16),
  #       linetype = c("solid", "solid", "blank")
  #       )
  #   ),
  #   shape = guide_legend(
  #     nrow = 3,
  #     byrow = TRUE,
  #     override.aes = list(
  #       linetype = c("blank", "blank", "blank"),
  #       color = c("Estimand" = hue_pal()(3)[2], "Selective" = hue_pal()(3)[3], "Naive" = "black")
  #     )
  #   )
  # ) +
  theme_bw() +
  theme(
    legend.spacing.y = unit(0.05, "in"),
    legend.margin = margin(t = 0.05, unit = "in"),
    legend.box.spacing = unit(0.1, "in"),
    legend.just = c("right", "top"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  facet_grid(
    sigma~selection_r,
    scales = "free_y",
    labeller = labeller(
      selection_r = function(x) paste0("r(s(X)) = ", x),
      sigma = function(x) paste0(TeX(r"( sigma = )"), x)
      ))
show(g)
ggsave(fname, width = 5.5, height = 4, unit = "in")
