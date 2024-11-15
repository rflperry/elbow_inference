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

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("input_file", nargs = 1, help = "File to be displayed")
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

args <- list()
args$input_file <- "data/sim_conf_ints_alpha=0.1_c=1_m=1_method=zg_n=50_p=10_rank=5_reps=1000_signal_alpha_frac=0.75.RData"

fname <- paste0("figures/Figure_", sub("\\.RData$", "-width_ratio.png", basename(args$input_file)))
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

theme_update(text = element_text(size = 10, family = "Times"))
load(args$input_file)

# selection_rule <- "zg"
# if (selection_rule == "zg") {
#   plot_sigmas <- c(0.1, 0.2)
# } else if (selection_rule == "elbow") {
#   plot_sigmas <- c(0.3, 2)
# }

theme_update(text = element_text(size = 10, family = "Times"))
color_map <- c("Choi et al. (2017)" = hue_pal()(3)[1], "Choi et al. (2017) [Bf]" = hue_pal()(3)[2], "Selective Inference" = hue_pal()(3)[3])
scale_map <- c("Choi et al. (2017)" = 16, "Choi et al. (2017) [Bf]" = 17, "Selective Inference" = 15)
linetype_map <- c("Selective Inference" = "solid", "Choi et al. (2017)" = "dashed")

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
    pve_ci_upper = mapply(get_signal_squared_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2,
    pve = pmin((signal / frob_norm)^2, 1),
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
    precision = 1 / sigma^2,
  ) %>%
  subset(
    method != "Choi et al. (2017) [Bf.]"
  ) %>%
  # subset(selection_r == rank) %>%
  # subset(sigma %in% plot_sigmas) %>%
  pivot_wider(
    id_cols = c(tested_k, rep, sigma, precision, selection_r),
    names_from = method,
    values_from = pve_ci_width
  ) %>%
  group_by(
    sigma, precision, tested_k
  ) %>%
  summarise(
    mean_pve_ci_diff = median(Selective / Unselective, na.rm = TRUE)
  )
head(plot_df)

# Scatter plot with different colors and shapes
g <- ggplot(
  plot_df,
  aes(x = tested_k, y = mean_pve_ci_diff, col = as.factor(sigma)),
) +
  geom_point() +
  labs(
    x = "Index k",
    y = "Width ratio",
    col = unname(TeX(r"( $\sigma$ )"))
  ) +
  theme_bw()
show(g)
ggsave(fname, width = 5.5, height = 2, unit = "in")
