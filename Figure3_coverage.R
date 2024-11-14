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
parser$add_argument("input_file", nargs=1, help="File to be displayed")
parser$add_argument("--sigma", 
                    type = "double", 
                    default = 0.2)
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

# args <- list()
# args$input_file <- "data/sim_conf_ints_alpha=0.1_c=1_m=1_method=zg_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData"
# args$sigma <- 0.2
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

fname <- paste0("figures/Figure3_", sub("\\.RData$", ".png", basename(args$input_file)))
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

theme_update(text = element_text(size=10, family="Times"))
load(args$input_file)

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    # pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper, val_k, val_sum)[1, ],
    # pve_ci_upper = mapply(get_signal_squared_ci, ci_lower, ci_upper, val_k, val_sum)[2, ],
    # pve = signal^2 / ( val_sum - val_k + signal^2),
    pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
    pve_ci_upper = mapply(get_signal_squared_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2,
    pve = pmin((signal / frob_norm )^2, 1)
  )
head(results_df)

plot_df <- results_df %>%
  subset(sigma == args$sigma) %>%
  mutate(
    method = case_when(
      method == "Selective inference" ~ "Selective",
      method == "Choi" ~ "Unselective",
    ),
    precision = 1/sigma^2,
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
    legend.spacing.y = unit(0.05, 'in'),
    legend.margin=margin(t=0.05, unit='in'),
    legend.box.spacing=unit(0.1, "in"),
    legend.just = c("right", "top"),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),
    legend.key.size=unit(2,"lines")
  )
show(g)
ggsave(fname, width = 5.5, height = 2, unit = "in")
