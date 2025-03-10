library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(scales)
library(latex2exp)
library(argparse)

source("./scripts/functions/hypothesis_tests.R")
source("./scripts/functions/randomized_test.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")
source("./scripts/functions/confidence_intervals.R")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("input_file", nargs = 1, help = "File to be displayed")
parser$add_argument("--sigmas",
  type = "double",
  nargs = "+",
  default = c(0.1, 0.2)
)
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

# args <- list()
# args$input_file <- "data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=50_p=10_rank=5_reps=1000.RData"

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

power_fname <- paste0("figures/FigureApp_", sub("\\.RData$", "-power_smoothed.png", basename(args$input_file)))

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

alpha <- 0.1

theme_update(text = element_text(size = 10, family = "Times"))

load(args$input_file)

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    method = case_when(
      method == "Selective inference" ~ "Selective", # "Selective Inference",
      method == "Choi" ~ "Unselective", # "Choi et al. (2017)",
      method == "Choi (bf.)" ~ "Choi et al. (2017) [Bf.]",
    ),
    precision = 1 / sigma^2,
    reject_null = as.numeric((p_value <= alpha)),
    signal_normalized = signal / sigma,
    pve = signal^2 / frob_norm^2
  ) %>% subset(
    method != "Choi et al. (2017) [Bf.]"
  )
head(results_df)

g <- ggplot(
  results_df %>% subset(
    method == "Selective"
    # & tested_k == 1
    # & rank == selection_r
    & tested_k <= 6
  ),
  aes(x = signal_normalized, y = reject_null, color = as.factor(tested_k))
) +
  geom_point() +
  # geom_smooth(method = "loess", se = FALSE, span=200) +
  # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
  geom_smooth(
    formula = "y ~ x", method = "glm",
    method.args = list(family = "quasibinomial"), se = F
  ) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  labs(
    # x = "Signal",
    x = TeX(r"( Normalized signal $( u_k^T\{X\} \Theta\,\,v_k\{X\} / \sigma\,)$ )"),
    y = "Selective power",
    col = "Tested\nindex k"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0.04, NA)) +
  theme(
    # legend.direction = "horizontal",
    # ylim(0, NA),
    # legend.position = "bottom",
    # legend.position = c()
    # legend.box = "vertical",
    legend.spacing.y = unit(0.05, "in"),
    legend.margin = margin(t = 0.05, unit = "in"),
    legend.box.spacing = unit(0.1, "in"),
    # legend.box.just = "left"
    legend.just = c("right", "top"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(color = guide_legend(ncol = 1))
  # facet_wrap(~ tested_k, scales = "free_y")
# show(g)
ggsave(power_fname, width = 6, height = 2.5, unit = "in")
