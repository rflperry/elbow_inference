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
parser$add_argument("input_file", nargs=1, help="File to be displayed")
parser$add_argument("--sigmas", 
                    type = "double", 
                    nargs = "+",
                    default = c(0.1, 0.2))
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

# args <- list()
# args$input_file <- "data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=50_p=10_rank=5_reps=1000.RData"

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

detection_fname <- paste0("figures/Figure5_", sub("\\.RData$", "-detection_probability.png", basename(args$input_file)))
power_fname <- paste0("figures/Figure5_", sub("\\.RData$", "-power.png", basename(args$input_file)))

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

alpha <- 0.1

theme_update(text = element_text(size=10, family="Times"))

load(args$input_file)

results_df <- results_df %>%
  # drop_na() %>%
  mutate(
    method = case_when(
      method == "Selective inference" ~ "Selective", # "Selective Inference",
      method == "Choi" ~ "Unselective", # "Choi et al. (2017)",
      method == "Choi (bf.)" ~ "Choi et al. (2017) [Bf.]",
    ),
    precision = 1/sigma^2,
    reject_null = as.numeric((p_value <= alpha)),
  ) %>% subset(
    method != "Choi et al. (2017) [Bf.]"
  )
head(results_df)

plot_df <- results_df %>%
  # drop_na() %>%
  subset(selection_r == rank) %>%
  group_by(
    method, precision, tested_k
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
ggsave(detection_fname, width = 2, height = 2, unit = "in")

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
ggsave(power_fname, width = 4, height = 2, unit = "in")

