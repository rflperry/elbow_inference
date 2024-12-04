library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(stringr)
library(scales)
library(latex2exp)
library(argparse)

source("./scripts/functions/confidence_intervals.R")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("input_files", nargs = "+", help = "Files to be displayed")
parser$add_argument("--sigma",
    type = "double",
    default = 0.2
)
print(commandArgs(trailingOnly = TRUE))
args <- parser$parse_args()

args <- list()
args$input_files <- c(
    "data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData",
    "data/sim_conf_ints_alpha=0.3_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData",
    "data/sim_conf_ints_alpha=0.5_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData",
    "data/sim_conf_ints_alpha=0.7_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData",
    "data/sim_conf_ints_alpha=0.9_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData"
)
args$sigma <- 0.1
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

fname <- paste0("figures/Figure_coverage_vs_alpha_", sub("\\.RData$", ".png", basename(args$input_files[1])))

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

theme_update(text = element_text(size = 10, family = "Times"))
# Load all dataframes in the list "input_files" and concatenate them by rows
data <- lapply(args$input_files, function(file) {
    load(file)
    results_df$alpha <- as.numeric(str_extract(file, "alpha=([0-9.]+)") %>% str_extract("[0-9.]+"))
    get(ls(pattern = "results_df"))
})
results_df <- bind_rows(data)

results_df <- results_df %>%
    # drop_na() %>%
    mutate(
        # pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper, val_k, val_sum)[1, ],
        # pve_ci_upper = mapply(get_signal_squared_ci, ci_lower, ci_upper, val_k, val_sum)[2, ],
        # pve = signal^2 / ( val_sum - val_k + signal^2),
        pve_ci_lower = mapply(get_signal_squared_ci, ci_lower, ci_upper)[1, ] / frob_ci_upper^2,
        pve_ci_upper = mapply(get_signal_squared_ci, ci_lower, ci_upper)[2, ] / frob_ci_lower^2,
        pve = pmin((signal / frob_norm)^2, 1)
    )
head(results_df)

plot_df <- results_df %>%
    subset(sigma == 0.1) %>%
    mutate(
        method = case_when(
            method == "Selective inference" ~ "Selective",
            method == "Choi" ~ "Unselective",
        ),
        precision = 1 / sigma^2,
        pve_covered = as.numeric((pve_ci_lower <= pve) & (pve <= pve_ci_upper)),
        frob_covered = as.numeric((frob_ci_lower <= frob_norm) & (frob_norm <= frob_ci_upper)),
        signal_covered = as.numeric((signal <= pve) & (pve <= pve_ci_upper)),
        pve_ci_width = pve_ci_upper - pve_ci_lower,
    ) %>%
    subset(method == 'Selective') %>%
    group_by(
        method, tested_k, alpha
    ) %>%
    summarise(
        mean_pve_ci_lower = median(pve_ci_lower, na.rm = TRUE),
        mean_pve = median(pve, na.rm = TRUE),
        mean_pve_ci_upper = median(pve_ci_upper, na.rm = TRUE),
        coverage = mean(pve_covered, na.rm = TRUE),
        frob_coverage = mean(frob_covered, na.rm = TRUE),
        mean_pve_width = median(pve_ci_width, na.rm = TRUE),
        se = sqrt(coverage * (1 - coverage) / n())
    )
head(plot_df)

g <- ggplot(
    plot_df,
    aes(
        x = 1 - alpha,
        y = coverage,
        fill = as.factor(tested_k),
        ymin = pmax(coverage - 2 * se, 0),
        ymax = pmin(coverage + 2 * se, 1)
    )
) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(width = 0.18, position = position_dodge()) +
    labs(
        x = TeX(r"( Coverage level $(1 - \,\alpha\,)$ )"),
        y = "Empirical coverage",
        fill = 'Tested k'
        ) +
    scale_x_continuous(
        breaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
        labels = c(0.1, 0.3, 0.5, 0.7, 0.9)
        ) +
    scale_y_continuous(limits = c(0, 1)) +
    coord_cartesian(ylim = c(0.04, NA)) +
    scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw() +
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
    )
show(g)
ggsave(fname, width = 6, height = 2, unit = "in")

