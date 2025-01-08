library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis) 
library(stringr)
library(scales)
library(latex2exp)

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Selection_figures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

source("./scripts/functions/selection_methods.R")

vals <- c(114.32788 , 95.63982,  85.51193 , 66.00243  ,61.09360  ,54.05573  ,47.87799 , 31.80371 , 31.01390,  19.73305)
r1 <- select_r_elbow(vals)[1] + 1

# Case 2
k <- 5
bounds <- c(52.38685, 61.7774)

t2 <- bounds[2] + 2
temp_vals2 <- vals
temp_vals2[k] <- t2
r2 <- select_r_elbow(temp_vals2)[1] + 1

t3 <- bounds[1] - 10
temp_vals3<- vals
temp_vals3[k] <- t3
r3 <- select_r_elbow(temp_vals3)[1] + 1

## Plotting
vals <- sqrt(vals)
bounds <- sqrt(bounds)
t2 <- sqrt(t2)
t3 <- sqrt(t3)
temp_vals2 <- sqrt(temp_vals2)
temp_vals3 <- sqrt(temp_vals3)

# regular
g <- ggplot(
  data=data.frame(Index=seq(1,10)[-k], Eigenvalue=vals[-k]),
  aes(x=Index, y=Eigenvalue)) +
  theme(text=element_text(size=10)) +
  scale_x_continuous(breaks= pretty_breaks()) +
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=vals[k-1], ymax=Inf, alpha=.2, fill='yellow') +
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=vals[k+1], alpha=.2, fill='yellow') +
  geom_point() +
  geom_segment(aes(x = k, y = bounds[1], xend = k, yend = bounds[2]), color='grey', alpha=1) +
  geom_segment(aes(x = k-0.5, y = bounds[1], xend = k+0.5, yend = bounds[1]), color='grey', alpha=1) +
  geom_segment(aes(x = k-0.5, y = bounds[2], xend = k+0.5, yend = bounds[2]), color='grey', alpha=1) +
  geom_point(aes(x = k, y = vals[k]), color='chartreuse3', size=2) +
  geom_segment(aes(x = r1, y = vals[r1], xend = r1-1, yend = vals[r1-1]), color='chartreuse3', linetype='dashed') +
  geom_segment(aes(x = r1, y = vals[r1], xend = r1+1, yend = vals[r1+1]), color='chartreuse3', linetype='dashed') +
  theme_bw() +
  labs(y = "Squared singular values")
# print(g)
ggsave(paste0('./figures/Figure2-elbow-scree_plot-case2-true.png'), width = 2, height = 2, unit = "in")

# lower
g <- ggplot(
  data=data.frame(Index=seq(1,10)[-k], Eigenvalue=vals[-k]),
  aes(x=Index, y=Eigenvalue)) +
  theme(text=element_text(size=10)) +
  scale_x_continuous(breaks= pretty_breaks()) +
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=vals[k-1], ymax=Inf, alpha=.2, fill='yellow') +
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=vals[k+1], alpha=.2, fill='yellow') +
  geom_point() +
  geom_point(aes(x = k, y = t2), color='red', size=2) +
  geom_segment(aes(x = r2, y = temp_vals2[r2], xend = r2-1, yend = temp_vals2[r2-1]), color='red', linetype='dashed') +
  geom_segment(aes(x = r2, y = temp_vals2[r2], xend = r2+1, yend = temp_vals2[r2+1]), color='red', linetype='dashed') +
  geom_segment(aes(x = k, y = bounds[1], xend = k, yend = bounds[2]), color='grey', alpha=1) +
  geom_segment(aes(x = k-0.5, y = bounds[1], xend = k+0.5, yend = bounds[1]), color='grey', alpha=1) +
  geom_segment(aes(x = k-0.5, y = bounds[2], xend = k+0.5, yend = bounds[2]), color='grey', alpha=1) +
  theme_bw() +
  labs(y = "Squared singular values")
# print(g)
ggsave(paste0('./figures/Figure2-elbow-scree_plot-case2-lower.png'), width = 2, height = 2, unit = "in")

# upper
g <- ggplot(
  data=data.frame(Index=seq(1,10)[-k], Eigenvalue=vals[-k]),
  aes(x=Index, y=Eigenvalue)) +
  theme(text=element_text(size=10)) +
  scale_x_continuous(breaks= pretty_breaks()) +
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=vals[k-1], ymax=Inf, alpha=.2, fill='yellow') +
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=vals[k+1], alpha=.2, fill='yellow') +
  geom_point() +
  geom_point(aes(x = k, y = t3), color='red', size=2) +
  geom_segment(aes(x = r3, y = temp_vals3[r3], xend = r3-1, yend = temp_vals3[r3-1]), color='red', linetype='dashed') +
  geom_segment(aes(x = r3, y = temp_vals3[r3], xend = r3+1, yend = temp_vals3[r3+1]), color='red', linetype='dashed') +
  geom_segment(aes(x = k, y = bounds[1], xend = k, yend = bounds[2]), color='grey', alpha=1) +
  geom_segment(aes(x = k-0.5, y = bounds[1], xend = k+0.5, yend = bounds[1]), color='grey', alpha=1) +
  geom_segment(aes(x = k-0.5, y = bounds[2], xend = k+0.5, yend = bounds[2]), color='grey', alpha=1) +
  theme_bw() +
  labs(y = "Squared singular values")
# print(g)
ggsave(paste0('./figures/Figure2-elbow-scree_plot-case2-upper.png'), width = 2, height = 2, unit = "in")
