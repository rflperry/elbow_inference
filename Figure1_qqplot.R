library(dplyr)
library(tidyr)

source("./scripts/functions/hypothesis_tests.R")
source("./scripts/functions/randomized_test.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Global null: run simulation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

n = 50
p = 10
sigma = 1
rank = 5
m = 0
eigen = TRUE
reps <- 10000

results <- c()

for (rep in 1:reps) {
  # try({
  # Simulate data
  sim <- simulate_matrix(n, p, sigma, rank, m, eigen=eigen)
  
  # Extract quantities
  duv <- svd(sim$obsv_mat)# - colMeans(sim$obsv_mat))
  # rand_duv <- get_randomized_svals(sim$obsv_mat, tau=tau, sigma=sigma)
  
  if (eigen) {
    vals <- duv$d^2
    # rand_vals <- rand_duv$d^2
  } else {
    vals <- duv$d
    # rand_vals <- rand_duv$d
  }
  
  # Perform selection
  # selection_event <- select_r_elbow(vals)
  # r <- selection_event[1]
  r <- select_r_elbow(vals)[1]
  
  for (k in seq(1, r)) {
    # Base quantities to save
    base_results <- c(n, p, sigma, m, rank, r, k)
    
    # Choi p-value
    choi_pvalue <- choi_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma, eigen=eigen)
    results <- c(results, c('Choi', c(rep, choi_pvalue, base_results)))
    
    # Selective inference p-value (R >= k)
    bounds <- compute_elbow_bounds(vals, k)
    
    si_pvalue <- si_test_pvalue(vals, k=k, n=n, p=p, sigma=sigma, bounds=bounds, eigen=eigen)
    results <- c(results, c('Selective inference', c(rep, si_pvalue, base_results)))
  }
  # })
}

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Process data and save
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

header <- c(
  "method", "rep", "pvalue", "n", "p", "sigma", "signal_strength",
  'rank', 'selection_r', 'tested_k')
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
tail(results_df)
print(nrow(results_df[complete.cases(results_df),]) / nrow(results_df))

save(results_df, file="data/elbow-power-global_null.RData")

#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 	Plot the figure
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis) 
library(stringr)
library(scales)
library(latex2exp)

load("data/elbow-power-global_null.RData")
results_df <- results_df %>% mutate(
  method = case_when(
    method == "Selective inference" ~ "Selective Inference",
    method == "Choi" ~ "Choi et al. (2017)",
  )
)

theme_update(text = element_text(size=10, family="Times"))
color_map <- c("Choi et al. (2017)" = hue_pal()(3)[1], "Choi et al. (2017) [Bf]" = hue_pal()(3)[2], "Selective Inference" = hue_pal()(3)[3])
scale_map <- c("Choi et al. (2017)" = 16, "Choi et al. (2017) [Bf]" = 17, "Selective Inference" = 15)
linetype_map <- c( "Selective Inference" = "solid", "Choi et al. (2017)" = "dashed")

g <- ggplot(
  data=results_df,
    # subset(tested_k <= 5) %>%
  # aes(x=pvalue, linetype=method, col=tested_k+2, group=interaction(method, tested_k))) +
  aes(x=1 - pvalue, linetype=method, col=as.factor(tested_k))) +
  stat_ecdf()+
  geom_abline(col="black", lty = 3) +
  ylab('Empirical quantiles') +
  xlab('Uniform distribution quantiles') +
  labs(linetype='Method', col='Tested') +
  # scale_color_viridis(option = "D")+
  # scale_color_viridis(discrete = TRUE, option =  "D") +
  scale_color_brewer(palette="Blues", type="seq") + #, limits=as.factor(seq(0, 8)), breaks=seq(1, 8)) +
  scale_linetype_manual(values = linetype_map) +
  theme_bw() +# legend.position="bottom", legend.box="vertical", 
  theme(legend.box="horizontal", text=element_text(size=10)) +
  guides(color=guide_legend(ncol=2))
print(g)
ggsave(paste0('./figures_FINAL/Figure1_qqplot_elbow.png'), width = 5, height = 2, unit = "in")

