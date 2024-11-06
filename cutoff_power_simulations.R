# Test power of sequential choi vs our SI approach

library(ggplot2)
library(dplyr)
library(tidyr)
# library(RMTstat)
# library(scales)

source("./scripts/functions/hypothesis_tests.R")
source("./scripts/functions/selection_methods.R")
source("./scripts/functions/simulations.R")


################# Simulation under global null  ###########################

# Repeat experiment
test_sval_trial <- function (simulator, cutoff, alpha=0.05) {
  X <- simulator()
  duv <- svd(X)
  svals <- duv$d
  r <- which.min(c(svals >= cutoff, 0)) - 1
  if (r == 0) {
    return(c(NA,NA,NA,NA))
  }
  choi_pvalue <- choi_test_pvalue(svals, r)
  si_pvalue <- si_test_pvalue(svals, cutoff)
  
  cov_evals <- svd(cov(X))$d
  muirhead_pvalue <- muirhead_test_pvalue(cov_evals, r, n, p)
  pseudo_pvalue <- pseudorank_test_pvalue(svals^2, r, n, p)
  # return(c(si_pvalue, choi_pvalue, muirhead_pvalue))
  return(c(si_pvalue, choi_pvalue, muirhead_pvalue, pseudo_pvalue))
}

n <- 20
p <- 10
sigma <- 1 #/ sqrt(n)
c <- 6# sqrt(n) 

runs <- 1000

dfs <- list()
i <- 1
for (c in c(6, 3, 2)) {
  sim.pvalues <- replicate(runs, test_sval_trial(function() {dgf(n, p, sigma)}, cutoff=c))
  plot_df <- data.frame(
    pvalue=c(t(sim.pvalues)), # Need to take transpose!
    method=c(
      rep('Selective Inference', ncol(sim.pvalues)),
      rep('Choi', ncol(sim.pvalues)),
      rep('Muirhead', ncol(sim.pvalues)),
      rep('Pseudo-rank', ncol(sim.pvalues))
    )
  )
  plot_df$cutoff <- c
  dfs[[i]] <- plot_df
  i <- i + 1
}

plot_df <- rbind(dfs[[1]], dfs[[2]], dfs[[3]])

ggplot(
  data=plot_df,
  aes(sample=pvalue, col=method))+
  geom_qq(
    distribution="qunif"
  )+
  geom_abline(col="black", lty = 2) +
  xlab('Uniform quantiles') +
  ylab('Sample quantiles') +
  facet_wrap(~ cutoff) +
  ggtitle(paste0("p-values under the global null, across cutoffs"))
ggsave(paste0('./figures/cutoff_rule/global_null_validity.png'), width = 8, height = 4, unit = "in")

################# Simulation under global null, with sequential  ###########################


# Repeat experiment
test_sval_trial <- function (simulator, cutoff, n, p, sigma, alpha=0.05) {
  X <- simulator()
  duv <- svd(X)
  svals <- duv$d
  selection_results <- select_r_cutoff(svals, cutoff)
  r <- selection_results[1]
  lb <- selection_results[2]
  ub <- selection_results[3]
  choi_pvalue <- (choi_test_pvalue(svals, r, n, p, sigma) <= alpha)
  si_pvalue <- (si_test_pvalue(svals, r, n, p, sigma, lb=lb, ub=ub) <= alpha)
  seq_choi_pvalue <- (choi_test_pvalue(svals, r, n, p, sigma) <= (alpha / (r+1)))
  return(c(si_pvalue, choi_pvalue, seq_choi_pvalue))
}

n <- 20
p <- 10
sigma <- 1 #/ sqrt(n)
c <- 6 # sqrt(n)

runs <- 1000

sim.rejections <- replicate(
  runs, test_sval_trial(
    function() {simulate_global_null(n=n, p=p, sigma=1)}, c=c, n, p, sigma
    )
  )

plot_df <- data.frame(
  rejects=c(t(sim.rejections)), # Need to take transpose!
  method=c(
    rep('Selective Inference', runs),
    rep('Choi', runs),
    rep('Corrected Choi', runs)
  )
)

plot_df %>% group_by(method) %>% summarise(t1_error = mean(rejects, na.rm=TRUE))

################# Simulation under non-null, power comparison  ###########################
n <- 50
p <- 10
sigma <- 1 #/ sqrt(n)
rank <- 5
m <- 1
# m_list <- c(0, 0.5, 1, 2, 5)# c(0, 0.2, 0.5, 0.8, 1, 1.3, 1.7, 2, 2.3, 2.6, 3)
sigma_list <- c(0.1, 0.2, 0.3, 0.5, 1, 2)
cutoff_quantiles <- c(0.1, 0.5, 0.9)

# Choose cutoffs as midpoints
cutoff_runs <- 100
cutoffs <- array(rep(NA, length(cutoff_quantiles)*length(sigma_list)), c(length(sigma_list), length(cutoff_quantiles)))
for (index in 1:length(sigma_list)) {
  sigma <- sigma_list[index]
  midpoints <- c()
  for (rep in 1:100) {
    midpoints <- c(midpoints, mean(svd(dgf(n, p, sigma, rank, m)$obsv_mat)$d[rank:(rank+1)]))
  }
  for (quant_i in 1:length(cutoff_quantiles) ) {
    cutoffs[index, quant_i] <- quantile(midpoints, cutoff_quantiles[quant_i])
  }
}

runs <- 1000
alpha <- 0.05

results <- c()

for (run in 1:runs) {
  for (i in 1:length(sigma_list)) {
    sigma <- sigma_list[i]
    sim <- dgf(n, p, sigma, rank, m)
    duv <- svd(sim$obsv_mat)
    
    svals <- duv$d
    for (k in 1:length(cutoff_quantiles)) {
      cutoff <- cutoffs[i, k]
      r <- which.min(c(svals >= cutoff, 0)) - 1
      
      # "leftover" signal
      projd_mat <- (duv$u[,r:p] %*% t(duv$u[,r:p])) %*% sim$mean_mat %*% (duv$v[,r:p] %*% t(duv$v[,r:p]))
      null_strength <- norm(projd_mat, type='2') # maximum singular value

      base_results <- c(n, p, sigma, m, cutoff, cutoff_quantiles[k], r, rank, null_strength)

      try ({
        si_pvalue <- si_test_pvalue(svals, cutoff, n, p, sigma)
        results <- c(results, 'Selective inference', base_results, si_pvalue)
        
        choi_pvalue <- choi_test_pvalue(svals, k=r, n, p, sigma)
        results <- c(results, 'Choi', base_results, choi_pvalue)
      })
    }
  }
}

header <- c(
  "method", "n_samples", "n_observations", "sigma", "signal_strength",
  "cutoff", "cutoff_quantile", 'selection_r', 'rank', 'null_strength', "pvalue")
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,2:length(header)] <- sapply(results_df[,2:length(header)], as.numeric)
head(results_df)

results_df %>%
  group_by(sigma) %>%
  summarise(
    count = n()
  )

save(results_df, file="data/cutoff_rule/cutoff-power_vs_sigma_comparison.RData")

summary <- results_df %>%
  mutate(
    detected_r = as.integer(selection_r == rank),
    selection_r = as.integer(selection_r),
    inverse_sigma = 1 / sigma,
  ) %>%
  mutate(rejected = as.integer(pvalue <= ifelse( method == 'Choi', alpha / (selection_r + 1), alpha ))) # %>%


# plot power vs signal in the null
g <- ggplot(
  summary,# %>% subset(detected_r == TRUE),# %>% group_by(signal_strength, method) %>% summarise(rejected = mean(rejected)),
  aes(x=null_strength, y=rejected, col=method)) +
  # aes(x=signal_strength, y=pvalue+1e-10, col=method)) +
  # geom_smooth(method='loess') +
  stat_summary_bin(fun = mean, bins=10, geom='line') +
  # geom_line() +
  # geom_point(size=1) +
  # xlab('Relative distance to cutoff') +
  xlab(expression('Signal in null: '~'||'~P[U['[r:p]']]~Theta~P[V['[r:p]']]~'||'[2])) +
  ylab('Power') +
  labs(color='Method') +
  # geom_hline(yintercept=0.05, linetype="dashed",color = "gray", size=0.5)
  # scale_x_discrete(breaks=seq(0, 1, 0.5)) +
  scale_x_continuous(trans='log10') +
  facet_wrap( ~ cutoff_quantile)
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-power_vs_null_strength.png'), width = 6, height = 3, unit = "in")

# plot pvalue vs signal strength
g <- ggplot(
  summary %>%
    subset(detected_r == 1) %>%
    group_by(signal_strength, method, cutoff_quantile) %>%
    summarise(pvalue = mean(pvalue)) %>%
    mutate(method = replace(method, method == "Corrected Choi", "Choi")),
  aes(x=signal_strength, y=pvalue)) +# +1e-10, col=method)) +
  # aes(x=signal_strength, y=pvalue+1e-10, col=method)) +
  # geom_smooth(method='loess') +
  # stat_summary_bin(fun = mean, bins=10, geom='line') +
  geom_line(aes(color=method))+
  geom_point(size=1) +
  xlab('Signal_strenth') +
  # xlab(expression('Signal in null: '~'||'~P[U['[r:p]']]~Theta~P[V['[r:p]']]~'||'[2])) +
  ylab('p value') +
  labs(color='Method') +
  # geom_hline(yintercept=0.05, linetype="dashed",color = "gray", size=0.5)
  # scale_x_discrete(breaks=seq(0, 1, 0.5)) +
  # scale_y_continuous(trans='pseudo_log') +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 1e-6)) +
  facet_wrap( ~ cutoff_quantile)
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-pvalue_vs_signal_strength.png'), width = 6, height = 3, unit = "in")

# plot pvalue ratio vs signal strength
g <- ggplot(
  summary %>%
    # group_by(signal_strength, cutoff_quantile, method) %>%
    # summarise(pvalue_ratio = mean(pvalue_ratio)),
    # summarise(denominator = mean(denominator), numerator = mean(numerator)) %>%
    mutate(method = replace(method, method == "Corrected Choi", "Choi")),
  # aes(x=signal_strength, y=denominator,col=method)) +# +1e-10, col=method)) +
  aes(x=null_strength, y=pvalue, col=method)) +
  # geom_smooth(method='loess') +
  stat_summary_bin(fun = mean, bins=10, geom='line') +
  # geom_line(aes(color=method))+
  # geom_point(size=1) +
  # xlab('Signal_strenth') +
  xlab(expression('Signal in null: '~'||'~P[U['[r:p]']]~Theta~P[V['[r:p]']]~'||'[2])) +
  ylab('p value ratio') +
  labs(color='Method') +
  # geom_hline(yintercept=0.05, linetype="dashed",color = "gray", size=0.5)
  # scale_x_discrete(breaks=seq(0, 1, 0.5)) +
  scale_y_continuous(trans='log10') +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 1e-6)) +
  facet_wrap( ~ cutoff_quantile)
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-pvalue_vs_null_strength.png'), width = 6, height = 3, unit = "in")

################3 orignal splots

summary <- results_df %>%
  mutate(
    detected_r = as.integer(selection_r == rank),
    selection_r = as.integer(selection_r),
    inverse_sigma = 1/sigma,
  ) %>%
  mutate(rejected = as.integer(pvalue <= ifelse( method == 'Corrected Choi', alpha / (selection_r + 1), alpha ))) %>%
  group_by(method, inverse_sigma, cutoff_quantile) %>%
  summarise(
    detection_probability = mean(detected_r),
    conditional_power = sum(detected_r * rejected) / sum(detected_r),
    conditional_pvalue = sum(detected_r * pvalue) / sum(detected_r),
    mean_pvalue = mean(pvalue))


# plot detection probability
g <- ggplot(
  summary %>% subset(method == 'Corrected Choi'),
  aes(x=signal_strength, y=detection_probability, group=1)) +
  geom_line(linetype = 'dashed') +
  geom_point(size=1) +
  xlab('Signal strength') +
  ylab('Detection probability') +
  scale_x_discrete(breaks=seq(0, 3, 0.5)) +
  facet_wrap(~ cutoff_quantile)
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-detection_probability_vs_signal_strength.png'), width = 6, height = 3, unit = "in")

# plot conditional power
g <- ggplot(summary, aes(x=inverse_sigma, y=conditional_power, group=method)) +
  geom_line(aes(color=method)) + #, linetype = 'dashed')+
  geom_point(size=1) +
  # xlab('Signal strength') +
  xlab("Inverse sigma") + 
  ylab('Conditional power') +
  labs(color='Method') +
  geom_hline(yintercept=0.05, linetype="dashed",color = "gray", size=0.5) +
  # scale_x_discrete(breaks=seq(0, 3, 0.5)) +
  facet_wrap( ~ cutoff_quantile) +
  scale_x_continuous(trans='log2')
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-detection_probability_vs_signal_strength.png'), width = 6, height = 3, unit = "in")

# plot pvalue
g <- ggplot(summary, aes(x=signal_strength, y=mean_pvalue, group=method)) +
  geom_line(aes(color=method), linetype = 'dashed')+
  geom_point(size=1) +
  xlab('Signal strength') +
  ylab('pvalue') +
  labs(color='Method') +
  scale_x_discrete(breaks=seq(0, 3, 0.5)) +
  facet_wrap( ~ cutoff_quantile) +
  scale_y_continuous(trans='log2')
plot(g)
ggsave(paste0('figures/power_comparison/pvalue_vs_signal_strength.png'), width = 6, height = 3, unit = "in")


########################### p value verification ############################

# Integration integrand, proportional to likelihood
integrand <- function(s_obsv, s_other, tr_mean, sigma=1) {
  return(
    mapply(function(s) {
      exp(-s**2 / 2 / sigma**2 + s * tr_mean / sigma**2) * s ** (n-p) * prod(abs(s_other**2 - s**2))
    },
    s_obsv
    )
  )
}

# Choi pvalue that the kth singular value is greater than zero
choi_test_pvalue <- function(svals, k, tr_mean, sigma=1) {
  bounds <- c(Inf, svals, 0)
  s_obs <- svals[k]
  numerator <- integrate(integrand, s_obs, bounds[k], svals[-k], tr_mean, sigma)$value
  denominator <- integrate(integrand, bounds[k+2], bounds[k], svals[-k], tr_mean, sigma)$value
  pvalue <- numerator / denominator
  if(pvalue < 0 | pvalue > 1){
    stop(paste0("Choi pvalue invalid ", pvalue))
  }
  return(c(pvalue, numerator, denominator))
}

# SI pvalue that the rth (smallest singular value greater than the cutoff)
# is greater than zero
si_test_pvalue <- function(svals, cutoff, tr_mean, sigma=1) {
  r <- which.min(c(svals >= cutoff, 0)) - 1
  if (r == 0) {
    return(NA)
  }
  bounds <- c(Inf, svals, 0)
  s_obs <- svals[r]
  numerator <- integrate(integrand, s_obs, bounds[r], svals[-r], tr_mean, sigma)$value
  denominator <- integrate(integrand, cutoff, bounds[r], svals[-r], tr_mean, sigma)$value
  pvalue <- numerator / denominator
  if(pvalue < 0 | pvalue > 1){
    stop(paste0("SI pvalue invalid ", pvalue))
  }
  return(c(pvalue, numerator, denominator))
}

# Signal
dgf <- function (n, p, sigma, rank, m){
  UV <- array(rnorm(n*p, sd=sigma), dim=c(n, p))
  duv <- svd(UV)
  U <- duv$u
  V <- duv$v
  svals <- c(rev(seq(1:rank)), rep(0, p-rank)) * m * sigma * (n*p)^(1/4)
  mean_mat <- U %*% diag(svals) %*% t(V)
  results <- {}
  results$duv <- duv
  results$mean_mat <- mean_mat
  results$obsv_mat <- mean_mat + array(rnorm(n*p, sd=sigma), dim=c(n, p))
  return(results)
}

n <- 50
p <- 10
sigma <- 1 #/ sqrt(n)
rank <- 5
m_list <- c(0, 0.5, 1, 2, 5)# c(0, 0.2, 0.5, 0.8, 1, 1.3, 1.7, 2, 2.3, 2.6, 3)
cutoff_quantiles <- c(0.1, 0.5, 0.9)

# Choose cutoffs as midpoints
cutoff_runs <- 100
cutoffs <- array(rep(NA, length(cutoff_quantiles)*length(m_list)), c(length(m_list), length(cutoff_quantiles)))
for (m_index in 1:length(m_list)) {
  m <- m_list[m_index]
  midpoints <- c()
  for (rep in 1:100) {
    midpoints <- c(midpoints, mean(svd(dgf(n, p, sigma, rank, m)$obsv_mat)$d[rank:(rank+1)]))
  }
  for (quant_i in 1:length(cutoff_quantiles) ) {
    cutoffs[m_index, quant_i] <- quantile(midpoints, cutoff_quantiles[quant_i])
  }
}

runs <- 500
alpha <- 0.05

results <- c()

for (run in 1:runs) {
  for (i in 1:length(m_list)) {
    m <- m_list[i]
    sim <- dgf(n, p, sigma, rank, m)
    duv <- svd(sim$obsv_mat)
    
    # duv$d <- sapply(
    #   seq(1,p),
    #   function(k) {
    #     (t(duv$u[,k]) %*% sim$obsv_mat %*% duv$v[,k])
    #   }
    # )
    # 
    # dec_order <- rev(order(duv$d))
    # duv$d <- duv$d[dec_order]
    # duv$u <- duv$u[,dec_order]
    # duv$v <- duv$v[,dec_order]
    # 
    # if (!sum(diff(duv$d) < 0) == (p-1)) {
    #   stop(paste0('Invalid ordering: ', duv$d))
    # }
    
    svals <- duv$d
    for (k in 1:length(cutoff_quantiles)) {
      cutoff <- cutoffs[i, k]
      r <- which.min(c(svals >= cutoff, 0)) - 1
      
      # "leftover" signal
      projd_mat <- (duv$u[,r:p] %*% t(duv$u[,r:p])) %*% sim$mean_mat %*% (duv$v[,r:p] %*% t(duv$v[,r:p]))
      null_strength <- norm(projd_mat, type='2') # maximum singular value
      
      tr_mean <- (t(duv$u[,r]) %*% sim$mean_mat %*% duv$v[,r])[1]
      
      base_results <- c(svals[rank+1], svals[rank], svals[rank-1], m, cutoff, cutoff_quantiles[k], r, rank, null_strength)
      
      pvals <- test_pvalues(svals, r)
      
      si_pvalue <- si_test_pvalue(svals, cutoff, tr_mean)
      results <- c(results, 'Selective inference', 'Full', base_results, si_pvalue)
      
      si_pvalue <- si_test_pvalue(svals, cutoff, 0)
      results <- c(results, 'Selective inference', 'Null', base_results, si_pvalue)
      
      si_pvalue <- si_test_pvalue(svals, cutoff, 1.1*tr_mean)
      results <- c(results, 'Selective inference', 'Extra', base_results, si_pvalue)
      
      choi_pvalue <- choi_test_pvalue(svals, k=r, tr_mean)
      results <- c(results, 'Choi', 'Full', base_results, choi_pvalue)
      
      choi_pvalue <- choi_test_pvalue(svals, k=r, 0)
      results <- c(results, 'Choi', 'Null', base_results, choi_pvalue)
      
      choi_pvalue <- choi_test_pvalue(svals, k=r, 1.1*tr_mean)
      results <- c(results, 'Choi', 'Extra', base_results, choi_pvalue)
    }
  }
}

header <- c("method", "likelihood", "sval_below", "target_sval", "sval_above", "signal_strength", "cutoff", "cutoff_quantile", 'selection_r', 'true_rank', 'null_strength', "pvalue", "numerator", "denominator")
results_df <- data.frame(
  data=t(array(results, dim = c(length(header), length(results)/length(header)))))
colnames(results_df) <- header
results_df[,3:length(header)] <- sapply(results_df[,3:length(header)], as.numeric)
head(results_df)

save(results_df, file="data/cutoff_rule/cutoff-pvalue_vs_signal_strength_comparison.RData")

summary <- results_df %>%
  mutate(
    detected_r = as.integer(selection_r == true_rank),
    selection_r = as.integer(selection_r),
    distance_to_cutoff = (target_sval - cutoff),# / (target_sval - sval_below),
    sval_separation = (target_sval - sval_below),# / (sval_above - sval_below),
  ) %>%
  mutate(rejected = as.integer(pvalue <= ifelse( method == 'Choi', alpha / (selection_r + 1), alpha ))) # %>%

# plot pvalue vs signal strength: 2-way facet across likelihoods null/alt
g <- ggplot(
  # summary %>%
  summary %>%
    subset(selection_r == (true_rank)) %>%
    group_by(signal_strength, method, cutoff_quantile, likelihood) %>%
    summarise(pvalue = mean(pvalue)) %>%
    mutate(method = replace(method, method == "Corrected Choi", "Choi")),
  aes(x=signal_strength, y=pvalue, col=method)) +# +1e-10, col=method)) +
  # aes(x=signal_strength, y=pvalue+1e-10, col=method)) +
  # geom_smooth(method='loess') +
  # stat_summary_bin(fun = mean, bins=30, geom='line') +
  geom_line(aes(color=method))+
  # geom_point(size=1) +
  xlab('Signal strength') +
  # xlab(expression('Signal in null: '~'||'~P[U['[r:p]']]~Theta~P[V['[r:p]']]~'||'[2])) +
  ylab('p value') +
  labs(color='Method') +
  # geom_hline(yintercept=0.05, linetype="dashed",color = "gray", size=0.5)
  # scale_x_discrete(breaks=seq(0, 1, 0.5)) +
  scale_y_continuous(trans='log', labels=function(x) format(x, nsmall = 1)) +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 1e-6)) +
  facet_wrap(likelihood ~ cutoff_quantile, scales="free_y") # 
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-pvalue_vs_signal_strength.png'), width = 12, height = 6, unit = "in")

# plot pvalue vs signal strength: 2-way facet scatterplots
temp_df <- summary %>%
  mutate(method = replace(method, method == "Corrected Choi", "Choi")) %>%
  filter(likelihood == 'Null', cutoff_quantile == 0.5) %>%
  pivot_wider(
    names_from = method, 
    values_from = pvalue,
    names_glue = "{method} pvalue"
  ) %>%
  rename('SI_pvalue' = 'Selective inference pvalue', 'Choi_pvalue' = 'Choi pvalue')

# temp_df <- temp_df[1:44964,]
# temp_df_si <- temp_df[rep(c(T, T, T, F, F, F), nrow(temp_df) / 6),]
# temp_df <- temp_df[rep(c(F, F, F, T, T, T), nrow(temp_df) / 6),]
# temp_df$SI_pvalue <- temp_df_si$SI_pvalue

temp_df_si <- temp_df[rep(c(T, F), nrow(temp_df) / 2),]
temp_df <- temp_df[rep(c(F, T), nrow(temp_df) / 2),]
temp_df$SI_pvalue <- temp_df_si$SI_pvalue


g <- ggplot(
  # summary %>%
  # summary %>%
  #   mutate(method = replace(method, method == "Corrected Choi", "Choi")) %>%
  #   filter(likelihood == 'Null', cutoff_quantile == 0.5) %>%
  #   pivot_wider(
  #     names_from = method, 
  #     values_from = pvalue,
  #     names_glue = "{method} pvalue"
  #   ),
  # temp_df,
  temp_df %>% mutate(
    SI_pvalue = ifelse(SI_pvalue <= 1e-10, 1e-10, SI_pvalue),
    Choi_pvalue = ifelse(Choi_pvalue <= 1e-10, 1e-10, Choi_pvalue),
  ),#%>% subset(signal_strength <= 3),
    # rename('SI pvalue' = 'Selective inference pvalue', 'Choi_'),
    # subset(selection_r == (true_rank)) %>%
    # group_by(signal_strength, method, cutoff_quantile, likelihood) %>%
    # summarise(pvalue = mean(pvalue)) %>%
  aes(x=Choi_pvalue, y=SI_pvalue, col=null_strength)) +# +1e-10, col=method)) +
  # aes(x=signal_strength, y=pvalue+1e-10, col=method)) +
  # geom_smooth(method='loess') +
  # stat_summary_bin(fun = mean, bins=30, geom='line') +
  # geom_line(aes(color=method))+
  geom_point() +
  geom_abline(slope=1, intercept=0, col='red', linetype='dashed') +
  # xlab('Signal strength') +
  # xlab(expression('Signal in null: '~'||'~P[U['[r:p]']]~Theta~P[V['[r:p]']]~'||'[2])) +
  # ylab('p value') +
  # labs(color='Method') +
  # geom_hline(yintercept=0.05, linetype="dashed",color = "gray", size=0.5)
  # scale_x_discrete(breaks=seq(0, 1, 0.5)) +
  # scale_y_continuous(trans='log', labels=function(x) format(x, nsmall = 1)) +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 1e-6)) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10')
  # facet_wrap(likelihood ~ cutoff_quantile)#, scales="free_y") # 
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-pvalue_vs_null_strength_scatter_clipped.png'), width = 5, height = 3, unit = "in")

# plot pvalue ratio vs signal strength
g <- ggplot(
  summary %>%
    group_by(signal_strength, cutoff_quantile, method) %>%
    # summarise(pvalue_ratio = mean(pvalue_ratio)),
    summarise(denominator = mean(denominator), numerator = mean(numerator)) %>%
    mutate(method = replace(method, method == "Corrected Choi", "Choi")),
  aes(x=signal_strength, y=denominator,col=method)) +# +1e-10, col=method)) +
  # aes(x=signal_strength, y=pvalue+1e-10, col=method)) +
  # geom_smooth(method='loess') +
  # stat_summary_bin(fun = mean, bins=10, geom='line') +
  geom_line(aes(color=method))+
  geom_point(size=1) +
  xlab('Signal_strenth') +
  # xlab(expression('Signal in null: '~'||'~P[U['[r:p]']]~Theta~P[V['[r:p]']]~'||'[2])) +
  ylab('p value ratio') +
  labs(color='Method') +
  # geom_hline(yintercept=0.05, linetype="dashed",color = "gray", size=0.5)
  # scale_x_discrete(breaks=seq(0, 1, 0.5)) +
  scale_y_continuous(trans='log10') +
  # scale_y_continuous(trans=scales::pseudo_log_trans(base = 10, sigma = 1e-6)) +
  facet_wrap( ~ cutoff_quantile)
plot(g)
ggsave(paste0('figures/cutoff_rule/cutoff-pvalue_vs_signal_strength.png'), width = 6, height = 3, unit = "in")
