# #
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 	ZG Rule
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #

# QQ-plot
# Rscript Figure1_qq_plot.R

# QQ-plot, estimated variance
# Rscript FigureApp_qq_plot.R

# Power: regular
# Rscript Figure4_power.R \
#     data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=50_p=10_rank=5_reps=1000.RData

# Power: regular, smoothed
# Rscript FigureApp_power_smoothed.R \
#     data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=50_p=10_rank=5_reps=1000.RData

# # Power: low rank, smoothed
# Rscript FigureApp_power_smoothed.R \
#     data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=50_p=10_rank=2_reps=1000.RData

# # Power: high rank, smoothed
# Rscript FigureApp_power_smoothed.R \
#     data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=50_p=10_rank=10_reps=1000.RData

# # Power: high dim, smoothed
# Rscript FigureApp_power_smoothed.R \
#     data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=100_p=50_rank=5_reps=1000.RData

# # Power: estimated variance
# Rscript Figure4_power.R \
#     data/sim_hypo_tests_alpha=0.1_m=1_method=zg_n=50_p=10_rank=5_reps=1000_var_est=TRUE.RData

# CI Widths: stratified
# Rscript FigureApp_confidence_interval_widths_stratified.R \
#     data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=1000_signal_alpha_frac=0.75.RData \
#     --sigmas 0.1 0.2

# # # CI Widths: regular
# Rscript Figure6_confidence_interval_widths.R \
#     data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=1000_signal_alpha_frac=0.75.RData \
#     --sigmas 0.1 0.2

# # # CI Widths: low rank
# Rscript Figure6_confidence_interval_widths.R \
#     data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=2_reps=1000_signal_alpha_frac=0.75.RData \
#     --sigmas 0.1 0.2

# # CI Widths: higher rank
# Rscript Figure6_confidence_interval_widths.R \
#     data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=10_reps=1000_signal_alpha_frac=0.75.RData \
#     --sigmas 0.1 0.2

# # CI Widths: higher dim
# Rscript Figure6_confidence_interval_widths.R \
#     data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=100_p=50_rank=5_reps=1000_signal_alpha_frac=0.75.RData \
#     --sigmas 0.1 0.2

# # CI Widths: regular
Rscript Figure6_confidence_interval_widths.R \
    data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=1000_signal_alpha_frac=0.75_var_est=TRUE.RData \
    --sigmas 0.1 0.2

# # CI coverage: across alphas with estimation variance
# Rscript Figure5_coverage_by_alpha.R \
#     "data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75_var_est=TRUE.RData" \
#     "data/sim_conf_ints_alpha=0.3_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75_var_est=TRUE.RData" \
#     "data/sim_conf_ints_alpha=0.5_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75_var_est=TRUE.RData" \
#     "data/sim_conf_ints_alpha=0.7_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75_var_est=TRUE.RData" \
#     "data/sim_conf_ints_alpha=0.9_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75_var_est=TRUE.RData"

# CI coverage: across alphas
# Rscript Figure5_coverage_by_alpha.R \
#     "data/sim_conf_ints_alpha=0.1_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData" \
#     "data/sim_conf_ints_alpha=0.3_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData" \
#     "data/sim_conf_ints_alpha=0.5_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData" \
#     "data/sim_conf_ints_alpha=0.7_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData" \
#     "data/sim_conf_ints_alpha=0.9_c=1_choi=FALSE_m=1_method=zg_mle=TRUE_n=50_p=10_rank=5_reps=10000_signal_alpha_frac=0.75.RData"

# #
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 	Elbow Rule
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #

# # CI Widths: regular
# Rscript Figure6_confidence_interval_widths.R \
#     data/sim_conf_ints_alpha=0.1_c=1_m=1_method=elbow_n=50_p=10_rank=5_reps=1000_signal_alpha_frac=0.75.RData \
#     --sigmas 0.15 0.25

# # Power: regular
# Rscript Figure4_power.R \
#     data/sim_hypo_tests_alpha=0.1_m=1_method=elbow_n=50_p=10_rank=5_reps=1000.RData
