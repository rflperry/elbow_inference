# Generate a timestamp
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

# #
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 	ZG Rule
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #

# fname="simulate_hypothesis_test"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     -n 50 \
#     -p 10 \
#     --rank 5 \
#     --reps 1000 \
#     --sigmas 0.1 0.2 0.4 0.5 0.7 1 \
#     2>&1


# # Power, higher dimensional
# fname="simulate_hypothesis_tests"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     -n 100 \
#     -p 50 \
#     --reps 1000 \
#     2>&1

# # Power, low rank
# fname="simulate_hypothesis_tests"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --rank 2 \
#     --reps 1000 \
#     2>&1

# # Power, high rank
# fname="simulate_hypothesis_tests"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --rank 10 \
#     --reps 1000 \
#     2>&1

# # CI widths, regular
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.1 0.2 \
#     --reps 1000 \
#     2>&1

# # CI, higher dimensional
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.05 0.1 0.2 \
#     -n 100 \
#     -p 50 \
#     --reps 1000 \
#     2>&1

# # CI, lower rank
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.05 0.1 0.2 \
#     --rank 2 \
#     --reps 1000 \
#     2>&1

# # CI, higher rank
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.05 0.1 0.2 \
#     --rank 10 \
#     --reps 1000 \
#     2>&1

# # CI, many alphas
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.1 \
#     --alpha 0.1 \
#     --reps 10000 \
#     2>&1
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --alpha 0.3 \
#     --reps 10000 \
#     2>&1
# Rscript "$fname.R" > "$log_file" \
#     --alpha 0.5 \
#     --reps 10000 \
#     2>&1
# Rscript "$fname.R" > "$log_file" \
#     --alpha 0.7 \
#     --reps 10000 \
#     2>&1
# Rscript "$fname.R" > "$log_file" \
#     --alpha 0.9 \
#     --reps 10000 \
#     2>&1

# Hypothesis tests, but with estimated variance
# fname="simulate_hypothesis_tests"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     -n 50 \
#     -p 10 \
#     --var_est \
#     --rank 5 \
#     --reps 1000 \
#     --sigmas 0.1 0.2 0.4 0.5 0.7 1 \
#     2>&1


# CI, many alphas but with estimated variance
fname="simulate_confidence_intervals"
log_file="logs/$fname"_"$timestamp.log"

Rscript "$fname.R" > "$log_file" \
    --sigmas 0.1 0.2 \
    --alpha 0.1 \
    --var_est \
    --reps 10000 \
    2>&1

# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.1 \
#     --alpha 0.3 \
#     --reps 10000 \
#     --var_est \
#     2>&1

# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.1 \
#     --alpha 0.5 \
#     --reps 10000 \
#     --var_est \
#     2>&1

# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.1 \
#     --alpha 0.7 \
#     --reps 10000 \
#     --var_est \
#     2>&1

# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.1 \
#     --alpha 0.9 \
#     --reps 10000 \
#     --var_est \
#     2>&1

# #
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 	Discrete elbow rule
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #


# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --method elbow \
#     --reps 2000 \
#     --sigmas 0.5 \
#     2>&1

# Elbow rule: CI widths
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.15 0.25 \
#     --method elbow \
#     --reps 1000 \
#     2>&1

# Elbow rule: power
# fname="simulate_hypothesis_tests"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --method elbow \
#     --sigmas 0.1 0.2 0.3 1 2 \
#     --reps 1000 \
#     2>&1