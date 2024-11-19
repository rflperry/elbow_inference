# Generate a timestamp
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

# Main experiments
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     -n 50 \
#     -p 10 \
#     --rank 5 \
#     --reps 1000 \
#     --sigmas 0.1 0.2 0.35 0.5 \
#     2>&1

# Elbow rule: coverage
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --method elbow \
#     --reps 10000 \
#     --sigmas 0.5 \
#     2>&1

# Elbow rule: CI widths
# fname="simulate_confidence_intervals"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --sigmas 0.15 0.5 \
#     --method elbow \
#     --reps 1000 \
#     2>&1

# Elbow rule: CI widths
# fname="simulate_hypothesis_tests"
# log_file="logs/$fname"_"$timestamp.log"
# Rscript "$fname.R" > "$log_file" \
#     --method elbow \
#     --sigmas 0.1 0.2 0.3 1 2 \
#     --reps 1000 \
#     2>&1