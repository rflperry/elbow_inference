# Generate a timestamp
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

fname="simulate_confidence_intervals"
# fname="simulate_hypothesis_tests"

# Define the log file with the timestamp in its name
log_file="logs/$fname"_"$timestamp.log"

# Run the R script and redirect both stdout and stderr to the log file
Rscript "$fname.R" > "$log_file" \
    -n 50 \
    -p 10 \
    --rank 5 \
    --reps 1000 \
    --sigmas 0.1 0.2 0.35 0.5 \
    2>&1
    # -n 100 \
    # -m 2 \
    # -p 50 \
    # --rank 10 \
    # --signal_alpha_frac 0.75 \