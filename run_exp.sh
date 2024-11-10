# Generate a timestamp
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

fname="Figure4_zg_confidence_interval_widths"

# Define the log file with the timestamp in its name
log_file="logs/$fname"_"$timestamp.log"

# Run the R script and redirect both stdout and stderr to the log file
Rscript "$fname.R" > "$log_file" \
    -n 100 \
    -p 50 \
    --rank 10 \
    2>&1