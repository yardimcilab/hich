#!/bin/bash

# Script name
slurm_script="setup_slurm.sh"


# Default number of cores
num_cores=24

# Copy the original arguments
original_args=("$@")

# Parse arguments for --cores
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores)
            num_cores=$2
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Function to add SBATCH -c line at the start of the file
add_sbatch_line_at_start() {
    temp_file=$(mktemp)
    echo "#SBATCH -c ${num_cores}" > "$temp_file"
    cat "$slurm_script" >> "$temp_file"
    mv "$temp_file" "$slurm_script"
}

# Check if the SLURM script has the SBATCH -c or equivalent line
if grep -qE "^#SBATCH[[:space:]]+(-c|--cores)" "$slurm_script"; then
    # Replace the line
    sed -i "/^#SBATCH[[:space:]]\+\(-c\|--cores\)/c\#SBATCH -c ${num_cores}" "$slurm_script"
else
    # Add the line at the start
    add_sbatch_line_at_start
fi

# Call sbatch with all original arguments
sbatch "$slurm_script" "${original_args[@]}"

