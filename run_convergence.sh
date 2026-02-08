#!/bin/bash

#SBATCH --job-name=forest_fire_blue_crystal  # Job name
#SBATCH --partition=teach_cpu  # Partition to submit to
#SBATCH --nodes=4  # Number of nodes
#SBATCH --ntasks-per-node=16  # Number of tasks per node
#SBATCH --cpus-per-task=1  # Number of CPUs per task
#SBATCH --time=01:00:00  # Maximum runtime
#SBATCH --output=convergence_%j.out  # Standard output file
#SBATCH --error=convergence_%j.err  # Standard error file
#SBATCH --mem-per-cpu=100M  # Memory per CPU

# Simulation parameters
N_value=100
p_values=(0.2 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8)
M_values=(10 20 30 40 50 60 70 80 90 100)
num_tasks=4
output_dir="convergence_results"
max_retries=3

# Create output directory if it doesn't exist
mkdir -p $output_dir
if [ $? -ne 0 ]; then
  echo "Error: Failed to create output directory '$output_dir'."
  exit 1
fi

# Check if the executable exists and is executable
if [ ! -x ./Forest_Fire ]; then
  echo "Error: Forest_Fire executable not found or not executable!"
  exit 1
fi

# Check if the output directory is writable
if [ ! -d $output_dir ] || [ ! -w $output_dir ]; then
  echo "Error: Output directory '$output_dir' is not writable or does not exist."
  exit 1
fi

# Initialize CSV file for results
csv_file="$output_dir/convergence_results.csv"
echo "M,p,avg_steps" > $csv_file

# Function to run a single simulation
run_simulation() {
  local M=$1
  local p=$2
  local retries=0
  local success=0

  while [ $retries -lt $max_retries ]; do
    output=$(srun --mpi=pmi2 -n $num_tasks ./Forest_Fire $N_value $p $M $output_dir 2>&1)
    if [ $? -eq 0 ]; then
      success=1
      break
    else
      echo "Error: srun failed for M=$M, p=$p on attempt $((retries+1))" >> $csv_file
      echo "$output" >> $csv_file
      retries=$((retries+1))
      sleep 5  # Wait before retrying
    fi
  done

  if [ $success -eq 0 ]; then
    echo "Error: srun failed for M=$M, p=$p after $max_retries attempts" >> $csv_file
    return
  fi

  steps=$(echo "$output" | grep -oP 'average_steps\s*=\s*\K[0-9.]+')
  if [ -z "$steps" ]; then
    echo "Error: Failed to extract steps for M=$M, p=$p" >> $csv_file
    return
  fi
  echo "$M,$p,$steps" >> $csv_file
}

# Run simulations for each combination of M and p values
for M in "${M_values[@]}"; do
  for p in "${p_values[@]}"; do
    run_simulation $M $p
  done
done

echo "All simulations completed."
