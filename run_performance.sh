#!/bin/bash
#SBATCH --job-name=forest_fire_performance
#SBATCH --partition=teach_cpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=performance_%j.out
#SBATCH --error=performance_%j.err
#SBATCH --mem-per-cpu=100M

# Parameters for the simulation
p=0.6
M=50
N_values=(50 100 500)
nodes_list=(1 2 4)
tasks_per_node_list=(1 2 4 8 16)
output_dir="performance_results"

# Create output directory if it doesn't exist
mkdir -p $output_dir

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

# Initialize CSV file for storing performance results
csv_file="$output_dir/performance_results.csv"
echo "N,nodes,tasks_per_node,avg_time" > $csv_file

# Function to run a single simulation
run_simulation() {
  local N=$1
  local nodes=$2
  local tasks_per_node=$3
  local num_processes=$((nodes * tasks_per_node))

  # Skip configurations where the number of processes exceeds the grid size
  if [ $num_processes -gt $N ]; then
    echo "Skipping configuration: num_processes=$num_processes exceeds grid size N=$N"
    return
  fi

  echo "Running performance simulation for N=$N, p=$p, M=$M, nodes=$nodes, tasks_per_node=$tasks_per_node, num_processes=$num_processes"

  # Measure the start time
  start_time=$(date +%s%N)
  # Run the simulation using srun
  if ! srun --mpi=pmi2 -N $nodes --ntasks-per-node=$tasks_per_node ./Forest_Fire $N $p $M $output_dir; then
    echo "Error: srun failed for N=$N, nodes=$nodes, tasks_per_node=$tasks_per_node" >> $csv_file
    echo "Error: srun failed for N=$N, nodes=$nodes, tasks_per_node=$tasks_per_node" >&2
    return
  fi
  # Measure the end time
  end_time=$(date +%s%N)
  elapsed_time=$((end_time - start_time))

  # Calculate the average time per simulation
  avg_time=$(echo "scale=6; $elapsed_time / ($M * 1000000000)" | bc)

  # Append the results to the CSV file
  echo "$N,$nodes,$tasks_per_node,$avg_time" >> $csv_file
}

# Run simulations for each combination of N, nodes, and tasks_per_node values
for N in "${N_values[@]}"; do
  for nodes in "${nodes_list[@]}"; do
    for tasks_per_node in "${tasks_per_node_list[@]}"; do
      run_simulation $N $nodes $tasks_per_node
    done
  done
done

echo "Performance simulations completed."