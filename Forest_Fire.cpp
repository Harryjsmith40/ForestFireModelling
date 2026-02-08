#include <iostream>
#include <vector>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <algorithm>

// Constants for cell states
const int EMPTY = 0;
const int TREE = 1;
const int BURNING = 2;

// Function to map global coordinates to a local index
int global_to_local_index(int global_row, int global_col, int local_start_row, int N) {
    return (global_row - local_start_row) * N + global_col;
}

// Function to map a local index to global row and column
void local_to_global_coords(int local_index, int local_start_row, int N, int& global_row, int& global_col) {
    global_row = local_start_row + local_index / N;
    global_col = local_index % N;
}

// Main simulation function
void ForestFireSimulation(int N, double p, int M, const std::string& output_dir) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Ensure buffers are properly initialized
    std::vector<int> send_buffer(N, 0);
    std::vector<int> recv_buffer(N, 0);

    // Check if buffers are properly initialized
    if (send_buffer.empty() || recv_buffer.empty()) {
        std::cerr << "Error: Null buffer pointer on rank " << rank << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < M; ++i) {
        std::fill(send_buffer.begin(), send_buffer.end(), rank);

        int prev_rank = (rank - 1 + size) % size;
        int next_rank = (rank + 1) % size;
        MPI_Sendrecv(send_buffer.data(), N, MPI_INT, next_rank, 0,
                     recv_buffer.data(), N, MPI_INT, prev_rank, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Default parameters
    int N = 100;
    double p = 0.6;
    int M = 50;
    std::string input_file = "";
    std::string output_dir = ".";

    // Parse command-line arguments
    if (argc > 1) {
        try {
            N = std::stoi(argv[1]);
            p = std::stod(argv[2]);
            M = std::stoi(argv[3]);
            if (argc > 4) {
                output_dir = argv[4];
            }
        } catch (const std::invalid_argument& e) {
            if (world_rank == 0) {
                std::cerr << "Error: Invalid command-line arguments. Using default values.\n";
            }
            N = 100;
            p = 0.6;
            M = 50;
            output_dir = ".";
        }
    }

    // Check if output directory is writable
    if (world_rank == 0) {
        std::ofstream test_file(output_dir + "/test_output.txt");
        if (!test_file.is_open()) {
            std::cerr << "Rank 0: Error: Output directory '" << output_dir << "' is not writable. Using current directory.\n";
            output_dir = ".";
        }
        test_file.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::srand(std::time(nullptr) + world_rank);

    // Initialize variables for performance measurement
    double total_steps = 0.0;
    double total_time = 0.0;
    bool fire_reached_bottom = false;
    bool local_fire_reached_bottom = false;

    std::chrono::time_point<std::chrono::steady_clock> start_time;
    std::chrono::time_point<std::chrono::steady_clock> end_time;
    std::chrono::duration<double> elapsed_seconds;

    double avg_steps = 0.0;
    double avg_time = 0.0;
    std::string filename;
    std::ofstream outfile;
    std::string timestamp;
    int runs_reached_bottom = 0;

    // Main simulation loop
    for (int run = 0; run < M; ++run) {
        std::stringstream ss;
        auto in_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
        timestamp = ss.str();

        // Determine the rows assigned to each process
        int my_start_row;
        int my_end_row;
        int rows_per_process = N / world_size;
        int remaining_rows = N % world_size;

        if (world_rank < remaining_rows) {
            my_start_row = world_rank * (rows_per_process + 1);
            my_end_row = my_start_row + rows_per_process + 1;
        } else {
            my_start_row = world_rank * rows_per_process + remaining_rows;
            my_end_row = my_start_row + rows_per_process;
        }

        int local_rows = my_end_row - my_start_row;
        int local_grid_size = local_rows * N;
        std::vector<int> local_grid(local_grid_size);

        // Initialize the grid
        if (input_file.empty()) {
            if (world_rank == 0) {
                std::vector<int> full_grid(N * N);
                for (int i = 0; i < N * N; ++i) {
                    full_grid[i] = (rand() < p * RAND_MAX) ? TREE : EMPTY;
                }
                for (int c = 0; c < N; ++c) {
                    if (full_grid[c] == TREE) {
                        full_grid[c] = BURNING;
                    }
                }
                int sendcounts[world_size];
                int displs[world_size];
                int current_disp = 0;
                for (int i = 0; i < world_size; ++i) {
                    int proc_rows = N / world_size;
                    int proc_remaining_rows = N % world_size;
                    if (i < proc_remaining_rows) {
                        proc_rows++;
                    }
                    sendcounts[i] = proc_rows * N;
                    displs[i] = current_disp;
                    current_disp += sendcounts[i];
                }
                MPI_Scatterv(full_grid.data(), sendcounts, displs, MPI_INT,
                             local_grid.data(), local_grid_size, MPI_INT,
                             0, MPI_COMM_WORLD);
            } else {
                int sendcounts[world_size];
                int displs[world_size];
                int current_disp = 0;
                for (int i = 0; i < world_size; ++i) {
                    int proc_rows = N / world_size;
                    int proc_remaining_rows = N % world_size;
                    if (i < proc_remaining_rows) {
                        proc_rows++;
                    }
                    sendcounts[i] = proc_rows * N;
                    displs[i] = current_disp;
                    current_disp += sendcounts[i];
                }
                MPI_Scatterv(NULL, sendcounts, displs, MPI_INT,
                             local_grid.data(), local_grid_size, MPI_INT,
                             0, MPI_COMM_WORLD);
            }
        } else {
            if (world_rank == 0) {
                std::ifstream file(input_file);
                if (!file.is_open()) {
                    std::cerr << "Rank 0: Error opening input file.\n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                std::vector<int> full_grid(N * N);
                for (int r = 0; r < N; ++r) {
                    for (int c = 0; c < N; ++c) {
                        if (!(file >> full_grid[r * N + c])) {
                            std::cerr << "Rank 0: Error reading input file. Check that the file has the correct format (N x N integers).\n";
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                    }
                }
                for (int c = 0; c < N; ++c) {
                    if (full_grid[c] == TREE) {
                        full_grid[c] = BURNING;
                    }
                }
                file.close();
                int sendcounts[world_size];
                int displs[world_size];
                int current_disp = 0;
                for (int i = 0; i < world_size; ++i) {
                    int proc_rows = N / world_size;
                    int proc_remaining_rows = N % world_size;
                    if (i < proc_remaining_rows) {
                        proc_rows++;
                    }
                    sendcounts[i] = proc_rows * N;
                    displs[i] = current_disp;
                    current_disp += sendcounts[i];
                }
                MPI_Scatterv(full_grid.data(), sendcounts, displs, MPI_INT,
                             local_grid.data(), local_grid_size, MPI_INT,
                             0, MPI_COMM_WORLD);
            } else {
                int sendcounts[world_size];
                int displs[world_size];
                int current_disp = 0;
                for (int i = 0; i < world_size; ++i) {
                    int proc_rows = N / world_size;
                    int proc_remaining_rows = N % world_size;
                    if (i < proc_remaining_rows) {
                        proc_rows++;
                    }
                    sendcounts[i] = proc_rows * N;
                    displs[i] = current_disp;
                    current_disp += sendcounts[i];
                }
                MPI_Scatterv(NULL, sendcounts, displs, MPI_INT,
                             local_grid.data(), local_grid_size, MPI_INT,
                             0, MPI_COMM_WORLD);
            }
        }

        // Simulation loop for each run
        bool simulation_done = false;
        local_fire_reached_bottom = false;
        std::vector<int> bottom_row_below(N);
        int global_burning_count;
        int steps = 0;
        start_time = std::chrono::steady_clock::now();
        while (!simulation_done) {
            std::vector<int> next_local_grid = local_grid;
            steps++;
            std::copy(local_grid.begin(), local_grid.end(), next_local_grid.begin());

            std::vector<int> top_row_above(N);
            std::vector<int> bottom_row_below(N);

            if (world_rank > 0) {
                MPI_Sendrecv(local_grid.data(), N, MPI_INT, world_rank - 1, 0,
                             top_row_above.data(), N, MPI_INT, world_rank - 1, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (world_rank < world_size - 1) {
                MPI_Sendrecv(&local_grid[(local_rows - 1) * N], N, MPI_INT, world_rank + 1, 0,
                             bottom_row_below.data(), N, MPI_INT, world_rank + 1, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Update the grid based on the fire spread rules
            for (int r = 0; r < local_rows; ++r) {
                for (int c = 0; c < N; ++c) {
                    int global_row = my_start_row + r;
                    int index = r * N + c;

                    if (local_grid[index] == TREE) {
                        bool neighbor_burning = false;
                        if (c > 0 && local_grid[index - 1] == BURNING) neighbor_burning = true;
                        if (c < N - 1 && local_grid[index + 1] == BURNING) neighbor_burning = true;
                        if (r > 0) {
                            if (local_grid[index - N] == BURNING) neighbor_burning = true;
                        } else if (world_rank > 0 && top_row_above[c] == BURNING) {
                            neighbor_burning = true;
                        }

                        if (r < local_rows - 1) {
                            if (local_grid[index + N] == BURNING) neighbor_burning = true;
                        } else if (world_rank < world_size - 1 && bottom_row_below[c] == BURNING) {
                            neighbor_burning = true;
                        }

                        if (neighbor_burning) {
                            next_local_grid[index] = BURNING;
                        }
                    } else if (local_grid[index] == BURNING) {
                        next_local_grid[index] = EMPTY;
                    }
                }
            }
            local_grid = next_local_grid;

            // Count the number of burning cells
            int local_burning_count = 0;
            for (int i = 0; i < local_grid_size; ++i) {
                if (local_grid[i] == BURNING) {
                    local_burning_count++;
                }
            }
            MPI_Allreduce(&local_burning_count, &global_burning_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            // Check if the simulation is done
            if (global_burning_count == 0) {
                simulation_done = true;
            }

            // Check if the fire reached the bottom row
            if (world_rank == world_size - 1) {
                for (int c = 0; c < N; ++c) {
                    if (local_grid[(local_rows - 1) * N + c] == BURNING) {
                        local_fire_reached_bottom = true;
                        break;
                    }
                }
            }
            MPI_Allreduce(&local_fire_reached_bottom, &fire_reached_bottom, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        }
        end_time = std::chrono::steady_clock::now();
        elapsed_seconds = end_time - start_time;
        double elapsed_time = elapsed_seconds.count();

        // Accumulate performance data
        total_steps += steps;
        total_time += elapsed_time;
        if (fire_reached_bottom) {
            runs_reached_bottom++;
        }
        local_fire_reached_bottom = false;
    }

    // Output results
    if (world_rank == 0) {
        std::stringstream filename_stream;
        filename_stream << output_dir << "/results_N" << N << "_p" << std::fixed << std::setprecision(2) << p
                        << "_M" << M << "_" << timestamp << ".txt";
        std::string filename = filename_stream.str();
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            std::cerr << "Rank 0: Error opening output file: " << filename << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        double avg_steps = total_steps / M;
        double avg_time = total_time / M;
        outfile << "average_steps = " << std::fixed << std::setprecision(6) << avg_steps << "\n";
        outfile << "average_elapsed_time = " << std::fixed << std::setprecision(6) << avg_time << "\n";
        outfile << "fire_reached_bottom = " << runs_reached_bottom << "/" << M << "\n";

        outfile << "performance_data:\n";
        outfile << "total_steps = " << total_steps << "\n";
        outfile << "total_time = " << total_time << "\n";
        outfile << "world_size = " << world_size << "\n";

        outfile.close();
        std::cout << "Rank 0: Wrote output to " << filename << "\n";

        // Output the average steps to standard output for the bash script to capture
        std::cout << "average_steps = " << std::fixed << std::setprecision(6) << avg_steps << "\n";
    }
    MPI_Finalize();
    return 0;
}
