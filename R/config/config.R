################################################################################
#   Title: Configuration for Causal Forest Simulation
#   Author: Michael Rosenbaum
################################################################################

# Simulation parameters
SIM_PARAMS <- list(
    n_cols = c(15, 100),
    n_rows = c(500, 1000, 2500, 5000, 10000), 
    n_iterations = 100,
    ate = c(0, 0.05, 0.25),
    train_fraction = 0.5
)

# File paths
PATHS <- list(
    seed_file = "data/seeds.csv",
    output_dir = "G:/My Drive/CAPP/2025_Q2_ECON 41300_Exp/Paper/",
    results_dir = "output/"
)

# GRF parameters
GRF_PARAMS <- list(
    num_trees = 2000,
    min_node_size = 5,
    honesty = TRUE,
    alpha = 0.05,
    w_hat = 0.5 # Experiment
) 