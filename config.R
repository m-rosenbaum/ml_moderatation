################################################################################
#   Title: Configuration for Causal Forest Simulation
#   Author: Michael Rosenbaum
################################################################################

# Simulation parameters
SIM_PARAMS <- list(
    n_cols = c(15, 100),
    n_rows = c(250, 500, 1000, 2500, 10000),
    n_iterations = 1000,
    constant_effects = c(FALSE, TRUE),
    train_fraction = 0.6
)

# File paths
PATHS <- list(
    seed_file = "seeds.csv",
    output_dir = "G:/My Drive/CAPP/2025_Q2_ECON 41300_Exp/Paper/",
    results_dir = "3_output/"
)

# GRF parameters
GRF_PARAMS <- list(
    num_trees = 2000,
    min_node_size = 5,
    honesty = TRUE,
    alpha = 0.05
) 