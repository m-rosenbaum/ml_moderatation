################################################################################
#   Title: Main Simulation Script
#   Author: Michael Rosenbaum
#   Notes: adapted from https://github.com/grf-labs/grf/tree/master/experiments/ijmpr
################################################################################
# 1. Set up environment
# 2. Run simulation



#######################################
# 1. Set up working environment
#######################################
# A. Clear workspace
# B. Load packages
# C. Load files


## A. Clear local files
rm(list = ls())

## B. Set up library
library(grf)
library(ggplot2)
library(tidyverse)
library(purrr)

## C. Load files
source("R/config/config.R")
source("R/sim/simulation_functions.R")
source("R/sim/data_generation.R")




#######################################
# 1. Main execution
#######################################
main <- function() {
    # Validate environment
    validate_environment()

    # Run simulations
    results <- expand.grid(
        n_cols = SIM_PARAMS$n_cols,
        n_rows = SIM_PARAMS$n_rows,
        ate = SIM_PARAMS$ate
    ) %>%
    pmap(run_simulation) # Map the cartesian product of results into simulation

    message("Simulations complete")
}

# Execute
main()
