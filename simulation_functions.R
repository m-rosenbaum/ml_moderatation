################################################################################
#   Title: Core Simulation Functions
#   Author: Michael Rosenbaum
################################################################################


#######################################
# 1. Simulation function
#######################################
run_simulation <- function(n_cols, n_rows, constant) {
    # Get seed and validate
    seed <- get_seed(n_cols, n_rows, constant)
    output <- gen_output_df()
    
    # Generate synthetic data
    data <- gen_data(n_obs = n_rows, n_cols = n_cols, seed = seed, constant = constant)
    
    # Calculate true correlations if heterogeneous effects
    if (!constant) {
        true_corrs <- calc_corrs(data$tau, data %>% select(X1:X5))
    }
    
    # Prepare data for simulation
    Y <- data$outcome
    W <- data$treatment
    X <- data[, -c(1, 2, 3)]  # Omit outcome, treatment, and tau
    
    # Run iterations
    for (i in (seed + 1):(seed + SIM_PARAMS$n_iterations)) {
        output <- sim_iter_cate(Y, X, W, output, i)
        log_progress(i, seed, n_cols, n_rows, constant)
    }
    
    # Save output to file
    save_results(output, n_cols, n_rows, constant)
    return(output)
}


#######################################
# 2. Helper functions for reproducibility and logging
#######################################
read_and_validate_seeds <- function() {
    # Read seeds file
    seeds <- read.csv(PATHS$seed_file)
    
    # Check columns correct
    required_cols <- c("n_cols", "n_rows", "constant", "seed")
    missing_cols <- setdiff(required_cols, names(seeds))
    if (length(missing_cols) > 0) {
        stop("Missing required columns in seeds file: ", 
             paste(missing_cols, collapse = ", "))
    }
    
    return(seeds)
}

get_seed <- function(n_cols, n_rows, constant) {
    seeds <- read_and_validate_seeds()
    
    seed <- seeds %>%
        filter(n_cols == !!n_cols, 
               n_rows == !!n_rows, 
               constant == !!constant) %>%
        pull(seed)
    
    if(length(seed) != 1) {
        stop("Could not find unique seed value for selected parameters")
    }
    
    return(seed)
}

log_progress <- function(current_iter, seed, n_cols, n_rows, constant) {
    # Reset to start at 0 
    iterations_completed <- current_iter - seed
    
    # Print a note on iteration 50
    if (iterations_completed %% 50 == 0) {
        current_time <- Sys.time()
        
        # Initialize or calculate time difference
        # TODO: Determine if messing with global environment is smart
        if (!exists("last_log_time", envir = .GlobalEnv)) {
            assign("last_log_time", current_time, envir = .GlobalEnv)
            time_diff <- "First batch"
        } else {
            time_diff <- difftime(current_time, get("last_log_time", envir = .GlobalEnv), units = "secs")
            assign("last_log_time", current_time, envir = .GlobalEnv)
        }
        
        # Printo progress
        cat(iterations_completed,
            " iterations done for N_COL = ", n_cols,
            " | N_ROW = ", n_rows, 
            " | Const ATE = ", constant,
            " | Time for 50 iterations: ", 
            if(is.character(time_diff)) time_diff else sprintf("%.2f secs", as.numeric(time_diff)),
            "\n"
        )
    }
}

save_results <- function(output, n_cols, n_rows, constant) {
    filename <- file.path(
        PATHS$output_dir,
        paste0(
            "sim_n_col_", n_cols, 
            "_n_row_", n_rows, 
            "_const_", constant, 
            ".csv"
        )
    )
    write.csv(output, filename)
} 

validate_environment <- function() {
    if (!file.exists(PATHS$seed_file)) {
        stop("Seed file not found")
    }
    if (!dir.exists(PATHS$output_dir)) {
        stop("Output directory not found")
    }
}