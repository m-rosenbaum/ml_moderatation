################################################################################
#   Title: Core Simulation Functions
#   Author: Michael Rosenbaum
################################################################################


#######################################
# 1. Simulation function
#######################################
run_simulation <- function(n_cols, n_rows, ate) {
    # Get seed and validate
    seed <- get_seed(n_cols, n_rows, ate)
    output <- gen_output_df()
    
    # Generate synthetic data
    data <- gen_data(n_obs = n_rows, n_cols = n_cols, seed = seed, ate = ate)
    
    # Calculate true correlations if heterogeneous effects
    if (ate != 0) {
        true_corrs <- calc_corrs(data$tau, data %>% select(X1:X5))
    }
    
    # Prepare data for simulation
    Y <- data$outcome
    W <- data$treatment
    X <- data[, -c(1, 2, 3)]  # Omit outcome, treatment, and tau
    
    # Run iterations
    for (i in (seed + 1):(seed + SIM_PARAMS$n_iterations)) {
        output <- sim_iter_cate(Y, X, W, output, i)
        log_progress(i, seed, n_cols, n_rows, ate)
    }
    
    # Save output to file
    save_results(output, n_cols, n_rows, ate)
    return(output)
}

sim_iter_cate <- function(Y, X, W, output, seed) {
    # Run a single simulation to estimate CATEs using the Wager procedure.
    #
    # Input:
    #   Y: Vector of outcomes
    #   X: Vector of Xs
    #   W: Vector of treatment assignments
    #   output: a n-1 x 5 tibble to add a simulation result to.
    #   seed: Seed used for this iteration of the simulation
    #
    # Returns: A N x 5 tibble documenting all current iterations of the
    #          simulation.

    ## Create sim prior to sample splitting
    set.seed(seed)

    # Split data into a train and test sample.
    train <- sample(nrow(X), 0.6 * nrow(X))
    test <- -train
    X.test <- X[test, ]

    ## 1. Estimate CATE forest
    # Fit a CATE function on training data.
    cate.forest <- causal_forest(X[train, ], Y[train], W[train], 
        W.hat = GRF_PARAMS$w_hat, 
        num.trees = GRF_PARAMS$num_trees,
        min.node.size = GRF_PARAMS$min_node_size,
        alpha = GRF_PARAMS$alpha,
        honesty = GRF_PARAMS$honesty
        )
    eval.forest <- causal_forest(X.test, Y[test], W[test], 
        W.hat = GRF_PARAMS$w_hat, 
        num.trees = GRF_PARAMS$num_trees,
        min.node.size = GRF_PARAMS$min_node_size,
        alpha = GRF_PARAMS$alpha,
        honesty = GRF_PARAMS$honesty
        )
    tau.hat.test <- predict(cate.forest, X.test)$predictions

    # *** Evaluate heterogeneity via TOC/AUTOC ***
    # Use eval.forest to form a doubly robust estimate of TOC/AUTOC.
    rate.cate <- rank_average_treatment_effect(
        eval.forest,
        tau.hat.test,
        q = seq(0.05, 1, length.out = 100)
    )

    # Get a 2-sided p-value Pr(>|t|) for RATE = 0 using a t-value.
    cate_e <- rate.cate$estimate
    cate_se <- rate.cate$std.err
    autoc_p <- 2 * pnorm(-abs(rate.cate$estimate / rate.cate$std.err))

    # Put in the errors in the sampling
    output <- add_sim_to_df(output,
                            seed = seed,
                            cate_e = cate_e,
                            cate_se = cate_se,
                            autoc_p = autoc_p)
    return(output)
}

#######################################
# 2. Helper functions for reproducibility and logging
#######################################
read_and_validate_seeds <- function() {
    # Read seeds file
    seeds <- read.csv(PATHS$seed_file)
    
    # Check columns correct
    required_cols <- c("n_cols", "n_rows", "ate", "seed")
    missing_cols <- setdiff(required_cols, names(seeds))
    if (length(missing_cols) > 0) {
        stop("Missing required columns in seeds file: ", 
             paste(missing_cols, collapse = ", "))
    }
    
    return(seeds)
}

get_seed <- function(n_cols, n_rows, ate) {
    seeds <- read_and_validate_seeds()
    
    seed <- seeds %>%
        filter(n_cols == !!n_cols, 
               n_rows == !!n_rows, 
               ate == !!ate) %>%
        pull(seed)
    
    if(length(seed) != 1) {
        stop("Could not find unique seed value for selected parameters")
    }
    
    return(seed)
}

log_progress <- function(current_iter, seed, n_cols, n_rows, ate) {
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
            " | ATE = ", ate,
            " | Time for 50 iterations: ", 
            if(is.character(time_diff)) time_diff else sprintf("%.2f secs", as.numeric(time_diff)),
            "\n"
        )
    }
}

save_results <- function(output, n_cols, n_rows, ate) {
    filename <- file.path(
        PATHS$output_dir,
        PATHS$results_dir,
        paste0(
            "sim_n_col_", n_cols, 
            "_n_row_", n_rows, 
            "_const_", ate, 
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

calc_corrs <- function(tau, columns) {
    #  Compute correlations between the outcome and the variables.
    #
    #  Input:
    #   - tau (vector): N x N vector of treatment effect
    #   - columns (list): List of N x N vectors to collect treatment effects for
    #
    #  Returns: Tibble of correlations
    return(sapply(columns, function(col) cor(tau, col)))
}
