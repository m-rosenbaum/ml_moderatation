################################################################################
#   Title: Data Generation Functions
#   Author: Michael Rosenbaum
################################################################################

source("config.R") 

gen_output_df <- function() {
    # Generates a tibble with appropriately typed and named columns to fill
    # with simulation results
    #
    # Inputs: None
    #
    # Returns: a 0x5 tibble


    df <- tibble(
        seed = integer(),
        autoc_p = double(),
        cate_e = double(),
        cate_se = double(),
        ks_p = double()
    )

    return(df)
}

# Add outcome to the DF with defaults as NA
add_sim_to_df <- function(df,
                            seed = NA_integer_,
                            autoc_p = NA_real_,
                            cate_e = NA_real_,
                            cate_se = NA_real_,
                            ks_p = NA_real_) {
    # Fills the output data frame with the results from this iteration
    # of the simulaion, adding in a correctly typed Null if the column
    # hasn't been implemented yet.
    #
    # Input:
    #   df (tibble): n x 5 set of simulation results
    #   5x value from the simulation with defaults as NUlls
    #
    # Returns: an n+1 x 5 tibble


    # Append the new row to the table
    df <- bind_rows(df, tibble(seed = seed,
                                autoc_p = autoc_p,
                                cate_e = cate_e,
                                cate_se = cate_se,
                                ks_p = ks_p)
                    )
    return(df)
}

gen_data <- function(n_obs, n_cols, seed, constant=TRUE) {
    # Create simulated data to test estimated CATE procedure that can
    # handle complex nonlinearities.
    #
    # Input
    #   n_obs (int): Number of data points to generate
    #   n_cols (int): Number of columns to generate
    #   seed (int): int to use as see dfor dataset creation
    #   constant (bool): If the ATE should be constant or of size 1
    #
    # Returns: tibble of N x M+2 datapoints with treatment information

    # First set seed
    set.seed(seed)

    # Generate uniform random numbers to decide distributions
    random_draws <- runif(n_cols)

    # Xs
    for (i in 1:n_cols) {
        draw <- random_draws[i]

        # Assign a distribution based on the draw value with mean 0 and sd 1
        if (draw < 1/6) {
            col_data <- rnorm(n_obs, mean = 0, sd = 1)
        } else if (draw < 2/6) {
            col_data <- rpois(n_obs, lambda = 1)
        } else if (draw < 3/6) {
            col_data <- runif(n_obs, min = -1, max = 1)
        } else if (draw < 4/6) {
            col_data <- rbeta(n_obs, shape1 = 2, shape2 = 5) * 2 - 1
        } else if (draw < 5/6) {
            col_data <- sample(0:2, n_obs, replace = TRUE)
        } else {
            col_data <- sin(rnorm(n_obs, 0, 1)) + log(abs(rnorm(n_obs, 0, 1)) + 1)
        }

        # Add the column to the Xs using dynamic naming
        if (i == 1) {
            # Initialize tibble with one col of specified length
            xs <- setNames(tibble(col_data), "X1")
        } else {
            xs <- xs %>%
                mutate(!!paste0("X", i) := col_data)
                # !! Evaluates the variable directly
                # := dynamically names for mutate
        }
    }

    # Constant treatment effect, centered at 1, independent of X
    rand <- runif(n_obs,0)
    if (constant) {
        tau <- 1
    } else {
        tau <-  runif(n_obs)*xs$X1*0.1 +
                runif(n_obs)*xs$X2*0.2 +
                runif(n_obs)*xs$X3*0.3 +
                runif(n_obs)*xs$X4*0.4 +
                runif(n_obs)*xs$X5*0.5 +
                1
    }

    # Outcome variable a complicated function of Xs with some noise
    if (n_cols >= 9) {
        outcome <-
            pull(xs, 1) +
            pull(xs, 4) +
            pull(xs, 7) * pull(xs, 2) +
            pull(xs, 9) +
            rnorm(n_obs, 0, 1)
    } else {
        # TODO: Figure out what to do if small Ns
        outcome <- rowSums(xs) + rnorm(n_obs, 0, 1)
    }

    # Create data frame
    df <- tibble(outcome, rand, xs, tau) %>%
        arrange(rand) %>%
        mutate(treatment = if_else(row_number() > (n_obs / 2), 1, 0),
                outcome = if_else(treatment == 1, outcome + tau, outcome)) %>%
        select(-c(rand)) %>%
        select(outcome, treatment, tau, starts_with("X")) # Reorder outcomes

    return(df)
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
    cate.forest <- causal_forest(X[train, ], Y[train], W[train])
    eval.forest <- causal_forest(X.test, Y[test], W[test])
    tau.hat.test <- predict(cate.forest, X.test)$predictions # Extract predictions

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
