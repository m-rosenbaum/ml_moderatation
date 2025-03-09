library(testthat)
library(dplyr)


#######################################
# 1. calc_mde
#######################################
test_that("calc_mde estimates MDEs correctly", {
    set.seed(1)
    outcome <- rnorm(100, 0, 1)

    mde <- calc_mde(outcome)
    expect_equal(mde, 0.505, tolerance=0.005) # From Stata
})


#######################################
# 2. gen_data
#######################################
test_that("gen_data produces correct number of taus below 0", {
  # Test with larger sample size for stable estimates
  seed <- 1
  n_obs <- 10000
  ate <- 0.05
  n_cols <- 15

  result <- gen_data(n_obs = n_obs, n_cols = n_cols, seed = seed, ate = ate)

  # Calculate proportion of negative tau_hat
  tau_hat <- result$tau[result$treatment==1]
  prop_negative <- mean(tau_hat < 0)

  # Test proportion of negative values is approximately ate within 0.5%
  expect_equal(prop_negative, ate, tolerance = 0.01)
})

test_that("gen_data produces ATEs that match epected", {
  # Test with larger sample size for stable estimates
  seed <- 1
  n_obs <- 1000
  ate <- 0.3
  n_cols <- 15

  # Select
  result <- gen_data(n_obs = n_obs, n_cols = n_cols, seed = seed, ate = ate)

  # Test treatment effect
  model <- lm(outcome ~ treatment, data = result)
  treatment_effect <- unname(coef(model)["treatment"])
  expect_equal(treatment_effect, mean(result$tau[result$treatment==1]), tolerance = 0.01)
})


test_that("gen_data maintains statistical properties across different parameters", {
  # Test with different parameter combinations
  test_cases <- list(
    list(n_obs = 500,  n_cols=15, seed=2, ate = 0),
    list(n_obs = 1000,  n_cols=15, seed=3, ate = 0.01),
    list(n_obs = 2500,  n_cols=15, seed=4, ate = 0),
    list(n_obs = 5000,  n_cols=15, seed=5, ate = 0.1),
    list(n_obs = 10000, n_cols=15, seed=6, ate = 0.25)
  )

  for(case in test_cases) {
    result <- do.call(gen_data,
                     c(case))

    tau_hat <- result$tau[result$treatment==1]
    prop_negative <- mean(tau_hat < 0)

    # Test proportion of negative values
    expect_equal(prop_negative, case$ate, tolerance = 0.05)

    # Test treatment effect
    model <- lm(outcome ~ treatment, data = result)
    treatment_effect <- unname(coef(model)["treatment"])
    expect_equal(treatment_effect, mean(result$tau[result$treatment==1]), tolerance = 0.01)
  }
})

test_that("gen_data statistical properties are reproducible", {
  n_obs <- 10000
  ate <- 0.3
  seed <- 5

  result1 <- gen_data(n_obs = n_obs, n_cols = 15, seed=seed, ate = ate)
  tau_hat1 <- result1$tau
  prop_negative1 <- mean(tau_hat1 < 0)
  model1 <- lm(outcome ~ treatment, data = result1)

  result2 <- gen_data(n_obs = n_obs, n_cols = 15, seed=seed, ate = ate)
  tau_hat2 <- result2$tau
  prop_negative2 <- mean(tau_hat2 < 0)
  model2 <- lm(outcome ~ treatment, data = result2)

  expect_equal(prop_negative1, prop_negative2)
  expect_equal(coef(model1), coef(model2))
})

test_that("gen_data creates correct dimensions", {
  n_cols <- 15
  seed <- 7
  ate <- 0

  # Basic case
  result <- gen_data(n_obs = 100, n_cols = 15, seed=seed, ate=ate)
  expect_equal(nrow(result), 100)
  expect_equal(ncol(result), 18) # 3 extra columns (treatment, Y, tau)

  # Larger case
  result_large <- gen_data(n_obs = 1000, n_cols = 50, seed=seed, ate=ate)
  expect_equal(nrow(result_large), 1000)
  expect_equal(ncol(result_large), 53) # 3 extra columns (treatment, Y, tau)
})

test_that("gen_data column names are correct", {
  result <- gen_data(n_obs = 100, n_cols = 20, seed=8, ate=0)
  expected_names <- paste0("X", 1:20)
  expect_equal(colnames(result[,4:ncol(result)]), expected_names)
})

test_that("gen_data generates numeric data", {
  result <- gen_data(n_obs = 200, seed=9, n_cols = 15, ate=0)
  expect_true(all(sapply(result, is.numeric)))
})

test_that("gen_data handles edge cases", {
  # Minimum allowed values
  result_min <- gen_data(n_obs = 3, n_cols = 11, seed=10, ate=0) # 2 needed for DF of 2 in qt()
  expect_equal(nrow(result_min), 3)
  expect_equal(ncol(result_min), 14)

  # Large values
  result_large <- gen_data(n_obs = 10000, n_cols = 100, seed=11, ate=0)
  expect_equal(nrow(result_large), 10000)
  expect_equal(ncol(result_large), 103)
})

test_that("gen_data produces expected range of values", {
  result <- gen_data(n_obs = 500, n_cols = 20, seed=12, ate=0)

  # Check if values are within expected range
  expect_true(all(result >= -10))
  expect_true(all(result <= 10))

  # Check if there's variation in the data
  expect_true(any(result != 0))
})

test_that("gen_data produces reproducible results with set.seed", {
  seed=13
  result1 <- gen_data(n_obs = 100, n_cols = 15, seed=13, ate=0)
  result2 <- gen_data(n_obs = 100, n_cols = 15, seed=13, ate=0)

  expect_equal(result1, result2)
})
