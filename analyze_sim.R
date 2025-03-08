################################################################################
#   Title: Simulating
#
#   Author: Michael Rosenbaum
#
#   Notes: adapted https://github.com/grf-labs/grf/tree/master/experiments/ijmpr
#
################################################################################
# 0. Working environment
# 1. Load data
# 2. Save results


#######################################
# 0. Working environment
#######################################
# 1. Set up working environment
# 2. Load packages
# 3. Helper functions


## 1. Set up working environment
rm(list = ls())


## 2. Set up library
library(ggplot2)
library(tidyverse)
library(stats)
library(stringr)


## 3. Helper functions
extract_metadata <- function(filename) {
    match <- stringr::str_match(basename(filename), "sim_n_col_(\\d+)_n_row_(\\d+)const_(TRUE|FALSE)")
    if (is.na(match[1])) return(NULL)
    return(tibble(
        n_col = as.integer(match[2]),
        n_row = as.integer(match[3]),
        const = as.logical(match[4])
    ))
}


analyze_csv_files <- function(directory) {
    files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
    results <- list()

    for (file in files) {
        metadata <- extract_metadata(file)
        if (is.null(metadata)) next

        data <- read.csv(file)
        if (nrow(data) == 1001) data <- data[-1001, ]
        if ("autoc_p" %in% colnames(data)) {
            autoc_p_values <- data$autoc_p
            below_threshold <- sum(autoc_p_values < 0.05, na.rm = TRUE)
            ks_test <- ks.test(autoc_p_values, "punif")

            results <- bind_rows(results, metadata %>% mutate(
                below_05 = below_threshold,
                ks_p_value = ks_test$p.value
            ))
        }
    }

    return(results)
}

#######################################
#
#######################################
# 1. Run Analysis


results <- analyze_csv_files("G:/My Drive/CAPP/2025_Q2_ECON 41300_Exp/Paper/Sim Results/")
print(results)
