# Moderation and Causal Forests
This repo contains code and output to simulate how causal forests perform under various conditions to understand when they can be usable as a tool to detect heterogeneity.

This work was conducted as part of Michael Rosenbaum's work for ECON 41300 A Course on Experimental Economics at the University of Chicago, taught by John List.

## 0 | Repo Structure

## 1 | Research Question

## 2 | Analysis Plan

### Data Generation
I generate data according to the following procedure:

- **Features (Xs)** are generated randomly that use an equal mix of 6 different types of variables:
    - Normal (mean=0, sd=1)
    - Poisson (lambda=1)
    - Uniform (-1 to 1)
    - Transformed Beta (scaled to -1 to 1)
    - Discrete values (0,1,2)
    - Non-linear combination: sin(normal) + log(|normal| + 1)
- **Treatment effect ($\tau$)** are generated randomly as well:
    - Depending on what I test, $\tau$ is fixed at 1 for all observations, or is estimated with 50% mediation from 5 Xs with either 5% or 25% of the taus below 0.
- **Outcome Generation (Y)**: 50% of the variation is explained by features. To do so, I combines a complex linear set of features (2 of which are included in estimating heterogeneous treatment effects) at 10% each and then 50% draws from a normal distribution, with a mean of 0 and SD 1. 
- **Random assignment** is done as simple random assignment with 50% treatment probability and is sorted by random number with a separate seed to ensure reproducibility. 

### Testing
I test two procedures:
- TODO: List's bagging procedure with num_cols
- Wager (2021)
- TODO: Haushofer et. al. (2022) https://www.nber.org/papers/w30138

I test the RATE test and KS test aproaches to detect heterogeneity using simulation over sample size. [Yadlowsky et. al.](https://arxiv.org/abs/2111.07966) estimate the power of the RATE test (2021):
> The statistical power of the hypothesis test described in Section A.4 depends on the strength
of the heterogeneous treatment effects, the sample size, the choice of weighting function, and
the choice of score.

I will hold the weighting function and choice of scoring constant -- using the defaults in the RATE package and then test varying sample size and how the heterogeneity is distributed. I also look at the number of columns. This leaves three types of tests:

1. Number of columns (`n_cols`): 15 or 1500
2. Heterogeneity (`constant`): Constant, 5% dif sign, 25% dif sign
3. Number of observations (`n_obs`): 500, 1000, 2500, 5000, 1000

#### Defining columns
I draw either 15 or 100 columns based on an observation in [Rehill](https://arxiv.org/pdf/2404.13356) that suggests most applied work uses fewer columns than the expected minimum in the grf package (2024):
> The number of variables used is also generally small with the median number being 17. Figure 4 shows the distribution of covariate counts.

Rehill also notes that this may have effects on out-of-sample prediction if the `grf` package is used:
> As an aside, the approach to sampling variables for use in the individual trees by default draws $\sqrt{p} + 20$ variables in the grf causal forest. This means for most papers in our sample, the causal forest is not so much a kind of modified random forest, but instead a modified bagged trees ensemble where different trees are fit on different subsets of the data but the full set of variables (Breiman 2001, 1996). This is not necessarily a problem, but it is worth keeping in mind that the decorrelation of trees via random variable selection is the key insight of the random forest and when this parameter is left at its default value, we may be hurting the out-of-sample fit of the ensemble.

Therefore, I include a higher dimensional dataset to see if there are differences in false negatives for heterogeneity detetcion for small datasets.

#### Defining heterogeneity
I create and test how three types of policy-relevant $\hat{\tau}$ distributions are detected by causal forest methods. To do so, I create the following:
1. Constant ${\tau_{ATE}}$ at 1 SD of the outcome to test if the causal forests methods are too specific and falsely identify heterogeneity when none exists.
2. 5% of observations have a $\tau_i$ with a different sign than the MDE with $\tau_i$ distributed normally around the MDE.
3. 25% of observations have a $\tau_i$ different sign than the MDE with $\tau_i$ distributed normally around the MDE.

This assumes a conservative lower bound for treatment effects that the experiment is designed to detect the minimum policy-relevant effect. This would be the worst case for detecting a type of heterogeneity.

### Testing $X$ corrleations
Since I directly generate ${\tau_i}$, I am able to precisely calculate the correlation of each $X$ with the $\tau$. I can then test what proportion of the true correlations are in the confidence interval of estimated correlations with $\hat{\tau}$ for each iteration of the simulation.

#### `grf` package set-up
The [grf package](https://grf-labs.github.io/grf/) suggests the following parameters as sensible defaults:

I use the following settings:
- Number of trees (`num_trees`): 2000
- Minimum leaf size (`min_node_size`) = 5
- honesty = TRUE
- alpha = 0.05

## Package Reference
Tibshirani, J., Athey, S., Sverdrup, E., and Wager, S. (2023). grf: Generalized Random Forests. R package version 2.4.0. https://grf-labs.github.io/grf/
Wickham H (2011). "testthat: Get Started with Testing." The R Journal, 3, 5-10. https://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf.
Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T., Miller, E., Bache, S., Müller, K., Ooms, J., Robinson, D., Seidel, D., Spinu, V., Takahashi, K., Vaughn, D., Wilke, C., Woo, K., Yutani, H. (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686. http://dx.doi.org/10.21105/joss.01686