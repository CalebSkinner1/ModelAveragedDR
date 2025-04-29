## ModelAveragedDR GitHub Page

This page hosts a simulation study on the Model Averaged Doubly Robust Estimator. I use the madr and BAC packages to compute the analysis.

## Abstract

The Doubly Robust (DR) estimator is consistent if either the outcome model or the propensity score model is specified correctly. If the set of confounders is unknown, then researchers may struggle to specify either model correctly. Recently, Cefalu et al. (2017) leverage advancements in causal Bayesian Model Averaging techniques to account for the uncertainty in model selection of the DR estimator. In this report, I detail the rational, theoretical results, and motivation of this Model Averaged DR (MA-DR) estimator and explore its advantages and disadvantages. In particular, I employ a three-part simulation study to explore the performance of the estimators in the presence of (1) an unknown confounding set, (2) a known confounding set with complex relationships, and (3) instrumental variables and covariates that are weakly associated with the outcome. I find the MA-DR estimator to be a powerful method for identifying the causal effect in real-world data applications.

## Computational Analysis

The Sim Function.R page contains the necessary packages to run the code. It also hosts the various wrapper functions I have written to streamline the simulations. It also includes several functions to improve analysis. The bulk of the analysis takes place in the following R Scripts.

# Compare Sim.R

This script reproduces the simulations in Cefalu et al. (2017). The main objective of these simulations is to demonstrate that the MA-DR estimator can identify the underlying set of confounding variables in the presence of a large set of confounders, and produce an unbiased estimate. I also compute the estimate from the BAC prior (Wang et al., 2012). The MA-DR and BAC estimators both produce strong results.

# Compare Complex Sim.R

This script analyzes the methods in the presence of a known set of confounders that each have complex relationships with the outcome and exposure. After adding extra terms, I find the model averaging methods are much better at approximating the true confounding relationship than the model selection methods.

# MA-DR vs BAC Sim.R

In the first two simulations, the MA-DR and BAC estimators perform similarly. I emphasize their differences with two scenarios. The MA-DR estimator performs better in the presence of instruments, and the BAC estimator performs better in the presence of covariates that are weakly associated with the outcome.
