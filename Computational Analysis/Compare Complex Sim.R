# try complex simulations with complex relationships between Y and A
# these simulations are essentially identical to that in Complex Sim.R, except add the BAC prior
# result as a competitor

source("ModelAveragedDR/Computational Analysis/Sim Functions.R")

p <- 3
n <- 500
sigma <- 1

# ~17 minutes
tic()
plan(multisession, workers = 4)

complex_sim_results <- future_map(1:500, ~complex_compare_sim(true_effect, n, p, sigma, include_BAC = TRUE),
                                  .progress = TRUE,
                                  .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
toc()