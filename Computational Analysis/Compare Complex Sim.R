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

complex_sim_results <- future_map(1:500, ~complex_sim(true_effect, n, p, sigma, include_BAC = TRUE),
                                  .progress = TRUE,
                                  .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
toc()

complex_sim_results <- complex_sim_results %>% rename_with(~str_replace(.x, "EXTRA", "X"))

complex_sim_results %>% replications_visual(1)

methods %>% factor(levels = c("GS", "MS-E", "MS-O", "MS-DR", "MA-DR", "BAC", "MS-E-X", "MS-O-X", "MS-DR-X", "MA-DR-X", "BAC-X"))



