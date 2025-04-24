# Try complex relationships between confounders and Y, A
# Reproduce Simulation from Cefalu et al (2016)

source("ModelAveragedDR/Computational Analysis/Sim Functions.R")

p <- 3
n <- 500
sigma <- 1

# ~8 minutes
tic()
plan(multisession, workers = 4)

complex_sim_results <- future_map(1:500, ~complex_sim(true_effect, n, p, sigma, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
toc()

complex_sim_results %>%
  mutate(across(everything(), ~.x - 1)) %>%
  group_by() %>%
  summarize(
    across(.cols = everything(), ~mean(abs(.x))) # mean absolute error
  )

complex_sim_results$truth %>% mean()

