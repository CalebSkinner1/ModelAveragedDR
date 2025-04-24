# Clearly, model averaging is valuable when the true confounding set is unknown.

source("ModelAveragedDR/Computational Analysis/Sim Functions.R")

# scenario 1 - instrumental variables
alpha_iv <- c(1, 2, 3, .5, .5, .5, 0, 0, 0, 0)
beta_iv <- c(0, 0, 0, .5, 1, 1.5, 0, 0, 0, 0)

p_iv <- length(alpha_iv)
n_iv <- 200
sigma_iv <- 2

# ~12 minutes
plan(multisession, workers = 4)
tic()
dr_bac_iv <- future_map(1:500, ~dr_bac_sim(true_effect, n_iv, p_iv, sigma_iv, alpha_iv, beta_iv),
                                  .progress = TRUE,
                                  .options = furrr_options(seed = TRUE)) %>%
  bind_rows()

toc()

dr_bac_iv %>% replications_l1_bias(1)
dr_bac_iv %>% replications_mean()
dr_bac_iv %>% replications_visual(1)

# scenario 2 - strong association with exposure, and weak association with outcome

alpha_e <- c(.5, 1, 1.5, 2, 2.5, 0, 0, 0, 0, 0)
beta_e <- c(.1, .2, .2, .2, .1, 0, 0, 0, 0, 0)

p_e <- length(alpha_e)
n_e <- 200
sigma_e <- 2

# ~12 minutes
plan(multisession, workers = 4)
tic()
dr_bac_e <- future_map(1:500, ~dr_bac_sim(true_effect, n_e, p_e, sigma_e, alpha_e, beta_e),
                        .progress = TRUE,
                        .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
toc()

dr_bac_e %>% replications_l1_bias(1)
dr_bac_e %>% replications_mean()
dr_bac_e %>% replications_visual(1)


