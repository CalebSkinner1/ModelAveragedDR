# Reproduce Simulation from Cefalu et al (2016)

source("ModelAveragedDR/Computational Analysis/Sim Functions.R")

# scenario 1 n = 200, p = 5, sigma^2 = 4 -------------------------------------------------------------------------
p1 <- 5
n <- 200
sigma1 <- 2

# four scenarios for PS and outcome models
# no confounding
alpha11 <- c(.4, .3, .2, .1, 0)
beta11 <- c(0, 0, 0, 0, 0)

# moderate confounding
alpha12 <- c(.5, .5, .1, 0, 0)
beta12 <- c(.5, 0, 1, .5, 0)

# strong predictors, weak treatment assoc
alpha13 <- c(.1, .1, 1, 1, 1)
beta13 <- c(2, 2, 0, 0, 0)

# strong confounding
alpha14 <- c(.5, .4, .3, .2, .1)
beta14 <- c(.5, 1, 1.5, 2, 2.5)

# 3-4 minutes
tic()
plan(multisession, workers = 4)

situation11 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma1, alpha11, beta11, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation12 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma, alpha12, beta12, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation13 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma, alpha13, beta13, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation14 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma, alpha14, beta14, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
toc()

# scenario 2 n = 200, p = 15, sigma^2 = 1 ------------------------------------------
# p = 100 is not feasible. I'm not sure how they did this in the paper. Running it once took 2 hours.

p2 <- 15
sigma2 <- 1

alpha21 <- c(.4, .3, .2, .1, 0, rep(0, p2-5))
beta21 <- c(0, 0, 0, 0, 0, rep(0, p2-5))

# moderate confounding
alpha22 <- c(.5, .5, .1, 0, 0, rep(0, p2-5))
beta22 <- c(.5, 0, 1, .5, 0, rep(0, p2-5))

# strong predictors, weak treatment assoc
alpha23 <- c(.1, .1, 1, 1, 1, rep(0, p2-5))
beta23 <- c(2, 2, 0, 0, 0, rep(0, p2-5))

# strong confounding
alpha24 <- c(.5, .4, .3, .2, .1, rep(0, p2-5))
beta24 <- c(.5, 1, 1.5, 2, 2.5, rep(0, p2-5))

plan(multisession, workers = 4)
tic()

situation21 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha21, beta21, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()

situation22 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha22, beta22, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()

situation23 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha23, beta23, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation24 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha24, beta24, include_BAC = FALSE),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
toc()




