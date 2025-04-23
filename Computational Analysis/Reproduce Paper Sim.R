# Reproduce Simulation from Cefalu et al (2016)

library("MASS")
library("tidyverse")
library("tidymodels")
library("madr")
library("tictoc")
library("furrr")

true_effect <- as.matrix(1)

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

paper_sim <- function(true_effect, n, p, sigma, alpha, # propensity
                      beta, # outcome
                      approximate)
{
  
  # generate confounders
  x <- mvrnorm(n, mu = rep(0, p), Sigma = diag(nrow = p))
  
  # generate exposure
  a <- rbinom(n, 1, prob = expit(x%*%alpha))
  
  # generate outcome
  y <- rnorm(n, a %*% true_effect + x %*% beta, sigma)
  
  invisible(capture.output(
    e_hat <- stepAIC(glm(factor(a) ~ ., data = tibble(a, x), family = binomial), direction = "both", k = log(n)) %>%
      fitted() %>% expit(),
    
    outcome_model <- stepAIC(lm(y ~ a, data = tibble(y, a, x)),
                             scope = list(
                               lower = as.formula(paste("y ~ a")),
                               upper = as.formula(paste("y ~ ."))),
                             direction = "both", k = log(n)),
    
    true_outcome_model <- suppressWarnings(lm(y ~ . - 1, data = tibble(y, a, as_tibble(x[,which(beta != 0)])))),
    
    y1_hat <- outcome_model %>%
      predict(tibble(y, a = 1, x)),
    
    y0_hat <- outcome_model %>%
      predict(tibble(y, a = 0, x)),
    
    dr_est <- mean((a*y - (a - e_hat)*y1_hat)/e_hat - ((1-a)*y + (a - e_hat)*y0_hat)/(1-e_hat)),
    prop_est <- mean(a*y/e_hat),
    outcome_est <- mean(y1_hat - y0_hat),
    true_outcome_est <- true_outcome_model$coefficients[[1]],
    
    madr_est <- madr(Y = y, X = a, U = as.matrix(x), M = 1000, cut = 1)$madr[[1]]
  ))
  
  tibble("truth" = true_outcome_est, "prop" = prop_est, "outcome" = outcome_est, "DR" = dr_est, "MADR" = madr_est)
}

# 3-4 minutes
tic()
plan(multisession, workers = 4)

situation11 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma1, alpha11, beta11),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation12 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma, alpha12, beta12),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation13 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma, alpha13, beta13),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation14 <- future_map(1:500, ~paper_sim(true_effect, n, p1, sigma, alpha14, beta14),
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

situation21 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha21, beta21),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()

situation22 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha22, beta22),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()

situation23 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha23, beta23),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
situation24 <- future_map(1:500, ~paper_sim(true_effect, n, p2, sigma2, alpha24, beta24),
                          .progress = TRUE,
                          .options = furrr_options(seed = TRUE)) %>%
  bind_rows()
toc()


# Functions for quick analysis of replications ----------------------------

replications_mean <- function(df){
  df %>% group_by() %>%
    summarize(across(everything(), ~mean(.x)))
}

replications_l1_bias <- function(df, true_effect){
  df %>%
    mutate(across(everything(), ~abs(.x - true_effect))) %>%
    group_by() %>%
    summarize(across(everything(), ~mean(.x)))
}



