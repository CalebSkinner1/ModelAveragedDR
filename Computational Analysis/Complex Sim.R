# Try complex relationships between confounders and Y, A
# Reproduce Simulation from Cefalu et al (2016)

library("MASS")
library("tidyverse")
library("tidymodels")
library("madr")
library("tictoc")
library("furrr")

true_effect <- as.matrix(1)
p <- 3
n <- 500
sigma <- 1

complex_sim <- function(true_effect, n, p, sigma, alpha, # propensity
                      beta) # outcome
  {
  
  # generate confounders
  x <- mvrnorm(n, mu = rep(0, p), Sigma = diag(nrow = p))
  
  # generate exposure
  a <- rbinom(n, 1, prob = expit(x[,1]*x[,2] + if_else(x[,3] > 1, 1, 0) + if_else(x[,3] > 2, 2, 0)))
  
  # generate outcome
  y <- rnorm(n, a %*% true_effect + .5*x[,1] + .5*x[,2] + if_else(x[,3] > 1, 1, 0) - x[,1]*x[,2] + if_else(x[,3] < -1, -1, 0), sigma)
  
  x_extra <- suppressMessages(as_tibble(x, .name_repair = "universal")) %>%
    rename_with(~str_replace(.x, "...", "x")) %>%
    mutate(
      x1_2 = x1^2,
      x2_2 = x2^2,
      x3_2 = x3^2,
      x1_x2 = x1*x2,
      x1_x3 = x1*x3,
      x2_x3 = x2*x3)
  
  invisible(capture.output(
    e_hat <- stepAIC(glm(factor(a) ~ ., data = tibble(a, x), family = binomial), direction = "both", k = log(n)) %>%
      fitted() %>% expit(),
    
    outcome_model <- stepAIC(lm(y ~ a, data = tibble(y, a, x)),
                             scope = list(
                               lower = as.formula(paste("y ~ a")),
                               upper = as.formula(paste("y ~ ."))),
                             direction = "both", k = log(n)),
    
    true_outcome_model <- lm(y ~ a + x[,1] + x[,2] + l1 + l2 + x[,1]:x[,2] - 1, data = tibble(y, a, x) %>%
                               mutate(l1 = if_else(x[,3] > 1, 1, 0),
                                      l2 = if_else(x[,3] < -1, 1, 0))),
    
    y1_hat <- outcome_model %>%
      predict(tibble(y, a = 1, x)),
    
    y0_hat <- outcome_model %>%
      predict(tibble(y, a = 0, x)),
    
    dr_est <- mean((a*y - (a - e_hat)*y1_hat)/e_hat - ((1-a)*y + (a - e_hat)*y0_hat)/(1-e_hat)),
    prop_est <- mean(a*y/e_hat),
    outcome_est <- mean(y1_hat - y0_hat),
    
    e_hat_extra <- stepAIC(glm(factor(a) ~ ., data = tibble(a, x_extra), family = binomial), direction = "both", k = log(n)) %>%
      fitted() %>% expit(),
    
    outcome_model_extra <- stepAIC(lm(y ~ a, data = tibble(y, a, x_extra)),
                             scope = list(
                               lower = as.formula(paste("y ~ a")),
                               upper = as.formula(paste("y ~ ."))),
                             direction = "both", k = log(n)),
    
    y1_hat_extra <- outcome_model %>%
      predict(tibble(y, a = 1, x_extra)),
    
    y0_hat_extra <- outcome_model %>%
      predict(tibble(y, a = 0, x_extra)),
    
    dr_est_extra <- mean((a*y - (a - e_hat_extra)*y1_hat_extra)/e_hat_extra - ((1-a)*y + (a - e_hat_extra)*y0_hat_extra)/(1-e_hat_extra)),
    prop_est_extra <- mean(a*y/e_hat_extra),
    outcome_est_extra <- mean(y1_hat_extra - y0_hat_extra),
    
    true_outcome_est <- true_outcome_model$coefficients[[1]],
    
    madr_est <- madr(Y = y, X = a, U = as.matrix(x), M = 1000, cut = 1)$madr[[1]],
    
    madr_est_extra <- madr(Y = y, X = a, U = as.matrix(x_extra), M = 1000, cut = 1)$madr[[1]]
  ))
  
  
  
  tibble("truth" = true_outcome_est, "prop" = prop_est, "outcome" = outcome_est,
         "DR" = dr_est, "MADR" = madr_est, "prop_EXTRA" = prop_est_extra, "outcome_EXTRA" = outcome_est_extra,
         "DR_EXTRA" = dr_est_extra, "MADR_EXTRA" = madr_est_extra)
         # , "MADR_EXTRA" = madr_est_extra)
}

# ~8 minutes
tic()
plan(multisession, workers = 4)

complex_sim_results <- future_map(1:500, ~complex_sim(true_effect, n, p, sigma),
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

