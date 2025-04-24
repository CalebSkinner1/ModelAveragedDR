# Functions for simulations

# load libraries
library("MASS")
library("BAC")
library("tidyverse"); theme_set(theme_minimal())
library("tidymodels")
library("madr")
library("tictoc")
library("furrr")
library("mvnfast")

true_effect <- as.matrix(1)

# Functions that produce simulations --------------------------------------

# this function reproduces the scenarios in the paper
paper_sim <- function(true_effect, n, p, sigma, alpha, # propensity
                      beta, # outcome
                      include_BAC = FALSE)
{
  # generate confounders
  x <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(nrow = p))
  
  # generate exposure
  a <- rbinom(n, 1, prob = expit(x%*%alpha))
  
  # generate outcome
  y <- rnorm(n, a %*% true_effect + x %*% beta, sigma)
  
  invisible(capture.output(
    e_hat <- MASS::stepAIC(glm(factor(a) ~ ., data = tibble(a, x), family = binomial), direction = "both", k = log(n)) %>%
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
    
    madr_est <- madr(Y = y, X = a, U = as.matrix(x), M = 1000, cut = 1)$madr[[1]],
    
    if(include_BAC){
      bac_chains <- 3
      bac_coefs <- BAC(Y = y, X = a, D = x, Nsims = 1000, chains = bac_chains)$coefs
      
      bac_est <- map(1:bac_chains, ~{
        bac_coefs["Outcome", .x, , "X"] %>% as.vector()
      }) %>% list_c() %>% mean()
      
      return(tibble("truth" = true_outcome_est, "prop" = prop_est, "outcome" = outcome_est,
                    "DR" = dr_est, "MADR" = madr_est, "BAC" = bac_est))
    }
  ))
  
  tibble("truth" = true_outcome_est, "prop" = prop_est, "outcome" = outcome_est,
         "DR" = dr_est, "MADR" = madr_est)
}

# this function simulates results under complex relationships between confounders and exposure and outcome
complex_sim <- function(true_effect, n, p, sigma, alpha, # propensity
                        beta,  # outcome
                        include_BAC = FALSE)
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
    
    madr_est_extra <- madr(Y = y, X = a, U = as.matrix(x_extra), M = 1000, cut = 1)$madr[[1]],
    
    if(include_BAC){
      bac_chains <- 3
      bac_coefs <- BAC(Y = y, X = a, D = x, Nsims = 1000, chains = bac_chains)$coefs
      
      bac_est <- map(1:bac_chains, ~{
        bac_coefs["Outcome", .x, , "X"] %>% as.vector()
      }) %>% list_c() %>% mean()
      
      bac_coefs_extra <- BAC(Y = y, X = a, D = x_extra, Nsims = 1000, chains = bac_chains)$coefs
      
      bac_est_extra <- map(1:bac_chains, ~{
        bac_coefs_extra["Outcome", .x, , "X"] %>% as.vector()
      }) %>% list_c() %>% mean()
      
      return(tibble("truth" = true_outcome_est, "prop" = prop_est, "outcome" = outcome_est,
                    "DR" = dr_est, "MADR" = madr_est, "BAC" = bac_est,
                    "prop_EXTRA" = prop_est_extra, "outcome_EXTRA" = outcome_est_extra,
                    "DR_EXTRA" = dr_est_extra, "MADR_EXTRA" = madr_est_extra, "BAC_EXTRA" = bac_est_extra))
    }
  ))
  
  tibble("truth" = true_outcome_est, "prop" = prop_est, "outcome" = outcome_est,
         "DR" = dr_est, "MADR" = madr_est,
         "prop_EXTRA" = prop_est_extra, "outcome_EXTRA" = outcome_est_extra,
         "DR_EXTRA" = dr_est_extra, "MADR_EXTRA" = madr_est_extra)
}

dr_bac_sim <- function(true_effect, n, p, sigma, alpha, # propensity
                       beta) # outcome
{
  
  # generate confounders
  x <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(nrow = p))
  
  # generate exposure
  a <- rbinom(n, 1, prob = expit(x%*%alpha))
  
  # generate outcome
  y <- rnorm(n, a %*% true_effect + x %*% beta, sigma)
  
  invisible(capture.output(
    madr_est <- madr(Y = y, X = a, U = as.matrix(x), M = 1000, cut = 1)$madr[[1]],
    
    bac_chains <- 3, 
    bac_coefs <- BAC(Y = y, X = a, D = x, Nsims = 1000, chains = bac_chains)$coefs,
    
    bac_est <- map(1:bac_chains, ~{
      bac_coefs["Outcome", .x, , "X"] %>% as.vector()
    }) %>% list_c() %>% mean()
  ))
  
  tibble("MADR" = madr_est, "BAC" = bac_est)
}


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

replications_visual <- function(df, true_effect){
  methods <- colnames(dr_bac_iv)
  
  dr_bac_iv %>%
    summarize_all(~list(mean = mean(.x), lower = quantile(.x, .025), upper = quantile(.x, .975)) %>% unlist()) %>%
    mutate(quantity = c("mean", "lower", "upper")) %>%
    pivot_wider(values_from = all_of(methods), names_from = quantity) %>%
    pivot_longer(cols = everything(), names_to = c("method", "quantity"), names_sep = "_", values_to = "value") %>%
    pivot_wider(names_from = quantity, values_from = value) %>%
    ggplot() +
    geom_point(aes(y = mean, x = method, color = method)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper, x = method, color = method), width = .2) +
    theme(legend.title = element_blank()) +
    geom_hline(yintercept = true_effect, linetype = "dotted") +
    labs(x = "Method", y = "Estimated Effect")
}



