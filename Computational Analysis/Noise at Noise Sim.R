# Simulation
library("MASS")
library("tidyverse")
library("tidymodels")
library("madr")
library("tictoc")
library("furrr")

# functions
expit <- function(logit) {
  exp(logit) / (1 + exp(logit))
}
logit <- function(expit) {
  log(expit / (1 + expit))
}


# Generate Data -----------------------------------------------------------

# unobserved values - let's say 3
samples <- 100
p_unobserved <- 3

u_rho <- 0
u <- mvrnorm(n = samples, mu = rep(0, 3), Sigma = matrix(c(1, rep(u_rho, 3), 1, rep(u_rho, 3), 1), nrow = p_unobserved))

# propensity
true_alpha <- c(1, 1, 1)

a <- rbinom(n = samples, size = 1, prob = expit(u %*% true_alpha))

true_effect <- as.matrix(1)

# outcome model
true_beta <- c(1, 1, 1)

y <- u %*% true_beta + a %*% true_effect + rnorm(samples, 0, sqrt(.1))

# generate observed values

generate_observed <- function(U, rho) {
  # U: matrix of unobserved confounders
  # rho: target correlation vector for each u_1 with generated vector
  
  U_std <- scale(U)                          # Standardize columns
  R <- cor(U_std)                            # Correlation matrix of confounders
  gamma <- solve(R, rho)                     # Solve R %*% gamma = rho
  norm_gamma2 <- sum(gamma^2)
  
  while (norm_gamma2 >= 1) {
    rho <- rho - 0.05
    U_std <- scale(U)                          # Standardize columns
    R <- cor(U_std)                            # Correlation matrix of confounders
    gamma <- solve(R, rho)                     # Solve R %*% gamma = rho
    norm_gamma2 <- sum(gamma^2)
  }
  
  z <- rnorm(nrow(U))                        # Independent noise
  x_std <- U_std %*% gamma + sqrt(1 - norm_gamma2) * z
  return(as.matrix(x_std))                  # Return as tibble
}

# observed values

# 2 observed
x <- map(1:2, ~{
  num <- .x
  generate_observed(u, rho = c(.5, .5, .5)) %>%
    as_tibble() %>% rename_with(~str_c("x", num))
  }) %>% bind_cols()
  # mutate(x3 = rnorm(samples, 0, 1)) # add third x that is just noise

# Modeling ----------------------------------------------------------------

# ideal confounders
e_hat_ideal <- logistic_reg() %>%
  fit(factor(a) ~ ., data = tibble(a, u), family = stats::binomial) %>%
  augment(tibble(a, u)) %>%
  select(.pred_1) %>% pull()

outcome_model_ideal <- linear_reg() %>%
  fit(y ~ ., data = tibble(y, a, u))

y1_hat_ideal <- outcome_model_ideal %>%
  augment(tibble(y, a = 1, u)) %>%
  select(.pred) %>% pull()

y0_hat_ideal <- outcome_model_ideal %>%
  augment(tibble(y, a = 0, u)) %>%
  select(.pred) %>% pull()

dr_est_ideal <- mean((a*y - (a - e_hat_ideal)*y1_hat_ideal)/e_hat_ideal - ((1-a)*y + (a - e_hat_ideal)*y0_hat_ideal)/(1-e_hat_ideal))
prop_est_ideal <- mean(a*y/e_hat_ideal)
outcome_est_ideal <- mean(y1_hat_ideal - y0_hat_ideal)

invisible(capture.output(
  madr_est_ideal <- madr(Y = y, X = a, U = as.matrix(u), M = 1000, cut = 1, enumerate = F, tau = NULL, two.stage = NULL)$madr[[1]]
))
c(prop_est_ideal, outcome_est_ideal, dr_est_ideal, madr_est_ideal)

# observed confounders
e_hat <- logistic_reg() %>%
  fit(factor(a) ~ ., data = tibble(a, x), family = stats::binomial) %>%
  augment(tibble(a, x)) %>%
  select(.pred_1) %>% pull()

outcome_model <- linear_reg() %>%
  fit(y ~ ., data = tibble(y, a, x))

y1_hat <- outcome_model %>%
  augment(tibble(y, a = 1, x)) %>%
  select(.pred) %>% pull()

y0_hat <- outcome_model %>%
  augment(tibble(y, a = 0, x)) %>%
  select(.pred) %>% pull()

dr_est <- mean((a*y - (a - e_hat)*y1_hat)/e_hat - ((1-a)*y + (a - e_hat)*y0_hat)/(1-e_hat))
prop_est <- mean(a*y/e_hat)
outcome_est <- mean(y1_hat - y0_hat)

invisible(capture.output(
  madr_est <- madr(Y = y, X = a, U = as.matrix(x), M = 1000, cut = 1, enumerate = F, tau = NULL, two.stage = NULL)$madr[[1]]
))
# generate extra psuedo confounders (w)

# generate 8 extra confounders
w <- map(1:8, ~{generate_observed(x, rho = c(.9, .9))}) %>% bind_cols() %>%
  rename_with(~str_replace(.x, "...", "w"))


invisible(capture.output(
  madr_extra_est <- madr(Y = y, X = a, U = as.matrix(bind_cols(x, w)), M = 1000, cut = 1, enumerate = F, tau = NULL, two.stage = NULL)$madr[[1]]
))


c("prop" = prop_est, "outcome" = outcome_est, "DR" = dr_est, "MADR" = madr_est, "MADR_extra" = madr_extra_est)


simulate_est <- function(samples = 100, p_unobserved = 3, x_observed = 2, w_extra = 8,
                         true_alpha = 1, true_beta = 1, true_effect = 1,
                         var = .1, x_rho = .5, u_rho = 0, w_rho = .9){
  
  # generate unknown confounders
  u <- mvrnorm(n = samples, mu = rep(0, p_unobserved), Sigma = matrix(c(1, rep(u_rho, p_unobserved), 1, rep(u_rho, p_unobserved), 1), nrow = p_unobserved))
  
  # generate exposures
  a <- rbinom(n = samples, size = 1, prob = expit(u %*% rep(true_alpha, p_unobserved) + rnorm(samples, 0, sqrt(var))))

  # generate outcomes
  y <- u %*% rep(true_beta, p_unobserved) + a %*% as.matrix(true_effect) + rnorm(samples, 0, sqrt(var))
  
  # observed values
  
  x <- map(1:x_observed, ~{
    num <- .x
    generate_observed(u, rho = rep(x_rho, p_unobserved)) %>% as_tibble() %>% rename_with(~str_c("x", num))}) %>% bind_cols()
  # mutate(x3 = rnorm(samples, 0, 1)) # add third x that is just noise
  
  # Modeling ----------------------------------------------------------------
  
  # ideal confounders
  e_hat_ideal <- logistic_reg() %>%
    fit(factor(a) ~ ., data = tibble(a, u), family = stats::binomial) %>%
    augment(tibble(a, u)) %>%
    select(.pred_1) %>% pull()
  
  outcome_model_ideal <- linear_reg() %>%
    fit(y ~ ., data = tibble(y, a, u))
  
  y1_hat_ideal <- outcome_model_ideal %>%
    augment(tibble(y, a = 1, u)) %>%
    select(.pred) %>% pull()
  
  y0_hat_ideal <- outcome_model_ideal %>%
    augment(tibble(y, a = 0, u)) %>%
    select(.pred) %>% pull()
  
  dr_est_ideal <- mean((a*y - (a - e_hat_ideal)*y1_hat_ideal)/e_hat_ideal - ((1-a)*y + (a - e_hat_ideal)*y0_hat_ideal)/(1-e_hat_ideal))
  prop_est_ideal <- mean(a*y/e_hat_ideal)
  outcome_est_ideal <- mean(y1_hat_ideal - y0_hat_ideal)
  
  invisible(capture.output(
    madr_est_ideal <- madr(Y = y, X = a, U = as.matrix(u), M = 1000, cut = 1, enumerate = F, tau = NULL, two.stage = NULL)$madr[[1]]
  ))
  
  
  # observed confounders
  e_hat <- logistic_reg() %>%
    fit(factor(a) ~ ., data = tibble(a, x), family = stats::binomial) %>%
    augment(tibble(a, x)) %>%
    select(.pred_1) %>% pull()
  
  outcome_model <- linear_reg() %>%
    fit(y ~ ., data = tibble(y, a, x))
  
  y1_hat <- outcome_model %>%
    augment(tibble(y, a = 1, x)) %>%
    select(.pred) %>% pull()
  
  y0_hat <- outcome_model %>%
    augment(tibble(y, a = 0, x)) %>%
    select(.pred) %>% pull()
  
  dr_est <- mean((a*y - (a - e_hat)*y1_hat)/e_hat - ((1-a)*y + (a - e_hat)*y0_hat)/(1-e_hat))
  prop_est <- mean(a*y/e_hat)
  outcome_est <- mean(y1_hat - y0_hat)
  
  invisible(capture.output(
    madr_est <- madr(Y = y, X = a, U = as.matrix(x), M = 1000, cut = .95, enumerate = F, tau = NULL, two.stage = NULL)$madr[[1]]
  ))
  
  # generate extra psuedo confounders (w)
  
  # generate 8 extra confounders
  w <- map(1:w_extra, ~{
    num <- .x
    generate_observed(x, rho = rep(w_rho, x_observed)) %>% as_tibble() %>% rename_with(~str_c("w", num))}) %>% bind_cols()
  
  invisible(capture.output(
    madr_extra_est <- madr(Y = y, X = a, U = as.matrix(bind_cols(x, w)), M = 1000, cut = .95, enumerate = F, tau = NULL, two.stage = NULL)$madr[[1]]
  ))
  
  tibble("prop_ideal" = prop_est_ideal, "outcome_ideal" = outcome_est_ideal, "dr_ideal" = dr_est_ideal, "madr_ideal" = madr_est_ideal,
         "prop" = prop_est, "outcome" = outcome_est, "DR" = dr_est, "MADR" = madr_est, "MADR_extra" = madr_extra_est)
}

sim_w8 <- map_dfr(1:100, ~simulate_est(w_rho = .8))

sim_w7 <- map_dfr(1:100, ~simulate_est(w_rho = .7))

sim_w6 <- map_dfr(1:100, ~simulate_est(w_rho = .7))

sim_w6 %>% group_by() %>% summarize(
  across(everything(), ~mean(.x))
)

