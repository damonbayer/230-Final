library(tidyverse)
library(here)
library(mvtnorm)

dat <- read.table(here("Data", "carcin.txt")) %>% 
  as_tibble() %>% 
  mutate(outcome = as_factor(outcome))

y <- dat$outcome
X <- model.matrix(outcome ~ female + treatment - 1, data = dat)

n <- nrow(X)
p <- ncol(X)

J <- nlevels(y)

n_samp <- 5000
n_burnin <- 1000

betas <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))
gammas <- matrix(NA, nrow = n_samp, ncol = J + 1, dimnames = list(NULL, str_c("gamma[", 0:J, "]")))

# Initialize with mle
mle_est <- MASS::polr(factor(outcome) ~ female + treatment, data = dat, weights = count, method = "probit")
betas[1,] <- mle_est$coefficients
gammas[1,] <- c(-Inf, mle_est$zeta, Inf)

xtx_inv <- solve(t(X) %*% X)

sample_Z <- function(X, y, beta_current, gamma_current) {
  n <- length(y)
  J <- nlevels(y)
  z <- rep(NA, n)
  
  mu <- as.vector(X %*% beta_current)
  
  for (j in 1:J) {
    current_j_tf <- y == j
    z[current_j_tf] <- rtruncnorm(n = 1,
                                  a = gamma_current[j],
                                  b = gamma_current[j+1],
                                  mean = mu[current_j_tf],
                                  sd = 1)
  }
  z
}

sample_gamma <- function(z, y, gamma_current) {
  J <- nlevels(y)
  gamma <- c(-Inf, 0, rep(NA, J - 2), Inf)
  
  for (j in 1:(J - 1)) {
    lower <- max(max(z[y == j]), gamma_current[j])
    upper <- min(min(z[y == (j + 1)]), gamma_current[j+2])
    
    if(is.na(lower) | is.na(upper)) {
      print(str_c("j: ", j))
      stop()
    }
    gamma[j+1] <- runif(1, lower, upper)
  }
  gamma
}

for (i in 2:n_samp) {
  z <- sample_Z(X, y, betas[i-1,], gammas[i-1,])
  gammas[i, ] <- sample_gamma(z, y, gammas[i-1, ])
  betas[i, ] <- rmvnorm(n = 1, as.vector(xtx_inv %*% t(X) %*% z), sigma = xtx_inv)  
}




trace_plot <- function(data) {
  data %>% 
    pivot_longer(cols = everything()) %>% 
    group_by(name) %>% 
    mutate(t = row_number()) %>% 
    ggplot(aes(x = t, y = value, color = name)) +
    facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
    geom_line() +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle("Trace Plot")
  
}

dist_plot <- function(data) {
  data %>% 
    pivot_longer(cols = everything()) %>% 
    ggplot(aes(x = value, fill = name)) +
    facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
    geom_density() +
    theme_bw() +
    theme(legend.position="none") + 
    ggtitle("Posterior Distributions")  
}

vars_post <- cbind(betas[(n_burnin+1):n_samp,],
                   gammas[(n_burnin+1):n_samp,2:J]) %>% 
  as_tibble()

colMeans(vars_post)
trace_plot(vars_post)
dist_plot(vars_post)

# TODO: Preductions