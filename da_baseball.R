library(tidyverse)
library(here)
library(mvtnorm)
library(truncnorm)

# Setup -------------------------------------------------------------------
dat <- read_csv(here("Data", "LAD_assessment_data.csv"),
                col_types = cols(batter_mlb_id = col_skip(), 
                                 pitcher_mlb_id = col_skip(), 
                                 venue_city = col_skip())) %>% drop_na()


X <- model.matrix(is_in_play ~ ., data = dat)
y <- dat$is_in_play

n <- nrow(X)
p <- ncol(X)

n_samp <- 5000
n_burnin <- 1000

# Function Definitions ----------------------------------------------------
generate_Z <- function(y, X, beta) {
  n <- length(y)
  mu <- X %*% beta
  
  z <- rep(NA, n)
  
  z[y == 0] <- rtruncnorm(n = 1,
                          b = 0,
                          mean = mu[y == 0],
                          sd=1)
  
  z[y == 1] <- rtruncnorm(n = 1,
                          a = 0,
                          mean = mu[y == 1],
                          sd = 1)
  
  z
}

# Initialization ----------------------------------------------------------
beta <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))
set.seed(230)
beta[1,] <- rnorm(p)
dat_z <- dat %>%
  select(-is_in_play) %>%
  mutate(z = generate_Z(y, X, beta[1,]))

sigma <- solve(t(X) %*% X)


# Gibbs Sampling ---------------------------------------------------------
set.seed(230)
timer <- rbenchmark::benchmark(
for (i in 2:n_samp){
  beta_z <- lm(z ~ ., data = dat_z)$coeff
  beta[i,] <- rmvnorm(n = 1, mean = beta_z, sigma = sigma)
  dat_z$z <- generate_Z(y, X, beta[i,])
}, replications = 1)



# Post processing ---------------------------------------------------------
# Throw away burn-in
beta_post <- beta[(n_burnin+1):n_samp,] %>% 
  as_tibble()

colMeans(beta_post)

probit_model <- glm(is_in_play ~ ., family = binomial(link = "probit"), data = dat)
summary(probit_model)

write_csv(beta_post,"~/Documents/STATS 230/230-Final/baseball_betas.rds")
# Analysis ----------------------------------------------------------------

source(here("plotting_functions.R"))
trace_plot(beta_post)
dist_plot(beta_post)
