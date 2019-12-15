library(tidyverse)
library(here)
library(truncnorm)
library(mvtnorm)
library(crch)

dat <- read_csv(here("Data", "breast-cancer-wisconsin.data"),
                col_names = c('Sample code number',
                              'Clump Thickness',
                              'Uniformity of Cell Size',
                              'Uniformity of Cell Shape',
                              'Marginal Adhesion',
                              'Single Epithelial Cell Size',
                              'Bare Nuclei',
                              'Bland Chromatin',
                              'Normal Nucleoli',
                              'Mitoses',
                              'Class'), na = '?') %>% 
  select(-'Sample code number') %>% 
  select(-"Bare Nuclei") %>% 
  mutate(Class = Class / 2 - 1)

X <- model.matrix(Class ~ ., data = dat)
y <- dat$Class

n <- nrow(X)
p <- ncol(X)


# Sample from a truncated t distribution -----------------------------
generate_Z <- function(y, X, beta,nu) {
  n <- length(y)
  mu <- X %*% beta
  
  z <- rep(NA, n)
  
  z[y == 0] <- rtt(n = 1,
                   right = 0,
                   location = mu[y == 0],
                   scale = 1,
                   df = nu)
  
  z[y == 1] <- rtt(n = 1,
                   left = 0,
                   location = mu[y == 1],
                   scale = 1,
                   df = nu)
  
  z
}


# Initailize Model ----------------------------
n_samp <- 6000
n_burnin <- 2000

nu <- 8

beta <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))

set.seed(230)
beta[1,] <- glm(Class ~ ., family = binomial(link = "probit"), data = dat)$coeff #start with probit

lambda <- rep(1,n)

dat_z <- dat %>%
  select(-Class) %>%
  mutate(z = generate_Z(y, X, beta[1,],nu))



timer <- rbenchmark::benchmark(
  for (i in 2:n_samp){
    W <- diag(lambda)
    sigma <- solve(t(X) %*% W %*% X)
    beta_z <- sigma %*% t(X) %*% W %*% dat_z$z
    beta[i,] <- rmvnorm(n=1, beta_z, sigma)

    lambda <- rgamma(n, shape = (nu+1)/2, scale = (nu + (dat_z$z - X %*% beta[i,])^2)/2)
    dat_z$z <- generate_Z(y, X, beta[i,],nu)
  }, replications = 1)

# Throw away burn-in
beta_post <- beta[(n_burnin+1):n_samp,] %>% 
  as_tibble()

# write_rds(beta_post, path = here("gibbs_samples_bc.rds"))

colMeans(beta_post/.634)

probit_model <- glm(Class ~ ., family = binomial(link = "probit"), data = dat)
logit_model <- glm(Class ~ ., family = "binomial", data = dat)

summary(logit_model)
source("plotting_functions.R")
trace_plot(beta_post/.634)
dist_plot(beta_post/.634)


write_rds(x = beta_post, path = "~/Documents/STATS 230/230-Final/Gibbs_logit_bc.rds")
