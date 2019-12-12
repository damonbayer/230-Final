#tobit model 

library(tidyverse)
library(here)
library(truncnorm)
library(mvtnorm)

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
  na.omit() %>% 
  select(-'Sample code number') %>% 
  mutate(Class = Class / 2 - 1)

X <- model.matrix(Class ~ ., data = dat)
y <- dat$Class

n <- nrow(X)
p <- ncol(X)


generate_Z <- function(y, X, beta,lambda) {
  n <- length(y)
  mu <- X %*% beta
  
  z <- rep(NA, n)
  
  z[y == 0] <- rtruncnorm(n = 1,
                          b = 0,
                          mean = mu[y == 0],
                          sd= 1/lambda[y == 0])
  
  z[y == 1] <- rtruncnorm(n = 1,
                          a = 0,
                          mean = mu[y == 1],
                          sd = 1/lambda[y == 1])
  
  z
}

n_samp <- 5000
n_burnin <- 1000


beta <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))

beta[1,] <- glm(Class ~ ., family = binomial(link = "probit"), data = dat)$coeff

augmented_data <- dat %>%
  select(-Class) %>%
  mutate(lambda = 1)

set.seed(230)
nu_choices <- c(4,8,16,32)
nu <- c(sample(nu_choices,1))

for (i in 2:n_samp){
  #generate z
  augmented_data$z <- generate_Z(y, X, beta[i-1,],augmented_data$lambda)
  # calculate sigma
  sigma <- solve(t(X) %*% diag(augmented_data$lambda) %*% X)
  #weight least squares
  beta_z <- lm(z ~ .-lambda, data = augmented_data, weights = lambda)$coeff
  #new betas
  beta[i,] <- rmvnorm(n = 1, mean = beta_z, sigma = sigma)
  # new lambda
  shape <- (nu[i-1] + 1)/2
  rate <- c(2/(nu[i-1] + (augmented_data$z - X %*% beta[i,])^2))
  augmented_data$lambda <- rgamma(n, shape = shape, rate = rate)
  nu[i] <- generate_nu(augmented_data$lambda, nu_choices)
}


c_function <- function(nu){
  1/gamma(nu/2)*((nu/2)^(nu/2))
}


generate_nu <- function(lambda, nu_choices,n){
  prior <- log(rep(1/length(nu_choices),length(nu_choices)))
  lik <- n*log(c_function(nu_choices)) + ((nu_choices/2)-1)*sum(log(lambda)) - (nu*sum(lambda))/2
  post <- prior + lik
  exp(post)
  ##still need to sample from post
}


