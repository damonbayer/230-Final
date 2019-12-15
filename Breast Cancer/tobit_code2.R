#tobit model 

library(tidyverse)
library(here)
library(truncnorm)
library(mvtnorm)

log_sum_exp <- function(values) {
  max_val <- max(values)
  max_val + log(sum(exp(values - max_val)))
}
dat <- read_csv(here("Breast Cancer", "breast-cancer-wisconsin.data"),
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
  select(-'Sample code number', -'Bare Nuclei') %>% 
  mutate(Class = Class / 2 - 1)

# dat <- sample_n(dat, 20)

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

c_func <- function(nu)   1 / (gamma(nu / 2) * (nu / 2)^(nu / 2))

n_samp <- 6000
n_burnin <- 2000

# nu_choices <- c(4, 8, 16, 32)
nu_choices <- 10^(2:5)
nu_prior <- rep(1 / length(nu_choices), length(nu_choices))

beta <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))
nu <- rep(NA, n_samp)


beta[1,] <- glm(Class ~ ., family = binomial(link = "probit"), data = dat)$coeff
nu[1] <- nu_choices[sample(x = length(nu_choices), size = 1, prob = nu_prior)]
# nu[1] <- nu_choices[1]

augmented_data <- dat %>%
  select(-Class) %>%
  mutate(lambda = 1)

set.seed(230)

rbenchmark::benchmark(
for (i in 2:n_samp){
  print(i)
  #generate z
  augmented_data$z <- generate_Z(y, X, beta[i-1,],augmented_data$lambda)
  
  W <- diag(augmented_data$lambda)
  
  # calculate sigma
  sigma <- solve(t(X) %*% W %*% X)
  
  #weight least squares
  beta_z <- as.vector(lm(z ~ . -lambda, data = augmented_data, weights = lambda)$coeff)
  
  # new betas
  beta[i,] <- rmvnorm(n = 1, mean = beta_z, sigma = sigma)
  
  # new lambdas
  shape <- (nu[i - 1] + 1) / 2
  
  scales <-  as.vector(2 / (nu[i - 1] + (augmented_data$z - X %*% beta[i,])^2))
  
  augmented_data$lambda <- vapply(scales, function(scale) rgamma(n = 1, shape = shape, scale = scale), FUN.VALUE = 1.0)
  
  log_pos_unnormal <- sapply(nu_choices, function(nu) log(nu_prior[nu_choices == nu]) + n * log(c_func(nu)) + sum((nu / 2 - 1) * log(augmented_data$lambda) - nu * augmented_data$lambda / 2))
  
  post_normalized <- exp(log_pos_unnormal - log_sum_exp(log_pos_unnormal))
  
  nu[i] <- nu_choices[sample(x = length(nu_choices), size = 1, prob = post_normalized)]
  # nu[i] <- nu[i-1]
}, replications = 0)

beta_post <- beta[(n_burnin+1):n_samp,] %>% 
  as_tibble()


write_rds(beta_post, here("Breast Cancer", "gibbs_logit_bc.rds"))

beta_post %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = value, fill = name)) +
  facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
  geom_density() +
  theme_bw() +
  theme(legend.position="none") + 
  ggtitle("Posterior Distributions")

beta_post %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  mutate(t = row_number()) %>% 
  ggplot(aes(x = t, y = value, color = name)) +
  facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
  geom_line() +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Trace Plot")


logit_model <- glm(Class ~ ., family = "binomial", data = dat)
glm(Class ~ ., family = "binomial", data = dat)
summary(logit_model)
colMeans(beta_post)
apply(beta_post, MARGIN = 2, median)

probit_model <- glm(Class ~ ., family = binomial(link = "probit"), data = dat)
probit_model$coefficients
colMeans(beta_post)
apply(beta_post, MARGIN = 2, median)

table(nu)
