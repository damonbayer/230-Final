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
  select(-'Sample code number') %>% 
  select(-"Bare Nuclei") %>% 
  mutate(Class = Class / 2 - 1)

X <- model.matrix(Class ~ ., data = dat)
y <- dat$Class

n <- nrow(X)
p <- ncol(X)

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

n_samp <- 4000
n_burnin <- 1000

beta <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))
set.seed(230)
beta[1,] <- rnorm(p)

dat_z <- dat %>%
  select(-Class) %>%
  mutate(z = generate_Z(y, X, beta[1,]))

sigma <- solve(t(X) %*% X)

timer <- rbenchmark::benchmark(
for (i in 2:n_samp){
  beta_z <- lm(z ~ ., data = dat_z)$coeff
  beta[i,] <- rmvnorm(n = 1, mean = beta_z, sigma = sigma)
  dat_z$z <- generate_Z(y, X, beta[i,])
}, replications = 1)

# Throw away burn-in
beta_post <- beta[(n_burnin+1):n_samp,] %>% 
  as_tibble()

write_rds(beta_post, path = here("gibbs_samples_bc.rds"))

colMeans(beta_post)

probit_model <- glm(Class ~ ., family = binomial(link = "probit"), data = dat)
logit_model <- glm(Class ~ ., family = "binomial", data = dat)

summary(probit_model)

beta_post %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = value, fill = name)) +
  facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
  geom_density() +
  theme_bw() +
  theme(legend.position="none") + 
  ggtitle("Posterior Distributions")

beta_post %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  mutate(t = row_number()) %>% 
  ggplot(aes(x = t, y = value, color = name)) +
  facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
  geom_line() +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Trace Plot")
