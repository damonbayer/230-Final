library(tidyverse)
library(here)
library(rstan)
rstan_options(auto_write = TRUE)


model_string <- '
  data {
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1> D;
  int<lower=1,upper=K> y[N];
  row_vector[D] x[N];
}
parameters {
  vector[D] beta;
  ordered[K-1] c;
}
model {
  vector[K] theta;
  for (n in 1:N) {
    real eta;
    eta = x[n] * beta;
    theta[1] = 1 - Phi(eta - c[1]);
    for (k in 2:(K-1))
      theta[k] = Phi(eta - c[k-1]) - Phi(eta - c[k]);
    theta[K] = Phi(eta - c[K-1]);
    y[n] ~ categorical(theta);
  }
}'



dat <- read.table(here("Data", "carcin.txt")) %>% 
  as_tibble() %>% 
  mutate(outcome = as_factor(outcome))

y <- dat$outcome
X <- model.matrix(outcome ~ female + treatment - 1, data = dat)

N <- nrow(X)
D <- ncol(X)
K <- nlevels(y)



set.seed(123456)
options(mc.cores = parallel::detectCores())  ## local multicore CPUs
dat2 = list(N=N, D = D, K=K, X=X,y=y)
model.stan = stan_model(model_code=model_string)
r= sampling(model.stan, dat2, chains = 4, iter = 4000, warmup =1000, init = "0")
write_rds(x = r, path = "~/Documents/STATS 230/230-Final/stan_samples_mp.rds")


