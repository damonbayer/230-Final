library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)


#Setup model string ------------------------------------------------
model_string <- '
data {
  int N; // number of obs 
  int K; // number of predictors
  
  int y[N]; // outcome
  row_vector[K] X[N]; // predictors
}
parameters {
  vector[K] beta;
}
model {
  beta ~ normal(0,100);
  for(n in 1:N) {
    y[n] ~ bernoulli(Phi_approx(dot_product(X[n],beta)));
  }
}
'


# Data Setup ------------------------------------------------
dat <- read_csv("~/Desktop/LAD_assessment_data.csv",
                col_types = cols(batter_mlb_id = col_skip(), 
                                 pitcher_mlb_id = col_skip(), 
                                 venue_city = col_skip())) %>% drop_na()


X <- model.matrix(is_in_play ~ ., data = dat)
y <- dat$is_in_play

N <-dim(X)[1]
K <- dim(X)[2]


#STAN setup and sampling ------------------------------------------------
set.seed(123456)
options(mc.cores = parallel::detectCores())  ## local multicore CPUs
dat_model = list(N=N, K=K, X=X,y=y)
model.stan = stan_model(model_code=model_string)
r_baseball = sampling(model.stan, dat_model, chains = 1, iter = 5000, warmup =1000, init = "0")
write_rds(x = r, path = "~/Documents/STATS 230/230-Final/stan_samples_baseball.rds")

#load("~/Documents/STATS 230/230-Final/stan_samples_baseball.rds")
summary(r_baseball)$summary
samples <- as.matrix(r_baseball)[,-10]


#Analysis ------------------------------------------------
source(here("plotting_functions.R"))
trace_plot(samples)
dist_plot(samples)
