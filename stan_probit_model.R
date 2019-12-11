library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)

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


data <- read_csv("~/Documents/STATS 230/230-Final/Data/breast-cancer-wisconsin.data", 
                col_names = c("ID", "thick","uni_size",
                              "uni_shape","adhesion",
                              "single_size","nuclei",
                              "chrom","nucleoli","mitoses",
                              "class"))
data$nuclei <- as.numeric(data$nuclei)
data$class <- data$class/2 - 1

y <- data$class
X <- model.matrix(~thick+uni_size+uni_shape+adhesion+single_size+chrom+nucleoli+mitoses,data=data)

N <-dim(X)[1]
K <- dim(X)[2]


set.seed(123456)
options(mc.cores = parallel::detectCores())  ## local multicore CPUs
dat = list(N=N, K=K, X=X,y=y)
model.stan = stan_model(model_code=model_string)
r = sampling(model.stan, dat, chains = 6, iter = 10000, warmup =2000, thin = 100)
write_rds(x = r, path = "~/Documents/STATS 230/230-Final/stan_samples.rds")
summary(r)


library(brms)
fit <- brm(class ~ thick+uni_size+uni_shape+adhesion+single_size+chrom+nucleoli+mitoses, 
           data = data, family = binomial(link = "probit"))
summary(fit)

