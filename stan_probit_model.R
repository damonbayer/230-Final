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
r = sampling(model.stan, dat, chains = 1, iter = 5000, warmup =1000, init = "0")
write_rds(x = r, path = "~/Documents/STATS 230/230-Final/stan_samples.rds")

#load("~/Documents/STATS 230/230-Final/stan_samples.rds")
summary(r)$summary
samples <- as.matrix(r)[,-10]



samples %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = value, fill = name)) +
  facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
  geom_density() +
  theme_bw() +
  theme(legend.position="none") + 
  ggtitle("Posterior Distributions")

samples %>% 
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

# library(brms)
# fit <- brm(class ~ thick+uni_size+uni_shape+adhesion+single_size+chrom+nucleoli+mitoses, 
#            data = data, family = binomial(link = "probit"))
# summary(fit)
# 
