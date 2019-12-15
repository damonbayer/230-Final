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
    y[n] ~ bernoulli_logit(dot_product(X[n],beta));
  }
}
'


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

N <-dim(X)[1]
K <- dim(X)[2]


set.seed(123456)
options(mc.cores = parallel::detectCores())  ## local multicore CPUs
dat = list(N=N, K=K, X=X,y=y)
model.stan = stan_model(model_code=model_string)
r = sampling(model.stan, dat, chains = 1, iter = 6000, warmup =2000, init = "0")
write_rds(x = r, path = "~/Documents/STATS 230/230-Final/stan_samples_logit.rds")

readRDS("~/Documents/STATS 230/230-Final/stan_samples_logit.rds")
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
