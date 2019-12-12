library(tidyverse)
library(here)
library(mvnfast)
# dat <- read.table( "https://www.ics.uci.edu/~dgillen/STAT211/Data/travel.txt" ) %>% 
#   as_tibble() %>%
#   mutate(travel = rep(0:3, 210)) %>% 
#   filter(mode==1) %>% 
#   select(travel, hinc, psize)

dat <- read.table(here("Data", "carcin.txt")) %>% 
  as_tibble() %>% 
  mutate(outcome = as_factor(outcome))

y <- dat$outcome
X <- model.matrix(outcome ~ female + treatment, data = dat)

n <- nrow(X)
p <- ncol(X)

j_max <- nlevels(y)
  
n_samp <- 5000
n_burnin <- 1000

betas <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))
gammas <- matrix(NA, nrow = n_samp, ncol = j_max + 1, dimnames = list(NULL, str_c("gamma[", 0:j_max, "]")))

# Initialize with mle
mle_est <- MASS::polr(factor(outcome) ~ female + treatment, data = dat, weights = count, method = "probit")
betas[1,] <- c(mle_est$zeta[1], mle_est$coefficients)
gammas[1,] <- c(-Inf, mle_est$zeta - mle_est$zeta[1], Inf)

xtx_inv <- solve(t(X) %*% X)

sample_Z <- function(X, y, beta_current, gamma_current) {
  n <- length(y)
  j_max <- max(as.numeric(y))
  z <- rep(NA, n)
  
  mu <- as.vector(X %*% beta_current)
  for (j in 1:j_max) {
    current_j_tf <- y == j
    z[current_j_tf] <- rtruncnorm(n = 1,
                                  a = gamma_current[j],
                                  b = gamma_current[j+1],
                                  mean = mu[current_j_tf])
  }
  z
}

sample_gamma <- function(z, y, gamma_current) {
  j_max <- max(as.numeric(y))
  gamma <- c(-Inf, 0, rep(NA, j_max - 2), Inf)
  
  for (j in 2:(j_max - 1)) {
    
    if (is.na(max(max(z[y == j]), gamma_current[j])) | is.na(min(min(z[y == (j + 1)]), gamma_current[j+1]))) {
      print(stc("j", j))
      stop()
    }
    gamma[j+1] <- runif(1, max(max(z[y == j]), gamma_current[j]), min(min(z[y == (j + 1)]), gamma_current[j+1]))
  }
  gamma
}

for (i in 2:n_samp) {
  set.seed(i)
  z <- sample_Z(X, y, betas[i-1,], gammas[i-1,])
  gammas[i, ] <- sample_gamma(z, y, gammas[i-1, ])
  betas[i, ] <- rmvnorm(n = 1, betas[i-1,], sigma = xtx_inv)  
}


print(i)

# beta_post <- betas[(n_burnin+1):n_samp,] %>% 
beta_post <- betas %>% 
  as_tibble() %>% 
  drop_na()


beta_post %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = value, fill = name)) +
  facet_wrap(. ~ name, scales = "free", labeller = label_parsed) +
  geom_density() +
  theme_bw() +
  theme(legend.position="none") + 
  ggtitle("Posterior Distributions")
