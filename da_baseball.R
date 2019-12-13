dat <- read_csv("~/Desktop/LAD_assessment_data.csv",
                col_types = cols(batter_mlb_id = col_skip(), 
                                 pitcher_mlb_id = col_skip(), 
                                 venue_city = col_skip())) %>% drop_na()


X <- model.matrix(is_in_play ~ ., data = dat)
y <- dat$is_in_play

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

n_samp <- 5000
n_burnin <- 1000

beta <- matrix(NA, nrow = n_samp, ncol = p, dimnames = list(NULL, colnames(X)))



set.seed(230)
beta[1,] <- rnorm(p)

dat_z <- dat %>%
  select(-is_in_play) %>%
  mutate(z = generate_Z(y, X, beta[1,]))

sigma <- solve(t(X) %*% X)

for (i in 2:n_samp){
  beta_z <- lm(z ~ ., data = dat_z)$coeff
  beta[i,] <- rmvnorm(n = 1, mean = beta_z, sigma = sigma)
  dat_z$z <- generate_Z(y, X, beta[i,])
}

# Throw away burn-in
beta_post <- beta[(n_burnin+1):n_samp,] %>% 
  as_tibble()

colMeans(beta_post)

probit_model <- glm(is_in_play ~ ., family = binomial(link = "probit"), data = dat)
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

