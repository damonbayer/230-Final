library(tidyverse)
library(truncnorm)
library(mvtnorm)
dat <- read_csv("Data/breast-cancer-wisconsin.data", 
                col_names = c("ID", "thick","uni_size",
                              "uni_shape","adhesion",
                              "single_size","nuclei",
                              "chrom","nucleoli","mitoses",
                              "class"))
dat$nuclei <- as.numeric(dat$nuclei)
dat$class <- as.factor(dat$class)
X <- as.matrix(cbind(1,dat[,c(2:6,8:10)]))
y <- dat$class

generate_Z <- function(y,X,beta){
  if(y == 4){
    Z <- rtruncnorm(1,a=0,mean = X %*% beta,sd=1)
  } 
  if(y == 2){
    Z <- rtruncnorm(1,b=0,mean = X %*% beta,sd=1)
  }
  return(Z)
}

n_samp <- 5000
n_burnin <- 1000

beta <- matrix(NA,nrow = n_samp, ncol = ncol(X))
beta[1,] <- rmvnorm(1,rep(0,ncol(X)))
Z <- unlist(lapply(1:nrow(X), function(k){generate_Z(y[k],X[k,],beta[1,])}))

for (i in 2:n_samp){
  beta_z <- lm(Z~X-1)$coeff
  sigma <- solve(t(X) %*% X)
  beta[i,] <- rmvnorm(1,mean=beta_z,sigma=sigma)
  Z <- unlist(lapply(1:nrow(X), function(k){generate_Z(y[k],X[k,],beta[i,])}))
}
colMeans(beta[(n_burnin+1):n_samp,])

probit_model <- glm(y~X-1,family = binomial(link = "probit"))
summary(probit_model)


beta[(n_burnin+1):n_samp,] %>% 
  as_tibble() %>% 
  `colnames<-`(colnames(X)) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = value, fill = name)) +
  facet_wrap(. ~ name, scales = "free") +
  geom_density() +
  theme_bw() +
  theme(legend.position="none")