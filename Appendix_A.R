library(rjags)
library(R2jags)


# Simulate synthetic data -------------------
set.seed(1234)
T <- 10   # Number of periods of observation
R <- 100  # Number of sites
p <- 0.7  # Detection probability
y <- matrix(data = NA,ncol = T,nrow = R)
lambda <- 20
N <- rpois(n = R, lambda)
for (t in 1:T){
  y[,t] <- rbinom(n = R,size = N, prob = p)
}


# Specify the model -------------------------
model_code = '
model{
for (i in 1:R) {
    N[i] ~ dpois(lambda)
        for (t in 1:T) {
            y[i,t] ~ dbin(p, N[i])
        }
}
# Priors
p~dunif(0, 1)
lambda~dgamma(1, 0.1)
}'
# Package data for R2jags
data_list <- list(R = R, T = T, y = y)
# Initial values
initial_values <- function(){list(N = apply(y, 1, max))}
# parameters to monitor
par_save <- c("p", "lambda")
# run model
model_run <- jags(data               = data_list,
                  inits              = initial_values,
                  parameters.to.save = par_save,
                  model.file         = textConnection(model_code),
                  n.chains           = 4,
                  n.iter             = 10000,
                  n.burn             = 5000)
