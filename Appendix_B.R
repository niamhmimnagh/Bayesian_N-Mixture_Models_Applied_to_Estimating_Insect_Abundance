library(R2jags)
library(clusterGeneration)
library(mvtnorm)
library(extraDistr)

# Simulate synthetic data ------------------
set.seed(1234)
R <- 40                    # number of sites
T <- 10                    # number of sampling occasions
S <- 4                     # number of species
K <- 2                     # number of years
theta <- 0.3               # probability of a site being unoccupied
p <- 0.8                   # probability of detection
phi <- runif(S, -0.5, 0.5) # autocorrelation coefficient

# Empty arrays to store the observed and latent counts
y <- array(NA, dim=c(R,T,S,K))
N <- lambda<-array(NA, dim=c(R,S,K))
# Occupancy array
Occ <- array(rbinom(R*S*K, size=1, prob=1-theta), dim=c(R,S,K))
# Scale matrix for wishart distribution
Omega <- diag(1, nrow=S, ncol=S)

covariance <- genPositiveDefMat(S, rangeVar=c(0.2, 1), 
                                covMethod="unifcorrmat")[["Sigma"]] 
correlation <- cov2cor(covariance)
# species-level MVN random effect
a <- rmvnorm(R, mean=rep(0,S), sigma=covariance)

# Generate the latent abundance, N[i,s,k]
for(i in 1:R){
    for(s in 1:S){
        # for year K = 1     
        lambda[i,s,1]<-exp(a[i,s])
        N[i,s,1] <- ifelse(Occ[i,s,1]==0, 0,
                           rtpois(1, lambda=lambda[i,s,1], a=0))
        # for year K > 1
        for(k in 2:K){
           lambda[i,s,k] <- exp(a[i,s]+ phi[s]*log(N[i,s,k-1]+1))   
           N[i,s,k] <- ifelse(Occ[i,s,k]==0, 0, 
                              rtpois(n=1, lambda=lambda[i,s,k], a=0))
        }
    }
}

# Generate the observed abundance, y[i,s,k]
  for(i in 1:R){
    for(t in 1:T){
      for(s in 1:S){
        for(k in 1:K){
          y[i,t,s,k] <- rbinom(1, size=N[i,s,k], prob=p)
        }
      }
    }
  }
  
# Specify the model -------------------------
model_code = 'model {
for(s in 1:S){
    for (i in 1:R) { 
        Occ[i,s,1] ~ dbern(1-theta)
        log(lambda[i,s,1]) <- a[i,s]
        C[i,s,1] ~ dpois(lambda[i,s,1])T(1,)  
        N[i,s,1] <- ifelse(Occ[i,s,1]==0, 0, C[i,s,1])
        for(k in 2:K){
            Occ[i,s,k] ~ dbern(1-theta)
            log(lambda[i,s,k]) <- a[i,s] + phi[s]*log(N[i,s,k-1]+1)
            C[i,s,k] ~ dpois(lambda[i,s,k])T(1,)  
            N[i,s,k] <- ifelse(Occ[i,s,k]==0, 0, C[i,s,k])
            }
        }
    }
for(i in 1:R){
    for(s in 1:S){
        for(t in 1:T){   
            for(k in 1:K){
                y[i,t,s,k] ~ dbin(p, N[i,s,k])
            }
        }
    }
}

# species-level random effect
for(i in 1:R){
    a[i,1:S] ~ dmnorm(rep(0,S), precision[,])
}
  
# Inter-species correlations
precision[1:S,1:S] ~ dwish(Omega[,], df)
covariance[1:S,1:S] <- inverse(precision[,])
for (s in 1:S){    
    for (s1 in 1:S){
        correlation[s,s1] <- covariance[s,s1]/sqrt(covariance[s,s]*covariance[s1, s1])
    }
}
 
# Priors
theta~dbeta(1,1)
p~dunif(0,1)
for(s in 1:S){  
    phi[s] ~ dnorm(0,0.01)
}
}
' # end of model specification

# Package data for R2jags
data_list <- list(R=R, y=y, T=T, S=S, K=K, Omega=diag(1, S), df=S+1)
# Initial values
initial_values <- function(){
 list(C = apply(y,c(1,3,4), max)+1, 
     Occ = apply(y, c(1,3,4), function(z) ifelse(any(z>0), 1, 0)))
}  
# parameters to monitor
par_save=c("correlation", "N", "p", "theta", "Occ")

# run model
model_run <- jags(data               = data_list,
                  inits              = initial_values,
                  parameters.to.save = par_save,
                  model.file         = textConnection(model_code),
                  n.chains           = 4,
                  n.iter             = 25000,
                  n.burn             = 10000)
