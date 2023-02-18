# Simulate synthetic data -------------------
set.seed(1234)
T <- 5  # Number of periods of observation
R <- 20  # Number of sites
p <- 0.1  # Detection probability
y <- matrix(data = NA,ncol = T,nrow = R)
lambda <- 2

N <- rpois(n = R, lambda)
for (t in 1:T){
  y[,t] <- rbinom(n = R,size = N, prob = p)
}

# Diagnostic for the N-mixture model --------
CovarianceDiagnostic<-function(y, T){
  ninj <- 0
  for(i in 1:(T-1)){
    for(j in (i+1):T){
      ninj<-cbind(ninj,y[,i]*y[,j])
    }
  }
  covDiag <- sum(colMeans(ninj))*2/(T*(T-1))-((sum(colMeans(y)))/T)^2
  print(paste0("Covariance Diagnostic: ", round(covDiag,4)))
}

# Run diagnostic
CovarianceDiagnostic(y, T)
