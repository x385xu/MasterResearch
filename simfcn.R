#=====================================================
library(netmeta)
source("generate_2arm.R")
#=====================================================
# Function to simulate dat
# X, d, se, phi, tau, ref are the sample arguments as in generate_2arm
# se has to be a vector 
# nsim: number of simulations
# return estimated phis and taus for additive and multiplicative models respectively
simfcn <- function(X, d, se, phi = 1, tau = 0, ref = "ref", nsim) {
  taus <- rep(0, nsim)
  phis <- rep(0, nsim)
  
  for (i in 1:nsim){
    simdata <- generate_2arm(X, d, se, phi, tau, ref)
    net <- netmeta(TE, seTE, treat1, treat2, studlab, data = simdata)
    taus[i]<-net$tau
    
    fitme <- lm(simdata$TE ~ X-1, weights = 1/(se^2))
    phis[i] <- (summary(fitme)$sigma)^2
  }
  
  # average, relative bias, MSE, variance
  taustats <- c(mean(taus), mean(taus-tau), mean((taus-tau)^2), var(taus))
  phistats <- c(mean(phis), mean(phis-phi), mean((phis-phi)^2), var(phis))
  
  results <- as.data.frame(matrix(c(taustats, phistats), nrow = 1))
  colnames(results) <- c("tau_mean", "tau_bias", "tau_mse", "tau_var", 
                         "phi_mean", "phi_bias", "phi_mse", "phi_var")
  
  return(results)
}


# nsim <- 10
# X <- matrix(c(rep(c(1,0),2), rep(c(0,1),2)), nrow = 4, byrow = TRUE)
# colnames(X) <- c("B", "C")
# se <- rep(0.25,4)
# d <- c(1,2)
# 
# phi1 <- 1
# tau1 <- 0.05
# results1 <- simfcn(X,d,se,phi1,tau1,ref="A",nsim)
# 
# phi2 <- 1.5
# tau2 <- 0
# results2 <- simfcn(X,d,se,phi2,tau2,ref="A",nsim)
