#=====================================================
library(netmeta)
source("generate_2arm.R")
source("simfcn.R")
#=====================================================
nsim <- 100
#=====================================================
# S1 AE:
# d1=0, d2=1, d3=2
# tau = c(0.05, 0.25, 1)
# phi = 1
# se(y) = c(0.25, 0.5, 1) for all studiesre
# # of 2 arm studies = c(2, 5, 10) for all contrasts 

d <- c(1,2)
ses  <- c(0.25, 0.5, 1)
n2s   <- c(2, 5, 10)

taus <- c(0.05, 0.25, 1)

results_list_s1 <- vector("list", 
                       length = length(taus)*length(ses)*length(n2s))
idx <- 1

for (tau_i in taus) {
  for (se_i in ses) {
    for (n2 in n2s) {
      X <- matrix(c(rep(c(1,0),n2), rep(c(0,1),n2)), nrow = n2*2, byrow = TRUE)
      colnames(X) <- c("B", "C")
      se_vec <- rep(se_i, nrow(X))
      
      out <- simfcn(X, d, se_vec,
                    phi = 1, tau = tau_i,
                    ref = "A", nsim)
      
      # store with labels
      results_list_s1[[idx]] <- data.frame(
        tau = tau_i,
        se  = se_i,
        n2  = n2,
        out,
        phi_theo = 1+(tau_i^2)/(se_i^2)
      )
      idx <- idx + 1
    }
  }
}

results_s1 <- do.call(rbind, results_list_s1)

#=====================================================
# S2 ME:
# d1=0, d2=1, d3=2
# tau = 0
# phi = c(1.5, 2, 5)
# se(y) = c(0.25, 0.5, 1) for all studies
# # of 2 arm studies = c(2, 5, 10) for all contrasts 

d <- c(1,2)
ses  <- c(0.25, 0.5, 1)
n2s   <- c(2, 5, 10)

phis <- c(1.5, 2, 5)

results_list_s2 <- vector("list", 
                       length = length(phis)*length(ses)*length(n2s))
idx <- 1

for (phi_i in phis) {
  for (se_i in ses) {
    for (n2 in n2s) {
      X <- matrix(c(rep(c(1,0),n2), rep(c(0,1),n2)), nrow = n2*2, byrow = TRUE)
      colnames(X) <- c("B", "C")
      se_vec <- rep(se_i, nrow(X))
      
      out <- simfcn(X, d, se_vec,
                    phi = phi_i, tau = 0,
                    ref = "A", nsim)
      
      # store with labels
      results_list_s2[[idx]] <- data.frame(
        phi = phi_i,
        se  = se_i,
        n2  = n2,
        out,
        tau_theo = se_i*sqrt(phi_i-1)
      )
      idx <- idx + 1
    }
  }
}

results_s2 <- do.call(rbind, results_list_s2)

