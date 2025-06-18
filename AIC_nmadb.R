library(netmeta)
library(nmadb)
library(dplyr)

# Compute AIC of datasets in nmadb given indices using
#   multiplicative and random/fixed effect models
# 
compute_AIC <- function(dat = dat_nmadb,
                        ind, 
                        additive = TRUE) {
  len <- length(ind)
  AIC_add <- numeric(len)
  AIC_fixed <- numeric(len)
  AIC_multiplicative <- numeric(len)
  
  for (j in seq_along(ind)) {
    i <- ind[j]
    
    net <- runnetmeta(dat_nmadb$recid[i])
    #observed TE
    theta <- net$TE
    m <- net$m
    n <- net$n
    V <- diag(net$seTE^2)
    
    # Multiplicative effect
    theta_me <- net$TE.nma.common
    phi <- as.numeric(1/(m-n+1) * t(theta-theta_me) %*% 
                        solve(V) %*% (theta-theta_me))
    logL_me <- -0.5*(m*log(2*pi)+
                        log(det(phi*V))+ 
                        t(theta-theta_me) %*% solve(phi*V) %*% (theta-theta_me))
    AIC_me[j] <- 2*n-2*logL_me

    if (additive == TRUE) {
      # additive effects
      tau_hat <- net$tau
      sig <- V+tau_hat^2*diag(m)
      # fitted TE 
      theta_ae <- net$TE.nma.random
      logL_ae <- -0.5*(m*log(2*pi)+
                         log(det(sig))+ 
                         t(theta-theta_ae) %*% solve(sig) %*% (theta-theta_ae))
      AIC_ae[j] <- 2*n-2*logL_ae
      
    } else {
      # fixed effect
      sig <- V
      theta_fe <- net$TE.nma.fixed
      logL_fe <- -0.5*(m*log(2*pi)+
                         log(det(V))+ 
                         t(theta-theta_fe) %*% solve(V) %*% (theta-theta_fe))
      AIC_fe[j] <- 2*(n-1)-2*logL_fe
    }
    
  }
  
  if (additive == TRUE ) {
    return(list(
      AIC_add  = AIC_ae,
      AIC_mul  = AIC_me,
      AIC_diff = AIC_me - AIC_ae
    ))
  } else {
    return(list(
      AIC_fixed  = AIC_fe,
      AIC_mul  = AIC_me,
      AIC_diff = AIC_me - AIC_fe
    ))
  }
}

#TEST
# AIC_rr_r <- AIC_nmadb(ind_rr_r)
# 
# hist(AIC_rr_r$AIC_diff, 
#      main = "AIC_add – AIC_mul (risk ratio)",
#      xlab = "AIC difference",
#      breaks = seq(floor(min(AIC_rr_r$AIC_diff)), ceiling(max(AIC_rr_r$AIC_diff)), by = 1))


