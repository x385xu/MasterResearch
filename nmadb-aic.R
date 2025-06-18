# newest version of netmeta does not have the function pairwise() anymore
#remove.packages("netmeta")
#install.packages("netmeta_2.9-0.tar.gz", repos = NULL, type = "source")
library(netmeta)

# nmadb is removed from cran, can find it in the archive: https://cran.r-project.org/src/contrib/Archive/nmadb/?C=D;O=A
#remove.packages("nmadb")
#install.packages("nmadb_1.2.0.tar.gz", repos = NULL, type = "source")
library(nmadb)
library(dplyr)
# Load NMA database catalog
dat_nmadb <- getNMADB()

# Extract relevant info
dat_nmadb <- dat_nmadb %>% 
  select(Record.ID, Title, First.Author, Year, 
         Number.of.Studies.., Number.of.Treatments, 
         Type.of.Outcome., Effect.Measure, Fixed.effect.or.Random.effects.model,
         Primary.Outcome, Description.of.the.outcome.s., Harmful.Beneficial.Outcome,
         dataset) %>%
  rename(recid = Record.ID) %>%
  mutate(Year = as.numeric(format(as.Date(Year, format="%Y-%m-%d"),"%Y")), .keep = "unused") 

nma_test <- runnetmeta(dat_nmadb$recid[2])


#=====RR=====================================================================================
# Extract dataset which use risk ratio, random effects model, and only have two-arm studies
rr <- which(dat_nmadb$Effect.Measure== "risk ratio")

ind_rr <- c()
for (i in rr) {
  nma <- tryCatch({runnetmeta(dat_nmadb$recid[i])}, error = function(e) {
    return(NULL)
  })
  
  # Skip if nma is NULL due to error
  if (!is.null(nma) && is.list(nma)) {
    if (nma$tau > 0 && all(nma$multiarm == FALSE)) {
      ind_rr <- c(ind_rr, i)
    }
  }
}
len_rr <- length(ind_rr)

# Loop through the extracted datasets and calculate the AIC difference
AICdiff_rr <- numeric(len_rr)
AIC_add_rr <- numeric(len_rr)
AIC_mul_rr <- numeric(len_rr)

for (j in 1:len_rr) {
  i <- ind_rr[j]
  net <- runnetmeta(dat_nmadb$recid[i])
  # Observed TE
  theta <- net$TE
  m <- net$m
  n <- net$n
  tau_hat <- net$tau
  # Diagonal matrix of observed var (multiplicative model)
  V <- diag(net$seTE^2)
  # var matrix for additive model
  sig <- V+tau_hat^2*diag(m)
  
  # fitted TE 
  theta_add <- net$TE.nma.random
  logL_add <- -0.5*(m*log(2*pi)+
                      log(det(sig))+ 
                      t(theta-theta_add) %*% solve(sig) %*% (theta-theta_add))
  
  theta_mul <- net$TE.nma.common
  phi <- as.numeric(1/(m-n+1) * t(theta-theta_mul) %*% 
                      solve(V) %*% (theta-theta_mul))
  logL_mul <- -0.5*(m*log(2*pi)+
                      log(det(phi*V))+ 
                      t(theta-theta_mul) %*% solve(phi*V) %*% (theta-theta_mul))
  
  # AIC = 2k-2log(L), k is the number of estimated parameters
  # d is n-1 dimensional + 1 (phi or tau) = n
  AIC_add <- 2*n-2*logL_add
  AIC_mul <- 2*n-2*logL_mul
  
  AIC_add_rr[j] <- AIC_add
  AIC_mul_rr[j] <- AIC_mul
  
  # AIC_add - AIC_mul
  AICdiff_rr[j] <- AIC_add-AIC_mul
}

hist(AICdiff_rr,
     main = "AIC_add - AIC_mul (risk ratio)",
     xlab = "AIC difference (additive - multiplicative)",
     breaks = seq(floor(min(AICdiff_rr)), ceiling(max(AICdiff_rr)), by = 1))

abline(v = 0, col = "red", lty = 2)



#======OR========================================================================
# Extract dataset which use odds ratio, random effects model, and only have two-arm studies
or <- which(dat_nmadb$Effect.Measure== "odds ratio")

ind_or <- c()
for (i in or) {
  nma <- tryCatch({runnetmeta(dat_nmadb$recid[i])}, error = function(e) {
    return(NULL)
  })
  
  # Skip if nma is NULL due to error
  if (!is.null(nma) && is.list(nma)) {
    if (nma$tau > 0 && all(nma$multiarm == FALSE)) {
      ind_or <- c(ind_or, i)
    }
  }
}

len_or <- length(ind_or)

# Loop through the extracted datasets and calculate the AIC difference
AICdiff_or <- numeric(len_or)
AIC_add_or <- numeric(len_or)
AIC_mul_or <- numeric(len_or)

for (j in 1:len_or) {
  i <- ind_or[j]
  net <- runnetmeta(dat_nmadb$recid[i])
  # Observed TE
  theta <- net$TE
  m <- net$m
  n <- net$n
  tau_hat <- net$tau
  # Diagonal matrix of observed var (multiplicative model)
  V <- diag(net$seTE^2)
  # var matrix for additive model
  sig <- V+tau_hat^2*diag(m)
  
  # fitted TE 
  theta_add <- net$TE.nma.random
  logL_add <- -0.5*(m*log(2*pi)+
                      log(det(sig))+ 
                      t(theta-theta_add) %*% solve(sig) %*% (theta-theta_add))
  
  theta_mul <- net$TE.nma.common
  phi <- as.numeric(1/(m-n+1) * t(theta-theta_mul) %*% 
                      solve(V) %*% (theta-theta_mul))
  logL_mul <- -0.5*(m*log(2*pi)+
                      log(det(phi*V))+ 
                      t(theta-theta_mul) %*% solve(phi*V) %*% (theta-theta_mul))
  
  # AIC = 2k-2log(L), k is the number of estimated parameters
  # d is n-1 dimensional + 1 (phi or tau) = n
  AIC_add <- 2*n-2*logL_add
  AIC_mul <- 2*n-2*logL_mul
  
  AIC_add_or[j] <- AIC_add
  AIC_mul_or[j] <- AIC_mul
  
  # AIC_add - AIC_mul
  AICdiff_or[j] <- AIC_add-AIC_mul
}

hist(AICdiff_or,
     main = "AIC_add - AIC_mul (odds ratio)",
     xlab = "AIC difference (additive - multiplicative)",
     breaks = seq(floor(min(AICdiff_or)), ceiling(max(AICdiff_or)), by = 1))

abline(v = 0, col = "red", lty = 2)
