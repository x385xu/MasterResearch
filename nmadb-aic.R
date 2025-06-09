
#newest version of netmeta does not have the function pairwise() anymore
#remove.packages("netmeta")
#install.packages("netmeta_2.9-0.tar.gz", repos = NULL, type = "source")
library(netmeta)


#nmadb is removed from cran, can find it in the archive: https://cran.r-project.org/src/contrib/Archive/nmadb/?C=D;O=A
install.packages("nmadb_1.2.0.tar.gz", repos = NULL, type = "source")
library(nmadb)

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

#==========================================================================================
# Extract dataset which use risk ratio, random effects model, and only have two-arm studies
riskratio <- which(dat_nmadb$Effect.Measure== "risk ratio")

index <- c()
for (i in riskratio) {
  nma <- tryCatch({
    runnetmeta(dat_nmadb$recid[i])
  }, error = function(e) {
    cat("  Failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Skip if nma is NULL due to error
  if (!is.null(nma) && is.list(nma)) {
    if (nma$tau > 0 && all(nma$multiarm == FALSE)) {
        index <- c(index, i)
    }
  }
}


# Loop through the extracted datasets and calculate the AIC difference
AICdiff <- c()

for (i in index) {
  nma_test <- runnetmeta(dat_nmadb$recid[i])
  theta <- nma_test$TE
  
  sig <- diag(nma_test$seTE.nma.random^2)
  theta_add <- nma_test$TE.nma.random
  logL_add <- -0.5*(log(2*pi)*det(sig)+ t(theta-theta_add) %*% solve(sig) %*% (theta-theta_add))
  
  V <- diag(nma_test$seTE.nma.common^2)
  theta_mul <- nma_test$TE.nma.common
  m <- nma_test$m
  n <- nma_test$n
  phi <- as.numeric(1/(m-n+1) * t(theta-theta_mul) %*% solve(V) %*% (theta-theta_mul))
  logL_mul <- -0.5*(log(2*pi)*det(phi*V)+ t(theta-theta_mul) %*% solve(phi*V) %*% (theta-theta_mul))
  
  # AIC = 2k-2log(L), k is the number of estimated parameters
  k <- n  # d is n-1 dimensional + 1 (phi or tau) = n
  AIC_add <- 2*k-2*logL_add
  AIC_mul <- 2*k-2*logL_mul
  
  # AIC_add - AIC_mul
  AICdiff <- c(AICdiff, 2*logL_mul - 2*logL_add)
}

hist(AICdiff, main = "Histogram of AIC difference (risk ratio)")



#==============================================================================
# Extract dataset which use odds ratio, random effects model, and only have two-arm studies
oddsratio <- which(dat_nmadb$Effect.Measure== "odds ratio")

index <- c()
for (i in oddsratio) {
  nma <- tryCatch({
    runnetmeta(dat_nmadb$recid[i])
  }, error = function(e) {
    cat("  Failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # Skip if nma is NULL due to error
  if (!is.null(nma) && is.list(nma)) {
    if (nma$tau > 0 && all(nma$multiarm == FALSE)) {
      index <- c(index, i)
    }
  }
}

# Loop through the extracted datasets and calculate the AIC difference
AICdiff <- c()

for (i in index) {
  nma_test <- runnetmeta(dat_nmadb$recid[i])
  theta <- nma_test$TE
  
  sig <- diag(nma_test$seTE.nma.random^2)
  theta_add <- nma_test$TE.nma.random
  logL_add <- -0.5*(log(2*pi)*det(sig)+ t(theta-theta_add) %*% solve(sig) %*% (theta-theta_add))
  
  V <- diag(nma_test$seTE.nma.common^2)
  theta_mul <- nma_test$TE.nma.common
  m <- nma_test$m
  n <- nma_test$n
  phi <- as.numeric(1/(m-n+1) * t(theta-theta_mul) %*% solve(V) %*% (theta-theta_mul))
  logL_mul <- -0.5*(log(2*pi)*det(phi*V)+ t(theta-theta_mul) %*% solve(phi*V) %*% (theta-theta_mul))
  
  # AIC = 2k-2log(L), k is the number of estimated parameters
  k <- n  # d is n-1 dimensional + 1 (phi or tau) = n
  AIC_add <- 2*k-2*logL_add
  AIC_mul <- 2*k-2*logL_mul
  
  # AIC_add - AIC_mul
  AICdiff <- c(AICdiff, 2*logL_mul - 2*logL_add)
}

hist(AICdiff, main = "AIC_add-AIC_mul (odds ratio)")

#==========================================================================================
# Take an example of random-effect model
nma_test <- runnetmeta(dat_nmadb$recid[2])

theta <- nma_test$TE

# Additive model
sig <- diag(nma_test$seTE.nma.random^2)
theta_add <- nma_test$TE.nma.random

logL_add <- -0.5*(log(2*pi)*det(sig)+ t(theta-theta_add) %*% solve(sig) %*% (theta-theta_add))

# Multiplicative model
V <- diag(nma_test$seTE.nma.common^2)
theta_mul <- nma_test$TE.nma.common

m <- nma_test$m
n <- nma_test$n
phi <- as.numeric(1/(m-n+1) * t(theta-theta_mul) %*% solve(V) %*% (theta-theta_mul))

logL_mul <- -0.5*(log(2*pi)*det(phi*V)+ t(theta-theta_mul) %*% solve(phi*V) %*% (theta-theta_mul))


# AIC = 2k-2log(L), k is the number of estimated parameters
k <- n  # d is n-1 dimensional + 1 (phi or tau) = n
AIC_add <- 2*k-2*logL_add
AIC_mul <- 2*k-2*logL_mul
