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

#nma_test <- runnetmeta(dat_nmadb$recid[2])

source("get_index_nmadb.R")
source("compute_AIC.R")

#=====risk ratio=========================================================================

ind_rr_a <- get_index_nmadb(dat = dat_nmadb,
                            measure = "risk ratio",
                            multiarm = FALSE,
                            fixed = FALSE)
AIC_rr_a <- compute_AIC(dat = dat_nmadb,ind_rr_a)

ind_rr_f <- get_index_nmadb(dat = dat_nmadb,
                            measure = "risk ratio",
                            multiarm = FALSE,
                            fixed = TRUE)
AIC_rr_f <- compute_AIC(dat = dat_nmadb,ind_rr_f)

par(mfrow = c(1, 2))
hist(AIC_rr_a$AIC_mul_add,
     main = "AIC difference - risk ratio",
     xlab = "AIC difference (multiplicative-additive)",
     breaks = seq(floor(min(AIC_rr_a$AIC_mul_add)), 
                  ceiling(max(AIC_rr_a$AIC_mul_add)), by = 1))
abline(v = c(-3,3), col = "red", lty = 2)

hist(AIC_rr_f$AIC_mul_fixed,
     main = "AIC difference - risk ratio",
     xlab = "AIC difference (multiplicative-fixed)",
     breaks = seq(floor(min(AIC_rr_f$AIC_mul_fixed)), 
                  ceiling(max(AIC_rr_f$AIC_mul_fixed)), by = 1))
abline(v = c(-3,3), col = "red", lty = 2)


#======odds ratio=======================================================================

ind_or_a <- get_index_nmadb(dat = dat_nmadb,
                            measure = "odds ratio",
                            multiarm = FALSE,
                            fixed = FALSE)
AIC_or_a <- compute_AIC(dat = dat_nmadb,ind_or_a)

ind_or_f <- get_index_nmadb(dat = dat_nmadb,
                            measure = "odds ratio",
                            multiarm = FALSE,
                            fixed = TRUE)
AIC_or_f <- compute_AIC(dat = dat_nmadb,ind_or_f)

par(mfrow = c(1, 2))

hist(AIC_or_a$AIC_mul_add,
     main = "AIC difference - odds ratio",
     xlab = "AIC difference (multiplicative-additive)",
     breaks = seq(floor(min(AIC_or_a$AIC_mul_add)), 
                  ceiling(max(AIC_or_a$AIC_mul_add)), by = 1))
abline(v = c(-3,3), col = "red", lty = 2)

hist(AIC_or_f$AIC_mul_fixed,
     main = "AIC difference - odds ratio",
     xlab = "AIC difference (multiplicative-fixed)",
     breaks = seq(floor(min(AIC_or_f$AIC_mul_fixed)), 
                  ceiling(max(AIC_or_f$AIC_mul_fixed)), by = 1))
abline(v = c(-3,3), col = "red", lty = 2)

#========mean difference =======================================================
ind_md_a <- get_index_nmadb(dat = dat_nmadb,
                            measure = "mean difference",
                            multiarm = FALSE,
                            fixed = FALSE)
AIC_md_a <- compute_AIC(dat = dat_nmadb,ind_md_a)

ind_md_f <- get_index_nmadb(dat = dat_nmadb,
                            measure = "mean difference",
                            multiarm = FALSE,
                            fixed = TRUE)
AIC_md_f <- compute_AIC(dat = dat_nmadb,ind_md_f)

par(mfrow = c(1, 2))

hist(AIC_md_a$AIC_mul_add,
     main = "AIC difference - mean difference",
     xlab = "AIC difference (multiplicative-additive)",
     breaks = seq(floor(min(AIC_md_a$AIC_mul_add)), 
                  ceiling(max(AIC_md_a$AIC_mul_add)), by = 1))
abline(v = c(-3,3), col = "red", lty = 2)

hist(AIC_md_f$AIC_mul_fixed,
     main = "AIC difference -
     mean difference",
     xlab = "AIC difference (multiplicative-fixed)",
     breaks = seq(floor(min(AIC_md_f$AIC_mul_fixed)), 
                  ceiling(max(AIC_md_f$AIC_mul_fixed)), by = 1))

abline(v = c(-3,3), col = "red", lty = 2)


#============================================
AIC_md_f$AIC_mul_fixed[5]
AIC_md_f$AIC_mul_add

net <- runnetmeta(dat_nmadb$recid[ind_md_f[5]])

theta <- net$TE
m <- net$m
n <- net$n
#multiplicative effect
V <- diag(net$seTE^2)
theta_m <- net$TE.nma.fixed
phi <- as.numeric(1/(m-n+1) * t(theta-theta_m) %*% 
                    solve(V) %*% (theta-theta_m))

logL_m <- -0.5*(m*log(2*pi)+
                   log(det(phi*V))+ 
                   t(theta-theta_m) %*% solve(phi*V) %*% (theta-theta_m))
AIC_m <- 2*n-2*logL_m

# fixed
sig <- V
theta_f <- net$TE.nma.fixed
logL_f <- -0.5*(m*log(2*pi)+
                   log(det(V))+ 
                   t(theta-theta_f) %*% solve(V) %*% (theta-theta_f))
AIC_f <- 2*(n-1)-2*logL_f

# additive
tau_hat <- net$tau
sig <- V+tau_hat^2*diag(m)
theta_a <- net$TE.nma.random
logL_a <- -0.5*(m*log(2*pi)+
                   log(det(sig))+ 
                   t(theta-theta_a) %*% solve(sig) %*% (theta-theta_a))
AIC_a <- 2*n-2*logL_a

#==plot===============================




