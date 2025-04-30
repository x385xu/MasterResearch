library(netmeta)
source("generate_2arm.R")

#=====================================================
# S1 AE:
# d1=0, d2=1, d3=2
# tau = c(0.05, 0.25, 1)
# se(y) = c(0.25, 0.5, 1) for all studies
# # of 2 arm studies = c(2, 5, 10) for all contrasts 

X <- matrix(c(rep(c(1,0),2), rep(c(0,1),2)), nrow = 4, byrow = TRUE)
colnames(X) <- c("B", "C")
X
se <- rep(0.25,4)
d <- c(-2,2)

phi1 <- 1
tau1 <- 0.05
data1 <- generate_2arm(X,d,se,phi1,tau1,ref="A")
data1

# ae model
fitae <- netmeta(TE, seTE, treat1, treat2, studlab, data = data1)
# estimate of tau
fitae$tau

# me model
fitme <- lm(data1$TE ~ X-1, weights = se^2)
# estimate of phi
(summary(fitme)$sigma)^2



