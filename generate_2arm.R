
# X: m*n-1 design matrix
# se: m-dimensional vector, standard errors
# d: (n-1)-dimensional vector, true effect sizes
# tau: scalar, additive parameter, can also handle m-dim vector values
# phi: scalar, multiplicative parameter
# ref: name of the reference treatment

# generate data that can be passed to netmeta() directly
# assume a fixed-effect model by default (tau = 0, phi = 1)
# output  treat1, treat2, studlab, seTE, TE

generate_2arm <- function(X, d, se, phi = 1, tau = 0, ref = "ref") {
  
  m <- nrow(X)   # number of pairwise comparisons
  n <- ncol(X)   # number of treatments
  
  treatments <- colnames(X)   # treatment names
  # apply(X,1,fun) apply fun to each row of X
  trt1 <- apply(X, 1, function(row) {
    treatments[row == 1]
  })
  trt2 <- rep(ref, m)
  nmadata <- as.data.frame(cbind(trt1, trt2))
  colnames(nmadata) <- c("treat1", "treat2")
  
  # study labels
  nmadata$studlab <- paste0("Study", 1:m)
  
  # standard errors 
  nmadata$seTE <- se
  
  # effect sizes
  theta <- X %*% d
  TE <- rnorm(m, mean = theta, sd = sqrt(phi * se^2 + tau^2))
  nmadata$TE <- TE
  
  return(nmadata)
}
