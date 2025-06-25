library(netmeta)
library(nmadb)
library(dplyr)


# Extract index of datasets in nmadb based on the measure,
#   whether it's multiarm, and whether it uses a random effect model
# measure: string, like "risk ratio" or "odds ratio"
# multiarm: TRUE = allow multiarm studies; FALSE = only include 2-arm studies
# fixed: TRUE = allow fixed-effect models; FALSE = only include additive effect models
get_index_nmadb <- function(dat = dat_nmadb,
                            measure, 
                            multiarm = FALSE, 
                            fixed = FALSE) {
  # Extract index with the according measure
  ind_measure <- which(dat_nmadb$Effect.Measure == measure)
  
  ind <- c()
  for (i in ind_measure) {
    # return null if runnetmeta cannot get the dataset
    net <- tryCatch({runnetmeta(dat_nmadb$recid[i])}, 
                    error = function(e) {return(NULL)})
    
    # Skip if net is not a list
    if (!is.list(net)) next
    
    has_multi <- any(net$multiarm) # is FALSE if every entry is FALSE (no multiarm studies)
    # Skip multiarm study when multi=FALSE
    if ((!multiarm) && has_multi) next
    
    is_fixed <- (net$tau == 0)
    # Skip fixed effect model when fixed=FALSE
    if ((!fixed) && is_fixed) next
    
    ind <- c(ind, i)
  }
  return(ind)
}


