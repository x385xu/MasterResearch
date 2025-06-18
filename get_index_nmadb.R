library(netmeta)
library(nmadb)
library(dplyr)


# Extract index of datasets in nmadb based on the measure,
#   whether it's multiarm, and whether it uses a random effect model
# measure: string, like "risk ratio" or "odds ratio"
# multiarm: TRUE = allow multiarm studies; FALSE = only include 2-arm studies
# fixed: TRUE = allow fixed-effect models
get_index_nmadb <- function(measure, 
                            multiarm = FALSE, 
                            random = FALSE) {
  # Extract index with the according measure
  ind_measure <- which(dat_nmadb$Effect.Measure== measure)
  
  ind <- c()
  for (i in ind_measure) {
    # return null if runnetmeta cannot get the dataset
    net <- tryCatch({runnetmeta(dat_nmadb$recid[i])}, 
                    error = function(e) {return(NULL)})
    
    # Skip if nma is NULL due to error
    if (is.null(net) && !is.list(net)) next
    
    has_multi <- any(net$multiarm) # is FALSE if every entry is FALSE
    if (multiarm != has_multi) next
    
    is_random <- (net$tau > 0)
    if (random != is_random) next
    
    ind <- c(ind, i)
  }
  return(ind)
}



#TEST
# # Load NMA database catalog
# dat_nmadb <- getNMADB()
# 
# # Extract relevant info
# dat_nmadb <- dat_nmadb %>% 
#   select(Record.ID, Title, First.Author, Year, 
#          Number.of.Studies.., Number.of.Treatments, 
#          Type.of.Outcome., Effect.Measure, Fixed.effect.or.Random.effects.model,
#          Primary.Outcome, Description.of.the.outcome.s., Harmful.Beneficial.Outcome,
#          dataset) %>%
#   rename(recid = Record.ID) %>%
#   mutate(Year = as.numeric(format(as.Date(Year, format="%Y-%m-%d"),"%Y")), .keep = "unused") 
# 
# 
# ind_rr_r <- get_index_nmadb(dat = dat_nmadb, 
#                             measure = "risk ratio", 
#                             multiarm = FALSE, 
#                             random = TRUE)
