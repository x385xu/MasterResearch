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
         Primary.Outcome, Description.of.the.outcome.s., 
         Harmful.Beneficial.Outcome, dataset) %>%
  rename(recid = Record.ID) %>%
  mutate(Year = as.numeric(format(as.Date(Year, format="%Y-%m-%d"),"%Y")), 
         .keep = "unused") 

#nma_test <- runnetmeta(dat_nmadb$recid[2])

source("get_index_nmadb.R")
source("compute_AIC.R")

#===============================================================================
#AIC_mul_add = AIC_mul - AIC_add,
#AIC_mul_fixed = AIC_mul - AIC_fixed
#ΔAIC=AIC_add-AIC_fixed
# ΔAIC < -3 -----additive-1
# -3 < ΔAIC < 3 ---fixed-2
# ΔAIC > 3 --------fixed-3


#---------or-------------------------------------------------------------------
ind_or <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "odds ratio",
                           multiarm = FALSE,
                           fixed = TRUE)
AIC_or <- compute_AIC(dat = dat_nmadb,ind_or)

diffs_or <- AIC_or$add - AIC_or$fixed
AICdiff_or_1 <- AIC_or$mul_add[which(diffs_or < -3)]
AICdiff_or_2 <- AIC_or$mul_fixed[which(abs(diffs_or) <= 3)]
AICdiff_or_3 <- AIC_or$mul_fixed[which(diffs_or > 3)]


hist(AICdiff_or_1,
     main = "AIC difference",
     breaks = seq(floor(min(AICdiff_or_1)), 
                  ceiling(max(AICdiff_or_1)), by = 1))
abline(v = c(-3,3), col = "red", lty = 2)

#-------------rr---------------------------------------------------------------
ind_rr <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "risk ratio",
                           multiarm = FALSE,
                           fixed = TRUE)
AIC_rr <- compute_AIC(dat = dat_nmadb,ind_rr)

diffs_rr <- AIC_rr$add - AIC_rr$fixed
AICdiff_rr_1 <- AIC_rr$mul_add[which(diffs_rr < -3)]
AICdiff_rr_2 <- AIC_rr$mul_fixed[which(abs(diffs_rr) <= 3)]
AICdiff_rr_3 <- AIC_rr$mul_fixed[which(diffs_rr > 3)]

#------------md----------------------------------------------------------------
ind_md <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "mean difference",
                           multiarm = FALSE,
                           fixed = TRUE)
AIC_md <- compute_AIC(dat = dat_nmadb,ind_md)

diffs_md <- AIC_md$add - AIC_md$fixed
AICdiff_md_1 <- AIC_md$mul_add[which(diffs_md < -3)]
AICdiff_md_2 <- AIC_md$mul_fixed[which(abs(diffs_md) <= 3)]
AICdiff_md_3 <- AIC_md$mul_fixed[which(diffs_md > 3)]

#===========plot===============================================================
library(ggplot2)

# build a data-frame with one row per observation
df <- bind_rows(
  tibble(value = AICdiff_or_1, category = "Additive better", measure = "OR"),
  tibble(value = AICdiff_or_2, category = "Similar",         measure = "OR"),
  tibble(value = AICdiff_or_3, category = "Fixed better",    measure = "OR"),
  tibble(value = AICdiff_rr_1, category = "Additive better", measure = "RR"),
  tibble(value = AICdiff_rr_2, category = "Similar",         measure = "RR"),
  tibble(value = AICdiff_rr_3, category = "Fixed better",    measure = "RR"),
  tibble(value = AICdiff_md_1, category = "Additive better", measure = "MD"),
  tibble(value = AICdiff_md_2, category = "Similar",         measure = "MD"),
  tibble(value = AICdiff_md_3, category = "Fixed better",    measure = "MD")
) %>%
  mutate(
    # specify the column and row order
    category = factor(category,
                      levels = c("Additive better",
                                 "Similar",
                                 "Fixed better")),
    measure  = factor(measure,  levels = c("OR", "RR", "MD")))


ggplot(df, aes(x = value)) + 
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "white") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  facet_grid(measure ~ category, scales = "free") +
  labs(x = expression(Delta~"AIC (Multiplicative-Additive or Fixed)"), 
       y = "Count") +
  theme_minimal(base_size = 12)






