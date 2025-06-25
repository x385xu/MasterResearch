# Files 
`nmadb_AIC.R`: File that computes the AIC values for additive random-effects model, fixed effect model and multiplicative effect model and plot AIC differences

`simulation draft.R`: where to run the simulation

`simulation2.md`: examples of simulation, additive random-effects model, fixed effect model and multiplicative effect model, also use it to do small tests

# Functions
## `compute_AIC()`
Function that that computes the AIC values for additive random-effects model, fixed effect model and multiplicative effect model given the indices of datasets in nmadb

Arguments: 

    dat=dat_nmadb: database    
    ind: indices
Value:

    ind: indices
    recid: recid in nmadb
    add: AIC under additive random-effects model
    fixed: AIC under fixed effect model
    mul: AIC under multiplicative effects model
    mul_add: AIC difference, multiplicative model - additive model
    mul_fixed: AIC difference, multiplicative model - fixed effect model

## `get_index_nmadb()` 
Get index of datasets in nmadb based on the measure

Arguments: 

     measure: string, like "risk ratio" or "odds ratio"
     multiarm: TRUE = allow multiarm studies; FALSE = only include 2-arm studies
     fixed: TRUE = allow fixed-effect models; FALSE = only include additive random effect models
Value:

    ind: a sequence of indices
## `Generate_2arm`
generate NMA data, only for 2-arm studies

