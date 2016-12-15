## compute cox regression on the coefficients for survival analysis
library(survival)
library(broom)
load('./nhanes_schema_merged_all_031816.Rdata')

MainTable <- bigData
VarDescription <- tabDesc
## create a new 'area' variable that is a combination of the SDMVPSU and SDMVSTRA
MainTable$area <- paste(MainTable$SDMVPSU, MainTable$SDMVSTRA, sep='_')

rm(bigData)
rm(tabDesc)

baseline.mod <- coxph(Surv(PERMTH_EXM, MORTSTAT) ~ 
                        I(scale(RIDAGEYR))  + I(scale(RIDAGEYR)^2) + RIAGENDR + scale(INDFMPIR) + factor(RIDRETH1, c(3, 1, 2, 4, 5)) + cluster(area), 
                      weights=MainTable$WTINT2YR, MainTable)

tidyBaseline <- tidy(baseline.mod, exponentiate = T)
tidyBaseline$HR <- sprintf('%.02f', tidyBaseline$estimate)
tidyBaseline$conf_int_string <- (sprintf('%.02f,%.02f', tidyBaseline$conf.low, tidyBaseline$conf.high))
write.csv(tidyBaseline[, c('term', 'HR', 'conf_int_string', 'p.value')], row.names = F, file='~/Downloads/baseline.csv')
