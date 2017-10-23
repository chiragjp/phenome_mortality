## examine the hospital utilization with respect to the predictors of mortality:
## 1.) does hospital use modify the associations?
## 2.) is hospital use correlated with the phenotypes?

library(survey, quietly = T)     # required for stratified analysis
options(survey.lonely.psu="adjust")
library(getopt, quietly = T)


spec <- matrix(c(
  'survey_number', 'n', 2, 'numeric'
),
  nrow=1, byrow=TRUE)

opt <- getopt(spec)

SDDSRVYR_INPUT <- NA
if(!is.null(opt$survey_number)) {
  SDDSRVYR_INPUT <- opt$survey_number
}

cat(sprintf('loading data...'))
load('./nhanes_schema_merged_all_031816.Rdata')
finalVariableList <- read.csv('FinalVariables_081517.csv', stringsAsFactors = F) # this was output from Nam's Supplementary Table
cat(sprintf('...done\n'))
MainTable <- bigData
VarDescription <- tabDesc
rm(bigData)
rm(tabDesc)

## create a new 'area' variable that is a combination of the SDMVPSU and SDMVSTRA
MainTable$area <- paste(MainTable$SDMVPSU, MainTable$SDMVSTRA, sep='_')

## code race and ethnicity
MainTable$black <- ifelse(MainTable$RIDRETH1 == 4, 1, 0)
MainTable$mexican <- ifelse(MainTable$RIDRETH1 == 1, 1, 0)
MainTable$other_eth <- ifelse(MainTable$RIDRETH1 == 5, 1, 0)
MainTable$other_hispanic <- ifelse(MainTable$RIDRETH1 == 2, 1, 0)
MainTable$white <- ifelse(MainTable$RIDRETH1 == 3, 1, 0)

# 
source('~/Dropbox (HMS)/RagGroup Team Folder/nhanes_quant_pheno/util.R')


###### encode hospital use and health utilization
## general health condition
#table(MainTable$HUQ010, MainTable$SDDSRVYR) # for A, B, C, D, E, F, G,H: General health condition
MainTable$general_health_condition <- ifelse(MainTable$HUQ010 < 7, MainTable$HUQ010, NA)
#table(MainTable$general_health_condition, MainTable$SDDSRVYR) # for A, B, C, D, E, F, G,H: General health condition


## number of times recieve healthcare
#table(MainTable$HUQ050, MainTable$SDDSRVYR) # for A, B, C, D, E, F, G: # times recieve healthcare over last year
#table(MainTable$HUQ051, MainTable$SDDSRVYR) # for H: # times recieve healthcare over last year
MainTable$num_times_health_care <- NA ;
MainTable$num_times_health_care <-ifelse(MainTable$SDDSRVYR <= 7, MainTable$HUQ050, MainTable$HUQ051)
MainTable[which(MainTable$num_times_health_care > 10), 'num_times_health_care'] <- NA
MainTable[which(MainTable$num_times_health_care > 5), 'num_times_health_care'] <- 5


# number of overnoight hospital visits
#table(MainTable$HUQ070, MainTable$SDDSRVYR) # for A: Overnight hospital patient in last year
#table(MainTable$HUD070, MainTable$SDDSRVYR) # for B: Overnight Hospital patient in last year
#table(MainTable$HUQ071, MainTable$SDDSRVYR) # for C,D,E,F,G,H: Overnight Hospital patient in last year
#table(MainTable$HUD080, MainTable$SDDSRVYR) # for A,C,D,E,F,G,H: Times overnite hospital patient last year
#table(MainTable$HUQ080, MainTable$SDDSRVYR) # for B: Times overnite hospital patient last year

MainTable$overnight_in_hospital <- ifelse(MainTable$SDDSRVYR >= 3, MainTable$HUQ071, ifelse(MainTable$SDDSRVYR == 1, MainTable$HUQ070, MainTable$HUD070))
#table(MainTable$overnight_in_hospital, MainTable$SDDSRVYR)

MainTable$times_overnight_in_hospital <- ifelse(MainTable$overnight_in_hospital == 2, 
                                                0, 
                                                ifelse(MainTable$SDDSRVYR == 2, 
                                                       MainTable$HUQ080,
                                                       MainTable$HUD080
                                                       )
                                                )
MainTable$times_overnight_in_hospital[which(MainTable$times_overnight_in_hospital >= 7)] <- NA
#table(MainTable$times_overnight_in_hospital, MainTable$SDDSRVYR)


### now try it out - linear associations between overnight stay and num times getting health care
## go through each variable per survey and pump out association

adjustmentVariables <- c( "RIAGENDR", 'RIDAGEYR', 'I(RIDAGEYR^2)' ,'INDFMPIR', 'black', 'mexican' ,'other_eth', 'other_hispanic')
adjustmentVariables2 <- c( "RIAGENDR", 'RIDAGEYR' ,'INDFMPIR', 'black', 'mexican' ,'other_eth', 'other_hispanic')
variables <- finalVariableList$var
surveyVariables <- c('SDMVPSU', 'SDMVSTRA', 'WTMEC2YR')
indVars <- c('num_times_health_care', 'times_overnight_in_hospital') 

FILEOUT_NAME <- 'utilization_association.Rdata'
if(!is.na(SDDSRVYR_INPUT)) {
  cat(sprintf('doing analysis on survey number %.0f\n',SDDSRVYR_INPUT))
  MainTable <- subset(MainTable, SDDSRVYR == SDDSRVYR_INPUT)
  FILEOUT_NAME <- sprintf('utilization_association_%.0f.Rdata', SDDSRVYR_INPUT)
}

## iterate over survey, indVar                               
newFrame <- data.frame()
for(ii in 1:length(variables)) {
  for(jj in 1:length(indVars)) {
    formul <- sprintf('I(log10(%s+1e-10)) ~ %s + %s', variables[ii],indVars[jj], paste(adjustmentVariables, collapse='+')) 
    cat(sprintf('%0.0f,%s\n', ii, formul))
    toKeep <- c(indVars[jj],variables[ii], adjustmentVariables2, surveyVariables)
    tab <- MainTable[, toKeep]
    tab <- tab[complete.cases(tab), ]
    if(nrow(tab) < 100) {
      cat(sprintf('not enough data %.0f\n', nrow(tab)))
      next
    } 
    ret <- analyze_linear_glm(as.formula(formul),tab, 'WTMEC2YR') # need to do for each survey
    if(!is.null(ret)) {
      ret$dep_var <- variables[ii]
      ret$ind_var <- indVars[jj]
      ret$varnameModel <- rownames(ret)
      newFrame <- rbind(newFrame, ret)    
    }
  }
}

utilizationFrame <- newFrame
save(utilizationFrame, adjustmentVariables, SDDSRVYR_INPUT, FILEOUT_NAME, finalVariableList, file=FILEOUT_NAME)
