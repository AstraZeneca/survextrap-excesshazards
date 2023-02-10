# Survival Extrapolation Incorporating General Population Mortality Using Excess 
# Hazard and Cure Models: a Tutorial

# Full Authors List
# Michael J Sweeting, Mark J Rutherford, Dan Jackson, Sangyu Lee
# Nicholas R Latimer, Robert Hettle, Paul C Lambert

# Code Authors
# Michael Sweeting, Mark Rutherford

# Contact
# michael.sweeting@astrazeneca.com

# Description of Project
# This tutorial demonstrates the use of excess hazard survival models for 
# survival extrapolation in Health Technology Assessment.

# Description of Code
# This code saves Stata datasets of the freely available German Breast Cancer 
# Study Group (GBCS) dataset and US lifetables extracted from the survival
# package.

# load required packages (and install first if not already installed)
list.of.packages <- c("tidyverse", "haven", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
for(pkg in list.of.packages){
  library(pkg, character.only = T)
}
# We require data from the condSURV package, which is now archived on CRAN, 
# so we install a previous version
if(!("condSURV" %in% installed.packages()[,"Package"]))
  devtools::install_version("condSURV",version="2.0.2")

# 1. gbcsCS dataset ----
data(gbcsCS, package="condSURV")
## Create age at diagnosis in days - used later for matching to expected rates
gbcsCS$agedays <- floor(gbcsCS$age * 365.25)
## Survival time in years
gbcsCS$survyrs <- gbcsCS$survtime / 365.25
## Diagnosis as a date variable
gbcsCS$diag <- as.Date(as.character(gbcsCS$diagdateb), "%d-%m-%Y")
## Create sex (assume all are female)
gbcsCS$sex <- factor("female", levels=c("male","female"))
## 2-level grade variable
gbcsCS$grade2 <- ifelse(gbcsCS$grade==3, "3", "1/2")
## Obtain attained age and attained calendar year in (whole) years.
## Label these attained_age_yr and attined_year for Stata dataset
gbcsCS <- gbcsCS %>% mutate(attained_age_yr = floor(age + survtime/365.25),
                            attained_year = lubridate::year(diag + survtime))
head(gbcsCS)


# 2. lifetables ----
# We will use the US lifetables that come with the survival package
# First, let's reshape US lifetable to be a tidy data.frame and convert rates to
# per person-year as our survival analysis time scale will be in years
survexp.us.df <- as.data.frame.table(survexp.us, responseName = "exprate") %>%
  mutate(exprate = 365.25 * exprate)
# Use names attained_age_yr and attained_year for matching in Stata
survexp.us.df$attained_age_yr <- as.numeric(as.character(survexp.us.df$age))
survexp.us.df$attained_year <- as.numeric(as.character(survexp.us.df$year))
head(survexp.us.df)

## Write Stata dataset for gbcs data
write_dta(gbcsCS, "stata/gbcs.dta")

## Write Stata dataset for US lifetable
write_dta(survexp.us.df, "stata/survexp.us.dta")
