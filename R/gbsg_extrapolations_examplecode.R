# Survival Extrapolation Incorporating General Population Mortality Using Excess 
# Hazard and Cure Models: a Tutorial

# Full Authors List
# Michael J Sweeting, Mark J Rutherford, Dan Jackson, Sangyu Lee, 
# Nicholas R Latimer, Robert Hettle, Paul C Lambert

# Code Authors
# Michael Sweeting, Mark Rutherford

# Contact
# michael.sweeting@astrazeneca.com

# Description of Project
# This tutorial demonstrates the use of excess hazard survival models for 
# survival extrapolation in Health Technology Assessment.

# Description of Code
# This code uses the freely available German Breast Cancer Study Group (GBCS) 
# dataset and fits a suite of seven parametric survival models that have been 
# recommended for consideration in health economic modelling. Excess hazard and 
# excess hazard cure models are fitted, and predictions of all-cause survival, 
# hazard, and restricted mean survival are obtained using the standsurv 
# functionality in the flexsurv R package.

# load required packages (and install first if not already installed)
list.of.packages <- c("flexsurv", "flexsurvcure", "dplyr", "tidyr", "ggplot2",
                      "gridExtra", "ggpubr", "devtools")
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
gbcsCS$sex <- factor("female")
## 2-level grade variable
gbcsCS$grade2 <- ifelse(gbcsCS$grade==3, "3", "1/2")
## Obtain attained age and attained calendar year in (whole) years
gbcsCS <- gbcsCS %>% mutate(attained.age.yr = floor(age + survtime/365.25),
                            attained.year = lubridate::year(diag + survtime))
head(gbcsCS)


# 2. lifetables ----
# We will use the US lifetables that come with the survival package
# First, let's reshape US lifetable to be a tidy data.frame and convert rates to
# per person-year as our survival analysis time scale will be in years
survexp.us.df <- as.data.frame.table(survexp.us, responseName = "exprate") %>%
  mutate(exprate = 365.25 * exprate)
survexp.us.df$age <- as.numeric(as.character(survexp.us.df$age))
survexp.us.df$year <- as.numeric(as.character(survexp.us.df$year))

# Now we merge in (left join) the US rates at the event times in the bc data
gbcsCS <- gbcsCS %>% left_join(survexp.us.df, by = c("attained.age.yr"="age", 
                                                     "attained.year"="year", 
                                                     "sex"="sex")) 
# Create a dataset containing Grade 1/2 only
gbcsCSLowGrade <- gbcsCS %>% filter(grade2=="1/2")


# 3. Model fitting to Good prognosis group ----
# Fit a standard parametric, excess hazard and mixture-cure model to the Good 
# prognosis group only.
models <- list()
# standard parametric
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  models[[dist]] <- flexsurvreg(Surv(survyrs, censdead)~1, 
                                data=gbcsCSLowGrade, dist=dist)
}
# excess hazard models
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  method <- ifelse(dist=="exp", "BFGS", "Nelder-Mead")
  models[[paste0(dist,".excesshazard")]] <- 
    flexsurvreg(Surv(survyrs, censdead)~1, 
                data=gbcsCSLowGrade, 
                dist=dist,
                bhazard=exprate, 
                method=method)
}

# excess hazard mixture-cure models
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  models[[paste0(dist,".excesshazardcure")]] <- 
    flexsurvcure(Surv(survyrs, censdead)~1, 
                 data=gbcsCSLowGrade, 
                 dist=dist,
                 bhazard=exprate,
                 method="Nelder-Mead")
}

# 4. AIC statisics ----
# For the excess hazard models the reported AIC is from a partial likelihood
# Therefore the AIC from an excess hazard model cannot be directly compared to
# an AIC from a standard parametric model
aic.tbl <- tibble()
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  aic.tbl <- aic.tbl %>% bind_rows(tibble(dist=dist, Method="param", 
                                          aic=models[[dist]]$AIC))
  # Excess hazard models
  dist.eh <- paste0(dist,".excesshazard")
  aic <- models[[dist.eh]]$AIC
  aic.tbl <- aic.tbl %>% bind_rows(tibble(dist=dist, Method="eh", aic=aic))
  
  # Excess hazard cure models
  dist.cure <- paste0(dist,".excesshazardcure")
  aic <- models[[dist.cure]]$AIC
  aic.tbl <- aic.tbl %>% bind_rows(tibble(dist=dist, Method="cure", aic=aic))
}
aic.tbl.wide <- aic.tbl %>% pivot_wider(names_from = Method, values_from = aic) 
aic.tbl.wide

# 5. Predicted and extrapolated all-cause survival and hazard ----
ss.surv <- tibble()
# We use the standsurv command in the flexsurv package to make the predictions
# Predicting all-cause survival over next 30-years in increments of 0.2 years
# Standard parametric model
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[dist]],
                           type="survival",
                           t=seq(0,30, by=0.2),
                           ci=F
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = FALSE) %>% 
    bind_cols(Cure = FALSE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
}
# Excess hazard model
# To get predictions of marginal all-cause survival, standsurv multiplies the 
# predicted relative survival function with the expected survival function for 
# each individual and then averages.
# We therefore must supply the expected ratetable for these calculations
# We also must scale from the ratetable time scale (days) to the regression 
# model time scale (years) using the scale.ratetable argument.
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[paste0(dist,".excesshazard")]],
                           type="survival",
                           t=seq(0,30, by=0.2),
                           ci=F,
                           rmap=list(sex = sex,
                                     year = diag,
                                     age = agedays
                           ),
                           ratetable = survexp.us,
                           scale.ratetable = 365.25,
                           newdata = gbcsCSLowGrade
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = TRUE) %>% 
    bind_cols(Cure = FALSE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
}
## Excess hazard cure model
for(dist in c("exp", "weibull", "gamma", "gompertz", "llogis", "lnorm",
              "gengamma")){
  ss.surv.new <- standsurv(models[[paste0(dist,".excesshazardcure")]],
                           type="survival",
                           t=seq(0,30, by=0.2),
                           ci=F,
                           rmap=list(sex = sex,
                                     year = diag,
                                     age = agedays
                           ),
                           ratetable = survexp.us,
                           scale.ratetable = 365.25,
                           newdata = gbcsCSLowGrade
  ) %>% 
    bind_cols(dist = dist) %>%
    bind_cols("Relative Survival" = TRUE) %>% 
    bind_cols(Cure = TRUE)
  ss.surv <- ss.surv %>% 
    bind_rows(ss.surv.new)
  ss.surv
}  
# standsurv can be re-run to predict
# all-cause hazard  using type="hazard"
# excess hazard (for excess hazard models) using type="excesshazard"
# RMST using type="rmst"
# Confidence intervals can be calculated from the delta method or bootstrapping
# see help("standsurv") for more details

# 6. Plot of extrapolated survival (faceted)----
####################################################################
ss.surv <- ss.surv %>% 
  mutate(Method = ifelse(`Relative Survival`==FALSE & Cure==FALSE,
                         "a) Standard parametric models",
                         ifelse(Cure==FALSE, 
                                "b) Excess hazard (no cure) models", "c) Excess hazard (cure) models")))
ss.surv$Method <- factor(ss.surv$Method, 
                         levels = c("a) Standard parametric models", 
                                    "b) Excess hazard (no cure) models", 
                                    "c) Excess hazard (cure) models"))
ggplot(ss.surv) + geom_line(aes(x=time,y=at1,color=dist)) + 
  facet_wrap(~ Method, ncol=1) + 
  theme_bw() +
  ylab("Survival") +
  xlab("Time (years)") + 
  theme(
    strip.background = element_rect(
      color="black", fill="white", size=0, linetype="solid"
    ),
    strip.text = element_text(hjust = 0, size=10),
    legend.position="bottom",
    legend.title=element_blank()
  )



# 7. Fitting a saturated model by treatment group----
################################################################################
models.sat <- list()
# We will demonstrate just using a standard parametric log-normal
# both the location and ancillary parameters are a function of grade2
# this is equivalent to fitting separate log-normal models to the two groups
models.sat[["lnorm"]] <- flexsurvreg(Surv(survyrs, censdead)~grade2, 
                                  anc=list(sdlog= ~grade2),
                                  data=gbcsCS, dist="lnorm")
# excess hazard model 
models.sat[["lnorm.excesshazard"]] <- 
  flexsurvreg(Surv(survyrs, censdead)~grade2, 
              anc=list(sdlog= ~grade2),
              data=gbcsCS, dist="lnorm",
              bhazard = exprate)
# excess hazard cure model
models.sat[["lnorm.excesshazardcure"]] <- 
  flexsurvcure(Surv(survyrs, censdead)~grade2, 
               anc=list(meanlog= ~grade2,
                        sdlog= ~grade2),
               data=gbcsCS, dist="lnorm",
               bhazard = exprate)

## 8. Predicted and extrapolated all-cause survival, hazard and hazard ratio----
################################################################################
ss.sat.surv <- tibble()
## All-cause survival
## Standard parametric (lognormal) model
## The attribute "standpred_at" gives predicted survival in long-format
ss.sat.surv.new <- attr(standsurv(models.sat[["lnorm"]],
                                       at=list(list(grade2="1/2"),
                                               list(grade2="3")),
                                       type="survival",
                                       t=seq(0,30, by=0.2),
                                       ci=F), "standpred_at") %>% 
    bind_cols(dist = "lnorm") %>%
    bind_cols("Relative Survival" = FALSE) %>% 
    bind_cols(Cure = FALSE)
ss.sat.surv <- ss.sat.surv %>% 
  bind_rows(ss.sat.surv.new)
## Excess hazard lognormal model
ss.sat.surv.new <- attr(standsurv(models.sat[["lnorm.excesshazard"]],
                                  at=list(list(grade2="1/2"),
                                          list(grade2="3")),
                                  type="survival",
                                  t=seq(0,30, by=0.2),
                                  ci=F,
                                  rmap=list(sex = sex,
                                            year = diag,
                                            age = agedays
                                  ),
                                  ratetable = survexp.us,
                                  scale.ratetable = 365.25,
                                  newdata = gbcsCS), "standpred_at") %>% 
    bind_cols(dist = "lnorm") %>%
    bind_cols("Relative Survival" = TRUE) %>% 
    bind_cols(Cure = FALSE)
ss.sat.surv <- ss.sat.surv %>% 
  bind_rows(ss.sat.surv.new)
## Excess hazard lognormal cure model
ss.sat.surv.new <- attr(standsurv(models.sat[["lnorm.excesshazardcure"]],
                                  at=list(list(grade2="1/2"),
                                          list(grade2="3")),
                                  type="survival",
                                  t=seq(0,30, by=0.2),
                                  ci=F,
                                  rmap=list(sex = sex,
                                            year = diag,
                                            age = agedays
                                  ),
                                  ratetable = survexp.us,
                                  scale.ratetable = 365.25,
                                  newdata = gbcsCS), "standpred_at") %>% 
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>% 
  bind_cols(Cure = TRUE)
ss.sat.surv <- ss.sat.surv %>% 
  bind_rows(ss.sat.surv.new)
  
## All-cause hazards and hazard ratio
ss.sat.haz <- ss.sat.HR <- tibble()
## Standard parametric lognormal
ss.sat.haz.run <- standsurv(models.sat[["lnorm"]],
                            at=list(list(grade2="1/2"),
                                    list(grade2="3")),
                            type="hazard",
                            contrast = "ratio",
                            atreference = 2,
                            t=seq(0,30, by=0.2),
                            ci=F)
ss.sat.haz.new <- attr(ss.sat.haz.run, "standpred_at")  %>% 
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = FALSE) %>% 
  bind_cols(Cure = FALSE)
ss.sat.HR.new <- attr(ss.sat.haz.run, "standpred_contrast")  %>% 
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = FALSE) %>% 
  bind_cols(Cure = FALSE)
ss.sat.haz <- ss.sat.haz %>% 
  bind_rows(ss.sat.haz.new)
ss.sat.HR <- ss.sat.HR %>%
  bind_rows(ss.sat.HR.new)

## Excess hazard lognormal model
ss.sat.haz.run <- standsurv(models.sat[["lnorm.excesshazard"]],
                            at=list(list(grade2="1/2"),
                                    list(grade2="3")),
                            type="hazard",
                            contrast="ratio",
                            atreference = 2,
                            t=seq(0,30, by=0.2),
                            ci=F,
                            rmap=list(sex = sex,
                                      year = diag,
                                      age = agedays
                            ),
                            ratetable = survexp.us,
                            scale.ratetable = 365.25,
                            newdata = gbcsCS)
ss.sat.haz.new <- attr(ss.sat.haz.run, "standpred_at")  %>% 
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>% 
  bind_cols(Cure = FALSE)
ss.sat.HR.new <- attr(ss.sat.haz.run, "standpred_contrast")  %>% 
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>% 
  bind_cols(Cure = FALSE)
ss.sat.haz <- ss.sat.haz %>% 
  bind_rows(ss.sat.haz.new)
ss.sat.HR <- ss.sat.HR %>%
  bind_rows(ss.sat.HR.new)

## Excess hazard lognormal cure model
ss.sat.haz.run <- standsurv(models.sat[["lnorm.excesshazardcure"]],
                            at=list(list(grade2="1/2"),
                                    list(grade2="3")),
                            type="hazard",
                            contrast = "ratio",
                            atreference = 2,
                            t=seq(0,30, by=0.2),
                            ci=F,
                            rmap=list(sex = sex,
                                      year = diag,
                                      age = agedays
                            ),
                            ratetable = survexp.us,
                            scale.ratetable = 365.25,
                            newdata = gbcsCS)
ss.sat.haz.new <- attr(ss.sat.haz.run, "standpred_at")  %>% 
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>% 
  bind_cols(Cure = TRUE)
ss.sat.HR.new <- attr(ss.sat.haz.run, "standpred_contrast")  %>% 
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>% 
  bind_cols(Cure = TRUE)
ss.sat.haz <- ss.sat.haz %>% 
  bind_rows(ss.sat.haz.new)
ss.sat.HR <- ss.sat.HR %>%
  bind_rows(ss.sat.HR.new)
  
# Restricted mean survival time at 30-years, with confidence intervals
# calculated using the delta method
ss.sat.rmst <- ss.sat.drmst <- tibble()
## Standard parametric lognormal
ss.sat.rmst.run <- standsurv(models.sat[["lnorm"]],
                            at=list(list(grade2="1/2"),
                                    list(grade2="3")),
                            type="rmst",
                            contrast = "diff",
                            atreference = 2,
                            t=30,
                            ci=T)
ss.sat.rmst.new <- attr(ss.sat.rmst.run, "standpred_at")  %>%
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = FALSE) %>%
  bind_cols(Cure = FALSE)
ss.sat.drmst.new <- attr(ss.sat.rmst.run, "standpred_contrast")  %>%
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = FALSE) %>%
  bind_cols(Cure = FALSE)
ss.sat.rmst <- ss.sat.rmst %>%
  bind_rows(ss.sat.rmst.new)
ss.sat.drmst <- ss.sat.drmst %>%
  bind_rows(ss.sat.drmst.new)

## Excess hazard lognormal model
ss.sat.rmst.run <- standsurv(models.sat[["lnorm.excesshazard"]],
                            at=list(list(grade2="1/2"),
                                    list(grade2="3")),
                            type="rmst",
                            contrast="diff",
                            atreference = 2,
                            t=30,
                            ci=T,
                            rmap=list(sex = sex,
                                      year = diag,
                                      age = agedays
                            ),
                            ratetable = survexp.us,
                            scale.ratetable = 365.25,
                            newdata = gbcsCS)
ss.sat.rmst.new <- attr(ss.sat.rmst.run, "standpred_at")  %>%
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>%
  bind_cols(Cure = FALSE)
ss.sat.drmst.new <- attr(ss.sat.rmst.run, "standpred_contrast")  %>%
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>%
  bind_cols(Cure = FALSE)
ss.sat.rmst <- ss.sat.rmst %>%
  bind_rows(ss.sat.rmst.new)
ss.sat.drmst <- ss.sat.drmst %>%
  bind_rows(ss.sat.drmst.new)

## Excess hazard lognormal cure model
ss.sat.rmst.run <- standsurv(models.sat[["lnorm.excesshazardcure"]],
                            at=list(list(grade2="1/2"),
                                    list(grade2="3")),
                            type="rmst",
                            contrast = "diff",
                            atreference = 2,
                            t=30,
                            ci=T,
                            rmap=list(sex = sex,
                                      year = diag,
                                      age = agedays
                            ),
                            ratetable = survexp.us,
                            scale.ratetable = 365.25,
                            newdata = gbcsCS)
ss.sat.rmst.new <- attr(ss.sat.rmst.run, "standpred_at")  %>%
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>%
  bind_cols(Cure = TRUE)
ss.sat.drmst.new <- attr(ss.sat.rmst.run, "standpred_contrast")  %>%
  bind_cols(dist = "lnorm") %>%
  bind_cols("Relative Survival" = TRUE) %>%
  bind_cols(Cure = TRUE)
ss.sat.rmst <- ss.sat.rmst %>%
  bind_rows(ss.sat.rmst.new)
ss.sat.drmst <- ss.sat.drmst %>%
  bind_rows(ss.sat.drmst.new)


## 9. Plot of extrapolated hazard ratios (faceted) ----
#######################################################
ss.sat.HR <- ss.sat.HR %>% 
  mutate(Method = ifelse(`Relative Survival`==FALSE & Cure==FALSE, 
                         "a) Standard parametric models",
                         ifelse(Cure==FALSE, 
                                "b) Excess hazard (no cure) models", 
                                "c) Excess hazard (cure) models")))
ss.sat.HR$Method <- factor(ss.sat.HR$Method, 
                           levels = c("a) Standard parametric models", 
                                      "b) Excess hazard (no cure) models", 
                                      "c) Excess hazard (cure) models"))
plot0 <- ggplot(ss.sat.HR) + 
  geom_line(aes(x=time,y=ratio, color=dist)) + 
  facet_wrap( ~ Method, ncol=1) + 
  ylab("Hazard Ratio (Grade 1/2 vs. Grade 3)") +
  xlab("Time (years)") + 
  scale_y_continuous(trans='log2',breaks=c(0.1, 0.3, 0.6, 1, 2, 4), 
                     lim=c(0.1,4)) +
  geom_hline(yintercept=1, alpha=0.4) +
  theme_bw() +
  theme(
    strip.background = element_rect(
      color="black", fill="white", size=0, linetype="solid"
    ),
    strip.text = element_text(hjust = 0, size=10),
    legend.position="bottom",
    legend.title=element_blank()
  )
plot0

# RMST Differences tables
tab1 <- ss.sat.drmst %>%
  filter(`Relative Survival`==FALSE & Cure==FALSE & !is.na(dist)) %>%
  arrange(difference) %>%
  mutate_if(is.numeric, ~sprintf("%.1f",.)) %>%
  mutate(`95% CI` = paste0("(",difference_lci,", ",difference_uci,")")) %>%
  select("dist","difference","95% CI") %>%
  rename("RMST \n difference"="difference")
table1 <- tableGrob(tab1, rows=NULL, theme = ttheme_minimal(base_size=7, padding = unit(c(2, 2), "mm")))
#
tab2 <- ss.sat.drmst %>%
  filter(`Relative Survival`==TRUE & Cure==FALSE & !is.na(dist)) %>%
  arrange(difference) %>%
  mutate_if(is.numeric, ~sprintf("%.1f",.)) %>%
  mutate(`95% CI` = paste0("(",difference_lci,", ",difference_uci,")")) %>%
  select("dist","difference","95% CI") %>%
  rename("RMST \n difference"="difference")
table2 <- tableGrob(tab2, rows=NULL, theme = ttheme_minimal(base_size=7, padding = unit(c(2, 2), "mm")))
#
tab3 <- ss.sat.drmst %>%
  filter(`Relative Survival`==TRUE & Cure==TRUE & !is.na(dist)) %>%
  arrange(difference) %>%
  mutate_if(is.numeric, ~sprintf("%.1f",.)) %>%
  mutate(`95% CI` = paste0("(",difference_lci,", ",difference_uci,")")) %>%
  select("dist","difference","95% CI") %>%
  rename("RMST \n difference"="difference")
table3 <- tableGrob(tab3, rows=NULL, theme = ttheme_minimal(base_size=7, padding = unit(c(2, 2), "mm")))
#
ggarrange(plot0,
          ggarrange(table1, table2, table3, ncol=1, nrow=4,heights=c(1,1,1,0.7)),
          ncol=2,
          widths = c(2,1))


