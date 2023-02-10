/*
Survival Extrapolation Incorporating General Population Mortality Using Excess 
Hazard and Cure Models: a Tutorial

Full Authors List
Michael J Sweeting, Mark J Rutherford, Dan Jackson, Sangyu Lee,
Nicholas R Latimer, Robert Hettle, Paul C Lambert

Code Authors
Michael Sweeting, Mark Rutherford

Contact
michael.sweeting@astrazeneca.com

Description of Project
This tutorial demonstrates the use of excess hazard survival models for 
survival extrapolation in Health Technology Assessment.

Description of Code
This code uses the freely available German Breast Cancer Study Group (GBCS) 
dataset and fits three parametric models (Weibull, Log-Normal, and Log-Logistic)
that have been recommended for consideration in health economic modelling. 
Excess hazard and excess hazard cure models are fitted using stpm2, and 
predictions of all-cause survival, hazard, and restricted mean survival 
are obtained using the standsurv package
*/

* Run the R script gbsg_stata_data.R first to get Stata data files for the GBCS
* data and expected survival from a US population

* To get all-cause survival and all-cause hazard from cure models we must use 
* the stexpect3 function to get expected survival. Place this ado file in your 
* personal adopath

* We use the stpm2 flexible parametric survival package to fit Weibull, 
* Log-logistic and Log-normal models together with their excess hazard 
* counterparts

* We use standsurv (on the ssc) to calculate marginal all-cause survival and 
* hazards
* ssc install stpm2
* ssc install standsurv /// installing standsurv from ssc


**SET UP GRAPH SCHEME**
grstyle init
grstyle set imesh, horizontal compact
grstyle set color Dark2,  n(8)
grstyle set color mono,  n(5) : histogram histogram_line
grstyle set intensity 80,: histogram

** Use the gbcs data set (saved using haven package in R)
use gbcs, clear

* Merge in US lifetable, keep only matches from master 
* (survexp.us.dta saved from haven package in R)
merge m:1 attained_age_yr attained_year sex using survexp.us.dta, keep(3)

* Generate grade 1/2 and grade 3 groups
gen grade_grp=grade==3
label define grade 0 "Grade 1/2" 1 "Grade 3"
label values grade_grp grade 

* stset the data
stset survtime censdead, scale(365.25)

* Kaplan-Meier plot by groups
sts graph, by(grade_grp) risktable(,order(1 "Grade 1/2" 2 "Grade 3") /// 
color("scheme p1")) risktable(,color("scheme p2") group(#2)) legend(off) ///
xtitle("Time since surgery (years)") ytitle("Survival Proportion") ///
name(KM_by_group, replace) title("")


* create some times for prediction (from 0 to 30 years, 601 data points)
range t30 0 30 601

* Standard parametric Weibull
* Log-likelihood differs from that reported by flexsurvreg in R because
* in flexsurvreg the sum of the log uncensored survival times is added to the 
* log-likelihood in survival models, to remove dependency on the time scale.
* The parameter estimates in streg below agree with flexsurvreg(dist="weibullPH")
* _cons == scale, p == shape
/*Call:
flexsurvreg(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    dist = "weibullPH")

Estimates: 
       est      L95%     U95%     se     
shape  1.69405  1.44261  1.98931  0.13887
scale  0.02035  0.01289  0.03211  0.00474

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -419.611, df = 2
AIC = 843.2221
*/
streg if grade_grp==0, dist(weibull)
* Alternatively, we can use a flexible parametric model on the log cumulative 
* hazard scale with 1 degree of freedom. This is equivalent to fitting a Weibull 
* PH model and gives the same log likeihood as above. The parameter estimates 
* from stpm2 correspond to the following in flexsurvreg:
* _cons = log(scale), _rcs1 = shape
stpm2 if grade_grp==0, scale(hazard) df(1) noorthog
* Using standsurv to get predicted all-cause survival
standsurv if grade_grp==0, at1(.) ///
	surv ///
	atvar(surv_weibull) /// new variable for all-cause survival prediction
	timevar(t30) 
* Using standsurv to get predicted all-cause hazard
standsurv if grade_grp==0, at1(.) ///
	hazard ///
	atvar(haz_weibull) /// new variable for all-cause hazard prediction
	timevar(t30) 

* Excess hazard Weibull model
* The parameter estimates from stpm2 correspond to the following in flexsurvreg:
* _cons = log(scale), _rcs1 = shape
/* Call:
flexsurvreg(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    bhazard = exprate, dist = "weibullPH")

Estimates: 
       est      L95%     U95%     se     
shape  1.75335  1.46284  2.10155  0.16205
scale  0.01615  0.00943  0.02765  0.00443

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -403.9079, df = 2
AIC = 811.8159
*/
stpm2 if grade_grp==0, scale(hazard) df(1) bhaz(exprate) noorthog
* Using standsurv to get predicted all-cause survival. We need to feed the 
* function the expected lifetable to get back to all-cause survival
standsurv if grade_grp==0, at1(.) ///
	surv ///
	atvar(surv_weibull_eh) /// new variable for all-cause survival prediction
	timevar(t30) ///
	 expsurv(using("survexp.us.dta") ///  Popmort file
                        agediag(age)        ///  Age at diagnosis in years in the dataset
                        datediag(diag) 	    ///  Date of diagnosis in the dataset
			pmage(attained_age_yr) 	    /// Age variable in population mortality file
			pmyear(attained_year) /// Calendar year variable in the population mortality file
                        pmother(sex)       ///  Other variables included in the popmort file
                        pmrate(exprate)       ///  Rate variable in the popmort file  
			pmmaxage(109)		     ///  Maximum age in the popmort file
			pmmaxyear(2014)       ///  Maximum year in the popmort 
                 )  
* Using standsurv to get predicted all-cause hazard. We need to feed the 
* function the expected lifetable to get back to all-cause hazard
standsurv if grade_grp==0, at1(.) ///
	hazard ///
	atvar(haz_weibull_eh) /// new variable for all-cause hazard prediction
	timevar(t30) ///
	 expsurv(using("survexp.us.dta") ///  Popmort file
                        agediag(age)        ///  Age at diagnosis in years in the dataset
                        datediag(diag) 	    ///  Date of diagnosis in the dataset
			pmage(attained_age_yr) 	    /// Age variable in population mortality file
			pmyear(attained_year) /// Calendar year variable in the population mortality file
                        pmother(sex)       ///  Other variables included in the popmort file
                        pmrate(exprate)       ///  Rate variable in the popmort file  
			pmmaxage(109)		     ///  Maximum age in the popmort file
			pmmaxyear(2014)       ///  Maximum year in the popmort 
                 )  

** Excess hazard cure Weibull model
* The parameter estimates from strsmix are close but not exact to the following 
* in flexsurvcure with dist = "weibullPH":
* pi = theta; ln_lambda = log(scale); ln_gamma = log(shape)
* The log-likelihoods however are identical
/*Call:
flexsurvcure(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    bhazard = exprate, dist = "weibullPH")

Estimates: 
       est     L95%    U95%    se    
theta  0.6564  0.5260  0.7669      NA
shape  2.2114  1.7760  2.7536  0.2474
scale  0.0356  0.0194  0.0656  0.0111

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -401.2316, df = 3
AIC = 808.4632
*/
strsmix if grade_grp==0, link(logistic) distribution(weibull) bhazard(exprate)
predict relsurv_weibull_ehcure, surv timevar(t30) 
predict ehaz_weibull_ehcure, hazard timevar(t30) 
* Get the expected survival and multiply by relative survival 
** NOTE; THIS WILL NOT GIVE THE CORRECT ANSWER IF RELATIVE SURVIVAL VARIES BY AGE (OR OTHER COVARIATES)
stexpect3 if grade_grp==0 using survexp.us.dta, ///
	agediag(age) /// Age at diagnosis in years in the dataset
	datediag(diag) /// Date of diagnosis in the dataset
	pmage(attained_age_yr) /// Age variable in population mortality file
	pmyear(attained_year) /// Calendar year variable in the population mortality file
	pmother(sex) /// Other variables included in the popmort file 
	pmrate(exprate)       ///  Rate variable in the popmort file  
	pmmaxage(109)		     ///  Maximum age in the popmort file
	pmmaxyear(2014)       ///  Maximum year in the popmort 
	every(0.05) /// calculate every 0.05 years. THIS NEEDS TO BE THE SAME EVERY AS IN t30
	maxt(30) // maximum follow-up time in years. THIS NEEDS TO BE THE SAME MAX AS IN t30
gen surv_weibull_ehcure = relsurv_weibull_ehcure * expsurv // all-cause survival
gen haz_weibull_ehcure = ehaz_weibull_ehcure + exphaz // all-cause hazard. T
	

		 
* Standard parametric Log-Logistic
* The parameter estimates in streg below agree with flexsurvreg(dist="llogis") 
* after transformation.
* exp(_cons) == scale, 1/gamma = shape
/*Call:
flexsurvreg(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    dist = "llogis")

Estimates: 
       est    L95%   U95%   se   
shape  1.854  1.583  2.173  0.150
scale  8.504  7.250  9.976  0.693

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -418.1352, df = 2
AIC = 840.2705
*/
streg if grade_grp==0, dist(llogis)
* Alternatively, we can use a flexible parametric model on the log odds scale 
* with 1 degree of freedom. This is equivalent to fitting a Log-logistic AFT 
* model and gives the same log likeihood as above. The parameter estimates from 
* stpm2 correspond to the following:
* exp(-_cons / _rcs1)= scale,  _rcs1 = shape 
stpm2 if grade_grp==0, scale(odds) df(1) noorthog
* Using standsurv to get predicted all-cause survival.
standsurv if grade_grp==0, at1(.) ///
	surv ///
	atvar(surv_llogis) /// new variable for all-cause survival prediction
	timevar(t30) 
* Using standsurv to get predicted all-cause hazard
standsurv if grade_grp==0, at1(.) ///
	hazard ///
	atvar(haz_llogis) /// new variable for all-cause survival prediction
	timevar(t30) 
	
* Excess hazard Log-Logistic model
* The parameter estimates from stpm2 (without orthogonalisation) correspond 
* to the following:
* exp(-_cons / _rcs1)= scale,  _rcs1 = shape 
/*
Call:
flexsurvreg(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    bhazard = exprate, dist = "llogis")

Estimates: 
       est     L95%    U95%    se    
shape   1.902   1.590   2.275   0.174
scale   9.101   7.588  10.915   0.844

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -402.7109, df = 2
AIC = 809.4217
*/
stpm2 if grade_grp==0, scale(odds) df(1) bhaz(exprate) noorthog
* Using standsurv to get predicted all-cause survival. We need to feed the 
* function the expected lifetable to get back to all-cause survival
standsurv if grade_grp==0, at1(.) ///
	surv ///
	atvar(surv_llogis_eh) /// new variable for all-cause survival prediction
	timevar(t30) ///
	 expsurv(using("survexp.us.dta") ///  Popmort file
                        agediag(age)        ///  Age at diagnosis in years in the dataset
                        datediag(diag) 	    ///  Date of diagnosis in the dataset
			pmage(attained_age_yr) 	    /// Age variable in population mortality file
			pmyear(attained_year) /// Calendar year variable in the population mortality file
                        pmother(sex)       ///  Other variables included in the popmort file
                        pmrate(exprate)       ///  Rate variable in the popmort file  
			pmmaxage(109)		     ///  Maximum age in the popmort file
			pmmaxyear(2014)       ///  Maximum year in the popmort 
                 )  
* Using standsurv to get predicted all-cause hazard. We need to feed the 
* function the expected lifetable to get back to all-cause hazard
standsurv if grade_grp==0, at1(.) ///
	hazard ///
	atvar(haz_llogis_eh) /// new variable for all-cause hazard prediction
	timevar(t30) ///
	 expsurv(using("survexp.us.dta") ///  Popmort file
                        agediag(age)        ///  Age at diagnosis in years in the dataset
                        datediag(diag) 	    ///  Date of diagnosis in the dataset
			pmage(attained_age_yr) 	    /// Age variable in population mortality file
			pmyear(attained_year) /// Calendar year variable in the population mortality file
                        pmother(sex)       ///  Other variables included in the popmort file
                        pmrate(exprate)       ///  Rate variable in the popmort file  
			pmmaxage(109)		     ///  Maximum age in the popmort file
			pmmaxyear(2014)       ///  Maximum year in the popmort 
	)  

** Excess hazard cure Log-logistic model
** Cannot implement using strsmix	

* Standard parametric Log-Normal
* The parameter estimates in streg below agree with flexsurvreg(dist="lnorm")
* _cons == meanlog, sigma == sdlog
/*
Call:
flexsurvreg(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    dist = "lnorm")

Estimates: 
         est     L95%    U95%    se    
meanlog  2.2325  2.0460  2.4190  0.0952
sdlog    1.0173  0.8796  1.1765  0.0755

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -415.7245, df = 2
AIC = 835.449
*/
streg if grade_grp==0, dist(lnorm)
* Alternatively, we can use a flexible parametric model on the probit scale 
* with 1 degree of freedom. This is equivalent to fitting a Log-normal AFT model
* and gives the same log likeihood as above. The parameter estimates from stpm2 
* correspond to the following:
*  -_cons/_rcs1 = meanlog, 1/_rcs1 = sdlog
stpm2 if grade_grp==0, scale(normal) df(1) noorthog
* Using standsurv to get predicted all-cause survival.
standsurv if grade_grp==0, at1(.) ///
	surv ///
	atvar(surv_lnorm) /// new variable for all-cause survival prediction
	timevar(t30) 
* Using standsurv to get predicted all-cause hazard
standsurv if grade_grp==0, at1(.) ///
	hazard ///
	atvar(haz_lnorm) /// new variable for all-cause survival prediction
	timevar(t30) 

* Excess hazard Log-Normal model
* The parameter estimates from stpm2 (without orthogonalisation) correspond to 
* the following:
*  -_cons/_rcs1 = meanlog, 1/_rcs1 = sdlog
/*
Call:
flexsurvreg(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    bhazard = exprate, dist = "lnorm")

Estimates: 
         est    L95%   U95%   se   
meanlog  2.285  2.077  2.494  0.106
sdlog    0.975  0.825  1.152  0.083

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -399.3806, df = 2
AIC = 802.7612
*/
stpm2 if grade_grp==0, scale(normal) df(1) bhaz(exprate) noorthog
* Using standsurv to get predicted all-cause survival. We need to feed the 
* function the expected lifetable to get back to all-cause survival
standsurv if grade_grp==0, at1(.) ///
	surv ///
	atvar(surv_lnorm_eh) /// new variable for all-cause survival prediction
	timevar(t30) ///
	 expsurv(using("survexp.us.dta") ///  Popmort file
                        agediag(age)        ///  Age at diagnosis in years in the dataset
                        datediag(diag) 	    ///  Date of diagnosis in the dataset
			pmage(attained_age_yr) 	    /// Age variable in population mortality file
			pmyear(attained_year) /// Calendar year variable in the population mortality file
                        pmother(sex)       ///  Other variables included in the popmort file
                        pmrate(exprate)       ///  Rate variable in the popmort file  
			pmmaxage(109)		     ///  Maximum age in the popmort file
			pmmaxyear(2014)       ///  Maximum year in the popmort 
                 )  
* Using standsurv to get predicted all-cause hazard. We need to feed the 
* function the expected lifetable to get back to all-cause hazard
standsurv if grade_grp==0, at1(.) ///
	hazard ///
	atvar(haz_lnorm_eh) /// new variable for all-cause hazard prediction
	timevar(t30) ///
	 expsurv(using("survexp.us.dta") ///  Popmort file
                        agediag(age)        ///  Age at diagnosis in years in the dataset
                        datediag(diag) 	    ///  Date of diagnosis in the dataset
			pmage(attained_age_yr) 	    /// Age variable in population mortality file
			pmyear(attained_year) /// Calendar year variable in the population mortality file
                        pmother(sex)       ///  Other variables included in the popmort file
                        pmrate(exprate)       ///  Rate variable in the popmort file  
			pmmaxage(109)		     ///  Maximum age in the popmort file
			pmmaxyear(2014)       ///  Maximum year in the popmort 
	)  


** Excess hazard cure lognormal model
* The parameter estimates from strsmix are close to those in flexsurvreg.
* pi = logit(theta); mu = meanlog; ln_sigma = log(sdlog)
* The log-likelihoods however are identical
/*Call:
flexsurvcure(formula = Surv(survyrs, censdead) ~ 1, data = gbcsCSLowGrade, 
    bhazard = exprate, dist = "lnorm")

Estimates: 
         est    L95%   U95%   se   
theta    0.593  0.385  0.773     NA
meanlog  1.443  1.006  1.881  0.223
sdlog    0.666  0.485  0.913  0.107

N = 525,  Events: 113,  Censored: 412
Total time at risk: 1975.181
Log-likelihood = -397.512, df = 3
AIC = 801.0241
*/
strsmix if grade_grp==0, link(logistic) distribution(lognormal) bhazard(exprate)
predict relsurv_lnorm_ehcure, surv timevar(t30) 
predict ehaz_lnorm_ehcure, hazard timevar(t30) 
* Get the expected survival and multiply by relative survival 
** NOTE; THIS WILL NOT GIVE THE CORRECT ANSWER IF RELATIVE SURVIVAL VARIES BY 
* AGE (OR OTHER COVARIATES)
gen surv_lnorm_ehcure = relsurv_lnorm_ehcure * expsurv // all-cause survival
gen haz_lnorm_ehcure = ehaz_lnorm_ehcure + exphaz // all-cause hazard.

		 
**********************************************************************************	
* Plot extrapolated survival
* THIS GIVES THE SAME PREDICTIONS AS IN R EXCEPT FOR THE EXCESS HAZARD LOG-NORMAL
line surv_weibull surv_llogis surv_lnorm ///
	surv_weibull_eh surv_llogis_eh surv_lnorm_eh ///
	surv_weibull_ehcure surv_lnorm_ehcure ///
	t30, sort ///
	lcolor(blue orange green blue orange green blue green) ///
	lpattern(solid solid solid dash dash dash dot dot) ///
	legend(order(1 "Weibull (standard parametric)" ///
	2 "Log-logistic (standard parametric)" ///
	3 "Log-Normal (standard parametric)" ///
	4 "Weibull (excess hazard model)" ///
	5 "Log-logistic (excess hazard model)" ///
	6 "Log-Normal (excess hazard model)" ///
	7 "Weibull (excess hazard cure model)" ///
	8 "Log-Normal (excess hazard cure model)")) ///
	name(surv_extrap,replace) xtitle("Time (years)") ///
	ylabel(,format(%3.2f)) ytitle("Survival")
	
* Plot extrapolated hazard
* THIS GIVES THE SAME PREDICTIONS AS IN R EXCEPT FOR THE EXCESS HAZARD LOG-NORMAL
preserve 
replace haz_weibull = . if haz_weibull >0.15
replace haz_llogis = . if haz_llogis >0.15
replace haz_lnorm = . if haz_lnorm >0.15
replace haz_weibull_eh = . if haz_weibull_eh>0.15
replace haz_llogis_eh = . if haz_llogis_eh >0.15
replace haz_lnorm_eh = . if haz_lnorm_eh > 0.15
replace haz_weibull_ehcure = . if haz_weibull_ehcure>0.15
replace haz_lnorm_ehcure = . if haz_lnorm_ehcure>0.15
line haz_weibull haz_llogis haz_lnorm ///
	haz_weibull_eh haz_llogis_eh haz_lnorm_eh ///
	haz_weibull_ehcure haz_lnorm_ehcure ///
	t30, sort ///
	lcolor(blue orange green blue orange green blue green) ///
	lpattern(solid solid solid dash dash dash dot dot) ///
	legend(order(1 "Weibull (standard parametric)" ///
	2 "Log-logistic (standard parametric)" ///
	3 "Log-Normal (standard parametric)" ///
	4 "Weibull (excess hazard model)" ///
	5 "Log-logistic (excess hazard model)" ///
	6 "Log-Normal (excess hazard model)" ///
	7 "Weibull (excess hazard cure model)" ///
	8 "Log-Normal (excess hazard cure model)")) ///
	name(haz_extrap,replace) xtitle("Time (years)") ///
	ylabel(0(0.05)0.15,format(%3.2f)) ytitle("Hazard (per person year)")
restore

