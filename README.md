# survextrap_excesshazards
![Maturity level-Prototype](https://img.shields.io/badge/Maturity%20Level-Prototype-red) 

This repository contains example R and Stata scripts to estimate parametric excess hazard and excess hazard cure survival models, and to derive predictions of all-cause survival, hazard and restricted mean survival. 
The repository contains two main scripts: 
 - ``R/gbsg_extrapolations_examplecode.R`` - Example R script to fit parametric survival models and obtain predictions.
 - ``stata/gbsg_extrapolations_examplecode.do`` - Example Stata do file to fit parametric survival models and obtain predictions;

Requirements
------------
Analyses were performed using R 4.1.0 and Stata/MP 17.0. 

The project requires the following R packages. The version numbers indicate the version of the packages 
that were used in the analysis. 
```
  ggpubr==0.4.0
  gridExtra==2.3
  condSURV==2.0.2
  ggplot2==3.3.6
  tidyr==1.2.0
  dplyr==1.0.9
  flexsurvcure==1.2.0 
  flexsurv==2.2
  survival==3.2-13 
```

Folder Contents
----------------
This folder contains the data files that was used in the analysis. The file descriptions are listed below:

```
|---- README.md :  readme file for the data folder
|
|---- R : This folder contains the R scripts 
|------- R/gbsg_extrapolations_examplecode.R : The demonstration R script for fitting parametric survival models and obtaining predictions
|------- R/gbsg_stata_data.R : An R script for generating Stata data files for the GBCS data and US population lifetable data
|
|---- stata:  This folder contains the Stata script.
|------- stata/gbsg_extrapolations_examplecode.do : The demonstration Stata do file for fitting parametric survival models and obtaining predictions
|---- ado: This folder contains additional ado files for use with the Stata script
|------- stexpect.ado : Calculates expected survival for a cohort of individuals using age at diagnosis, calendar year, 
sex and other matching variables in the popmort file
```

Contact
--------
michael.sweeting@astrazeneca.com
