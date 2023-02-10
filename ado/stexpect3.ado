*! version 0.35 2019-08-15

/*
drop fuinteval?
checks variation in indweights if using by?
more labels
check for period analysis
*/

program define stexpect3, rclass
	version 15.0
    syntax using/                                               ///
	       [if] [in]                                            ///
		   ,                                                    ///
		   agediag(varname)                                     ///
		   datediag(varname)                                    ///
		   [                                                    ///  
		   by(varlist)                                          ///
		   every(real 0.5)                                      ///
		   indweights(varname)                                  ///
		   pmage(string)                                        ///
		   pmother(string)                                      ///
		   pmrate(string)                                       ///
		   pmyear(string)                                       ///
           pmmaxage(real 99)                                    ///
           pmmaxyear(real 100000)                               ///
		   maxt(real 0)                                         ///
		   verbose                                              ///
		   ]

    st_is 2 analysis
    marksample touse
    qui replace `touse' = 0  if _st==0 | _st==. 	

	
	if `maxt' == 0 {
		summ t, meanonly
		local maxt `r(maxt)'
	}
	

************
// Checks //	
************

    confirm var `agediag' `datediag'
    confirm var `pmother'
    foreach var in `agediag' `datediag' `pmother' {
       qui count if missing(`var') & `touse'
       if `r(N)' >0 {
           di as error "`var' has missing values"
       }
    }
	
    if "`by'"         != "" confirm var `by'    


	
*******************	
// popmort file ///
*******************
	if "`pmage'" == "" local pmage _age
	if "`pmyear'" == "" local pmyear _year
	if "`pmrate'" == "" local pmrate rate
	qui describe using "`using'", varlist short
	local popmortvars `r(varlist)'
	foreach var in `pmage' `pmyear' `pmother' `pmrate' {
		local varinpopmort:list posof "`var'" in popmortvars
		if !`varinpopmort' {
			di "`var' is not in popmort file"
			exit 198
		}
	}
	
    local popmortfile "`using'"
	
	mata: expsurv()

	
end

version 15.0
mata
void function expsurv() {
// Read in options
  //fuinterval = strtoreal(st_local("fuinterval"))
  verbose = st_local("verbose") != ""
  pmage = st_local("pmage")
  pmyear = st_local("pmyear")
  pmother = tokens(st_local("pmother"))
  pmrate= st_local("pmrate")
  pmmaxage = strtoreal(st_local("pmmaxage"))
  pmmaxyear = strtoreal(st_local("pmmaxyear"))
  pmvars = (pmage,pmyear,pmother,pmrate)
  maxt = strtoreal(st_local("maxt"))
  every = strtoreal(st_local("every"))
  hasby= st_local("by") != ""
  popmortfile = st_local("popmortfile")

// popmort file as a view
  if(verbose) printf("Reading in popmort file\n")
  stata("preserve")
  stata("use " + popmortfile + ", clear")
  popmort = st_data(.,pmvars,.)
  stata("restore")

// store popmort file in array
  Npmvars = 2 :+ cols(pmother)
  pm = asarray_create("real",Npmvars)
  Nrows_pm = rows(popmort)
  ratevarcol = Npmvars + 1
  for(i=1;i<=Nrows_pm;i++) {
	asarray(pm,popmort[i,(1..Npmvars)],popmort[i,ratevarcol])
  }

// read in relevant data
  touse     = st_local("touse")
 
  datediag_varname   = st_local("datediag")
  datediag_varlabel = st_varlabel(datediag_varname)
  datediag  = st_data(.,datediag_varname,touse)

  agediag_varname = st_local("agediag")
  agediag_varlabel = st_varlabel(agediag_varname)
  agediag   = st_data(.,agediag_varname,touse)
  Nobs = rows(agediag)
  
  
  // by variables
  if(hasby) {
	byvars = st_local("by")
	by = st_data(.,byvars,touse)
	bylevels = uniqrows(by)
	Nbylevels = rows(bylevels)
	Nbyvars = cols(by)
  }
  else {
	byvars = J(1,0,"")
	by = J(Nobs,0,1)
	bylevels = 1
	Nbylevels = 1
	Nbyvars = 1
  }
  

  Npmother = cols(pmother)
  pmothervars = J(Nobs,Npmother,.)
  for(i=1;i<=cols(pmother);i++) {
	pmothervars[,i] = st_data(.,pmother[1,i],touse)
  }

//  unique_t = uniqrows(select(t,d)) 
//  Nunique_t = rows(unique_t)
  dots = J(Nobs,1,.)
  
  // col1: attained age. col2: attained year
  data = (dots,dots,pmothervars)  

// create matrix of yearly rates
  fumax_int = ceil(maxt)
  t_rates_start = range(0,fumax_int:-every,every)
  t_rates_stop  = range(every,fumax_int,every)
  N_t_rates = rows(t_rates_start)

// loop over intervals to calulate expected survival
  printf("Calculating expected survival\n")
  exprates = J(Nobs,(N_t_rates),.)
  pmmaxyear_vec = J(Nobs,1,pmmaxyear)
  pmmaxage_vec =  J(Nobs,1,pmmaxage)
  for(j=1;j<=N_t_rates;j++) {
	//midpoint of interval
  	data[,1] = rowmin((floor(agediag  :+ t_rates_start[j]:+every:/2),pmmaxage_vec))
  	data[,2] = rowmin((year(datediag  :+ (t_rates_start[j]:+every:/2)*365.24),pmmaxyear_vec))
 	for(i=1;i<=Nobs;i++) {
      exprates[i,j] = asarray(pm,(data[i,.]))
  	}
  }
  
  hazard_contrib = exprates:*every
  expsurv_i = J(Nobs,N_t_rates,.)
  for(i=1;i<=Nobs;i++) {
	expsurv_i[i,] = exp(-runningsum(hazard_contrib[i,]))
  }
  meansurv = mean(expsurv_i)
  
  meanhaz = J(1,N_t_rates,.)
  for(j=1;j<=N_t_rates;j++) {
	meanhaz[1,j] =  mean(exprates[,j]:*expsurv_i[,j]):/meansurv[1,j]
  }
  newvars = ("t_exp","expsurv","exphaz")
  (void) st_addvar("double", newvars)
  
  
  st_store((1..(rows(t_rates_stop)+1))',newvars,((0\t_rates_stop),(1\meansurv'),(.\meanhaz')))
 
}

end










