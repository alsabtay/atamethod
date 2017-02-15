atamethod <- function(X, sRange=c(1,400,1006,1007,1899,2129), parP="opt", parQ="opt", modelType="A", optimPar="optPQ", errType="sMAPE", fcast="single")
{

##atamethod(X, sRange="full", parP="opt", parQ="opt", modelType = "A", optimPar = "both", errType = 'sMAPE', fcast="single")
# X     	: nxm matrix and ts object
# sRange	: The selected series of M3
				# "full"       : all series are selected
				# 1:3	       : the 1st, 2nd and 3rd series are selected 
				# c(1,5,8,11)  : the 1st, 5th, 8th and 11th series are selected
				##    	 	 n : n <= length(X) - a single series is selected 
# parP 		: Level parameter 
				## 		 "opt" : p has all values from 1 to length(X)
				## 		   x:n : p has n values from x to n (n <= length(X))
				## c(1,5,8,11) : p has four values which are 1, 5, 8 and 11.
				##    	 	 n : n <= length(X) - p can take only one value
# parQ 		: Trend parameter
				## 		 "opt" : q has all values from 1 to length(X)
				## 		   x:n : q has n values from x to n (n <= length(X))
				## c(1,5,8,11) : q has four values which are 1, 5, 8 and 11.
				##    	 	 n : n <= length(X) - q can take only one value
# modelType	: 
			    ##  'A' for additive model
			    ##  'M' for multiplicative model
# optimPar	: selection of best model  
				##  'both'  	: Fits both ATA(p*,0) and ATA(p,1) (this option should be used when combines forecasts are needed or model selection is needed)
				##  'optPQ'  	: Fits ATA(p,q) where both parameters are optimized simultaneously (for p >= q)
				##  'pStar'     : Fits ATA(p*,0) where p* is the optimum value of p for q = 0
				##  'pOne'      : Fits ATA(p,1) where p is optimized for q = 1
				##  'pStarQ'    : Fits ATA(p*,q) where q is optimized for p = p*
# errType	: error type  
				##  'MAE'
				##  'MAPE'
				##  'sMAPE'				
				##  'MSE'
# fcast		: forecast model  
				##  'single' :  when the options ' optPQ ', 'pStar', 'pOne' or 'pStarQ' are chosen for optimPar, this option should be used to obtain forecasts from the chosen model  
				##  'select' :  a simple model selection of the two models ATA(p*,0) and ATA(p, 1) is carried out based on in-sample sMAPE (for this option optimPar should be 'both')
				##  'comb'   :  a simple average of the forecasts from the two models ATA(p*,0) and ATA(p, 1) is used as a forecast (for this option optimPar should be 'both')


## if (class(X)!="ts"){
##		argPar_1 <- tkmessageBox(title = "Data: Insufficient Parameters in ATA Method", message = "The data set must be time series object (ts) \rATA Method will terminate!", icon = "error", type = "ok")
##	    error_argPar_1 <- as.character(argPar_1)
##		return("The data set must be time series object (ts) \rATA Method was terminated!")
##	}

	WD <- getwd()
	
	if (!is.null(nrow(parQ))){
		Qlen <- nrow(parQ)
	}else {
		Qlen <- length(parQ)
	}
	if (!is.null(nrow(parP))){
		Plen <- nrow(parP)
	}else {
		Plen <- length(parP)
	}
	if (!is.null(nrow(sRange))){
		sRlen <- nrow(sRange)
	}else {
		sRlen <- length(sRange)
	}
	if (!is.null(nrow(modelType))){
		MTlen <- nrow(modelType)
	}else {
		MTlen <- length(modelType)
	}
	if (!is.null(nrow(optimPar))){
		optimlen <- nrow(optimPar)
	}else {
		optimlen <- length(optimPar)
	}
	if (!is.null(nrow(errType))){
		errlen <- nrow(errType)
	}else {
		errlen <- length(errType)
	}
	if (!is.null(nrow(fcast))){
		fclen <- nrow(fcast)
	}else {
		fclen <- length(fcast)
	}
	
	if (min(parP) == 0 | (class(parP) =="character" & Plen > 1)){
		argPar_2 <- tkmessageBox(title = "p: Insufficient Parameters in ATA Method", message = "p value must be greater than 0 and integer. \rATA Method will terminate!", icon = "error", type = "ok")
		error_argPar_2 <- as.character(argPar_2)
		return("p value must be greater than 0 and integer. ATA Method was terminated!")
	}
	if (class(sRange) == "character" & sRlen > 1){
		argPar_9 <- tkmessageBox(title = "sRange: Insufficient Parameters in ATA Method", message = "sRange value must be greater than 0 and integer. \rATA Method will terminate!", icon = "error", type = "ok")
		error_argPar_9 <- as.character(argPar_9)
		return("p value must be greater than 0 and integer. ATA Method was terminated!")
	}
	if ((class(parQ) == "character" & Qlen > 1)){
		argPar_3 <- tkmessageBox(title = "q: Insufficient Parameters in ATA Method", message = "q value must be greater than 0 and integer \rATA Method will terminate!", icon = "error", type = "ok")
		error_argPar_3 <- as.character(argPar_3)
		return("q value must be greater than 0 and integer. ATA Method was terminated!")
	}
	if (Qlen == 1 & Plen == 1){
		if (parQ > parP){
			argPar_8 <- tkmessageBox(title = "q: Insufficient Parameters in ATA Method", message = "p must be greater than q (p>q>=0)! \rATA Method will terminate!", icon = "error", type = "ok")
			error_argPar_8 <- as.character(argPar_8)
			return("p must be greater than q (p>q>=0)!. ATA Method was terminated!")
		}
	}
	if (Qlen > 1 & Plen >= 1){
		if (min(parQ) > max(parP)){
			argPar_8 <- tkmessageBox(title = "q: Insufficient Parameters in ATA Method", message = "p must be greater than q (p>q>=0)! \rATA Method will terminate!", icon = "error", type = "ok")
			error_argPar_8 <- as.character(argPar_8)
			return("p must be greater than q (p>q>=0)!. ATA Method was terminated!")
		}
	}
	if ((modelType != "A" & modelType != "M") | !is.character(modelType) | MTlen > 1){	
		argPar_4 <- tkmessageBox(title = "Model Type: Insufficient Parameters in ATA Method", message = "Model Type value must be string. \r    A for additive or (M for multiplicative) \rATA Method will terminate!", icon = "error", type = "ok")
		error_argPar_4 <- as.character(argPar_4)
		return("Model Type value must be string. A for additive or (M for multiplicative). ATA Method was terminated!")
	}
	if ((optimPar != "both" & optimPar != "optPQ" & optimPar != "pStar" & optimPar != "pOne" & optimPar != "pStarQ" ) | !is.character(optimPar) | optimlen > 1){		
		argPar_5 <- tkmessageBox(title = "Optimization Technique: Insufficient Parameters in ATA Method", message = "Model Type value must be string and it must get one value. \r  (both, optPQ or pStar or pOne or pStarQ)  \rATA Method will terminate!", icon = "error", type = "ok")
		error_argPar_5 <- as.character(argPar_5)
		return("Model Type value must be string. optPQ or both or pStar or pOne or pStarQ. ATA Method was terminated!")
	}	
	if ((errType != "MAE" & errType != "MAPE" & errType != "sMAPE" & errType != "MSE") | !is.character(errType) | errlen > 1){		
		argPar_6 <- tkmessageBox(title = "Error Selection Technique: Insufficient Parameters in ATA Method", message = "Error Type value must be string and it must get one value. \r  (MAE or MAPE or sMAPE or MSE) \rATA Method will terminate!", icon = "error", type = "ok")
		error_argPar_6 <- as.character(argPar_6)
		return("Model Type value must be string. MAE or MAPE or sMAPE or MSE. ATA Method was terminated!")
	}
	if ((fcast != "single" & fcast != "select" & fcast != "comb" ) | !is.character(fcast) | fclen > 1){		
		argPar_7 <- tkmessageBox(title = "Forecast Technique: Insufficient Parameters in ATA Method", message = "Forecast Type value must be string and it must get one value. \r  (single or select or comb) \rATA Method will terminate!", icon = "error", type = "ok")
		error_argPar_7 <- as.character(argPar_7)
		return("Model Type value must be string (single or select or comb). ATA Method was terminated!")
	}
	if (fcast == "comb"){
		fcast = "comb"
		parQ = 0:1
		parP = "opt"
		optimPar <- "both"
	}else if (fcast == "select"){
		parQ <- 0:1
		parP = "opt"
		optimPar <- "both"
	}else if (Qlen == 1 & Plen == 1 & parQ[1] != "opt" & parP[1] != "opt"){
		optimPar <- "optPQ"
		fcast = "single"
	}else if ((Qlen > 1 & Plen == 1) | (Qlen == 1 & Plen > 1)){
		optimPar <- "optPQ"
		fcast = "single"
	}else{
	}

	ptm <- proc.time()
	M3Results <<- matrix(NA,0,29)
	colnames(M3Results) <<- c("Series", "N", "Nh", "Seasonality", "Category", "Frequency","in_sample","out_sample","q", "p","F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14","F15","F16","F17","F18","F19")

		optAllResults <- atamethod.func(M3, srange = sRange ,optParP = parP, optParQ = parQ, mdlType = modelType, optoptimPar = optimPar, opterrType = errType, optfcast = fcast)
	
	WD <- getwd()
	subDir <- "ATA Method"
	ifelse(!dir.exists(file.path(WD, subDir)), dir.create(file.path(WD, subDir)), FALSE)	  
	txtWrite <- paste("write.csv(M3Results, file=\"",WD,"/",subDir,"/","M3Results.csv\")",sep = "")	
	eval(parse(text= txtWrite))
	ataResultOutputCSV(sRange,WD,subDir)
	gc()						
	proc.time() - ptm
}