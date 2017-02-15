atamethod.func <- function(X, srange, optParP, optParQ, mdlType, optoptimPar, opterrType, optfcast)
{

# X     	: nxm matrix and ts object
# srange	: The selected series of M3
				# "full"       : all series are selected
				# 1:3	       : the 1st, 2nd and 3rd series are selected 
				# c(1,5,8,11)  : the 1st, 5th, 8th and 11th series are selected
				##    	 	 n : n <= length(X) - a single series is selected 
# optparP 		: Level parameter 
				## 		 "opt" : p has all values from 1 to length(X)
				## 		   x:n : p has n values from x to n (n <= length(X))
				## c(1,5,8,11) : p has four values which are 1, 5, 8 and 11.
				##    	 	 n : n <= length(X) - p can take only one value
# optparQ 		: Trend parameter
				## 		 "opt" : q has all values from 1 to length(X)
				## 		   x:n : q has n values from x to n (n <= length(X))
				## c(1,5,8,11) : q has four values which are 1, 5, 8 and 11.
				##    	 	 n : n <= length(X) - q can take only one value
# mdlType	: 
			    ##  'A' for additive model
			    ##  'M' for multiplicative model
# optoptimPar	: selection of best model  
				##  'both'  	: Fits both ATA(p*,0) and ATA(p,1) (this option should be used when combines forecasts are needed or model selection is needed)
				##  'optPQ'  	: Fits ATA(p,q) where both parameters are optimized simultaneously (for p >= q)
				##  'pStar'     : Fits ATA(p*,0) where p* is the optimum value of p for q = 0
				##  'pOne'      : Fits ATA(p,1) where p is optimized for q = 1
				##  'pStarQ'    : Fits ATA(p*,q) where q is optimized for p = p*
# opterrType	: error type  
				##  'MAE'
				##  'MAPE'
				##  'sMAPE'				
				##  'MSE'
# optfcast		: forecast model  
				##  'single' :  when the options ' optPQ ', 'pStar', 'pOne' or 'pStarQ' are chosen for optimPar, this option should be used to obtain forecasts from the chosen model  
				##  'select' :  a simple model selection of the two models ATA(p*,0) and ATA(p, 1) is carried out based on in-sample sMAPE (for this option optimPar should be 'both')
				##  'comb'   :  a simple average of the forecasts from the two models ATA(p*,0) and ATA(p, 1) is used as a forecast (for this option optimPar should be 'both')

	if (srange[1] == "full"){
		colSeries <- length(M3)
	}else{
		colsRange <- as.data.frame(srange)
		colSeries <- nrow(colsRange)
	}
	maxfHorizon <- 18
	ssindex <- read.csv("M3_Season.csv")
	for (ll in 1:colSeries){
		if (srange[1]=="full"){
			k <- ll
		}else{
			k <- as.numeric(colsRange[ll,])
		}
		Xh <- M3[[k]]$xx
		lenXh <- length(Xh)
		VarSeriesName <- M3[[k]]$sn 
		obsvNum <- M3[[k]]$n
		fHorizon <- M3[[k]]$h
		tscategory <- M3[[k]]$type
		obsPeriod <- M3[[k]]$period
		is.season <- ifelse (nrow(ssindex[ssindex["Series"]==M3[[k]]$sn,])>0,1,0)
		if (mdlType=="A" & is.season ==1){
			desX <- decompose(M3[[k]]$x, type = c("multiplicative"))
			# desX <- decompose(M3[[k]]$x, type = c("additive"))
			X <- seasadj(desX)
			lenX <- length(X)
			freqXh <- cycle(M3[[k]]$xx)
			freqSeas <- M3[[k]]$period
			SeasIndex <<- desX$figure
			SeasActual <<- desX$seasonal
		}else if (mdlType=="M" & is.season ==1){
			desX <- decompose(M3[[k]]$x, type = c("multiplicative"))
			X <- seasadj(desX)
			lenX <- length(X)
			freqXh <- cycle(M3[[k]]$xx)
			freqSeas <- M3[[k]]$period
			SeasIndex <<- desX$figure
			SeasActual <<- desX$seasonal
		}else{
			X <- M3[[k]]$x
			lenX <- length(X)
			freqXh <- cycle(M3[[k]]$xx)
			freqSeas <- M3[[k]]$period
			if (mdlType=="A"){
				SeasActual <<- rep(1,times=lenX)
				SeasFcast <- rep(1,times=lenXh)
			}else if (mdlType=="M"){
					SeasActual <<- rep(1,times=lenX)
					SeasFcast <- rep(1,times=lenXh)
			}else{
			}		
		}
		optParP <- as.data.frame(optParP)
		optParQ <- as.data.frame(optParQ)
		if (optParP[1,] == "opt"){
			parPstart <- 1
			parPend <- lenX
		}else {
			parPstart <- 1
			parPend <- nrow(optParP)
		}
		if (optParQ[1,] == "opt"){
			parQstart <- 1
			parQend <- lenX
		}else {
			parQstart <- 1
			parQend <- nrow(optParQ)
		}
		varATA <- paste(VarSeriesName,"ATAComponents",sep="")
		assign(varATA,matrix(NA,0,12), envir = .GlobalEnv)
		eval(parse(text= paste(varATA," <<- as.data.frame(",varATA,")",sep="")))
		eval(parse(text= paste("colnames(",varATA,") <<- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
		for (jj in parQstart:parQend){
			if (optParQ[1,]=="opt"){
				qh <- as.numeric(jj)-1
			}else{
				qh <- as.numeric(optParQ[jj,])
			}
			for (kk in parPstart:parPend){		
				if (optParP[1,]=="opt"){
					ph <- as.numeric(kk)
				}else{
					ph <- as.numeric(optParP[kk,])
				}
				if (ph==0){
					tempCalc <- as.data.frame(matrix(0,1,12))
					eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
					eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
				}else if (ph<qh){
					tempCalc <- as.data.frame(matrix(0,1,12))
					eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
					eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))	
				}else {
					if (mdlType=="A"){	
						coefph <- abs(ph/lenX)
						coefqh <- abs(qh/lenX)
						ATA <- ata.sub.func(X, ph, qh, lenX-1, VarSeriesName, mdlType, opterrType, SeasActual)
						T_1 <- ATA$trend
						S_1 <- ATA$level
						S <- coefph * X[lenX] + (1-coefph)*(S_1 + T_1)
						T <- coefqh * (S-S_1) + (1-coefqh) * (T_1)
						if (is.season ==1){
							SIValue <- SeasIndex[as.numeric(freqXh[1])]
						}else{
							SIValue <- 1
						}
						if (opterrType=="MAE"){	
							E_ErrType <- abs(((S + T) * SIValue) - Xh[1])
						}else if (opterrType=="MAPE"){
							E_ErrType <- abs((((S + T) * SIValue) - Xh[1])/Xh[1])*100
						}else if (opterrType=="sMAPE"){	
							E_ErrType <- (abs(((S + T) * SIValue) - Xh[1])/(abs((S + T) * SIValue) + abs(Xh[1])))*200
						}else if (opterrType=="MSE"){	
							E_ErrType <- (((S + T) * SIValue) - X[lenX])^2
						}else {
						}
						tempCalc <- matrix(c(1, qh, ph, lenX, coefqh, coefph, X[lenX], ((S + T) * SIValue), SIValue, S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
						eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
						eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
						for (h in 1:(lenXh-1)){
							if (is.season ==1){
								SIValue <- SeasIndex[as.numeric(freqXh[h+1])]
							}else{
								SIValue <- 1
							}
							if (opterrType=="MAE"){	
								E_ErrType <- abs(((S + ((h+1)*T)) * SIValue) - Xh[h+1])
							}else if (opterrType=="MAPE"){
								E_ErrType <- abs((((S + ((h+1)*T)) * SIValue) - Xh[h+1])/Xh[h+1])*100
							}else if (opterrType=="sMAPE"){	
								E_ErrType <- (abs(((S + ((h+1)*T)) * SIValue) - Xh[h+1])/(abs((S + ((h+1)*T)) * SIValue) + abs(Xh[h+1])))*200
							}else if (opterrType=="MSE"){	
								E_ErrType <- (((S + ((h+1)*T)) * SIValue) - Xh[h+1])^2
							}else {
							}
						tempCalc <- matrix(c(1, qh, ph, h, NA, NA, Xh[h], ((S + ((h+1)*T)) * SIValue), SIValue, NA, NA, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
						eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
						eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
						}
						hplus <- (freqXh[lenXh]%%max(freqXh))+1
						if (is.season ==1){
							SIValue <- SeasIndex[as.numeric(hplus)]
						}else{
							SIValue <- 1
						}
						tempCalc <- matrix(c(1, qh, ph, lenXh, NA, NA, Xh[lenXh], ((S + ((lenXh+1)*T)) * SIValue), SIValue, NA, NA, NA), nrow=1, ncol=12, byrow = TRUE)
						eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
						eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
					}else {
						if (mdlType=="M"){	
							coefph <- abs(ph/lenX)
							coefqh <- abs(qh/lenX)
							ATA <- ata.sub.func(X, ph, qh, lenX-1, VarSeriesName, mdlType, opterrType, SeasActual)
							T_1 <- ATA$trend
							S_1 <- ATA$level
							S <- coefph * X[lenX] + (1-coefph)*(S_1 * T_1)
							T <- coefqh * (S/S_1) + (1-coefqh) * (T_1)
							if (is.season ==1){
								SIValue <- SeasIndex[as.numeric(freqXh[1])]
							}else{
								SIValue <- 1
							}
							if (opterrType=="MAE"){	
								E_ErrType <- abs((S * T * SIValue) - Xh[1])
							}else if (opterrType=="MAPE"){
								E_ErrType <- abs(((S * T * SIValue) - Xh[1])/Xh[1])*100
							}else if (opterrType=="sMAPE"){	
								E_ErrType <- (abs((S * T * SIValue) - Xh[1])/(abs(S * T * SIValue) + abs(Xh[1])))*200
							}else if (opterrType=="MSE"){	
								E_ErrType <- ((S * T * SIValue) - Xh[1])^2
							}else {
							}
							tempCalc <- matrix(c(1, qh, ph, lenX, coefqh, coefph, X[lenX], (S * T * SIValue), SIValue, S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
							eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
							eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
							for (h in 1:(lenXh-1)){
								if (is.season ==1){
									SIValue <- SeasIndex[as.numeric(freqXh[h+1])]
								}else{
									SIValue <- 1
								}	
								if (opterrType=="MAE"){	
									E_ErrType <- abs((S * (T^(h+1)) * SIValue) - Xh[h+1])
								}else if (opterrType=="MAPE"){
									E_ErrType <- abs(((S * (T^(h+1)) * SIValue) - Xh[h+1])/Xh[h+1])*100
								}else if (opterrType=="sMAPE"){	
									E_ErrType <- (abs((S * (T^(h+1)) * SIValue) - Xh[h+1])/(abs(S * (T^(h+1)) * SIValue) + abs(Xh[h+1])))*200
								}else if (opterrType=="MSE"){	
									E_ErrType <- ((S * (T^(h+1)) * SIValue) - Xh[h+1])^2
								}else {
								}
								tempCalc <- matrix(c(1, qh, ph, h, NA, NA, Xh[h], (S * (T^(h+1)) * SIValue), SIValue, NA, NA, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
								eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
								eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
							}
							hplus <- (freqXh[lenXh]%%max(freqXh))+1
							if (is.season ==1){
								SIValue <- SeasIndex[as.numeric(hplus)]
							}else{
								SIValue <- 1
							}
							tempCalc <- matrix(c(1, qh, ph, lenXh, NA, NA, Xh[lenXh], (S * (T^(lenXh+1)) * SIValue), SIValue, NA, NA, NA), nrow=1, ncol=12, byrow = TRUE)
							eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
							eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
						}
					}
				}	
			}	
		}
		if (!is.null(nrow(optParQ))){
			lenQ <- nrow(optParQ)
		}else {
			lenQ <- length(optParQ)
		}
		if (!is.null(nrow(optParP))){
			lenP <- nrow(optParP)
		}else {
			lenP <- length(optParP)
		}
		if (optoptimPar=="optPQ" & lenQ == 1 & lenP == 1 & optParP[1,] != "opt" & optParQ[1,] != "opt"){
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"q\"]==optParQ,]", sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"p\"]==optParP,]", sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<-as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"Results) <<- c(\"Sample\",\"q\", \"p\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\", as.character(opterrType))",sep=""))) 
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}else if (optoptimPar == "optPQ" & ((lenQ > 1 & lenP == 1) | (lenQ == 1 & lenP > 1))){
			eval(parse(text= paste("PreOptError","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]",sep="")))
			eval(parse(text= paste("PreOptError","<- PreOptError[,c(\"q\",\"p\",\"",opterrType,"\")]",sep="")))
			PreOptError <- as.data.frame(cbind(PreOptError[,c(1,2)],as.data.frame(as.numeric(as.character(PreOptError[,3])))))
			colnames(PreOptError) <- c("q","p", "ErrType") 
			OptError <- melt(PreOptError, varnames=c("q","p"),measure.vars=c("ErrType"),na.rm = FALSE)
			OptMeanError <- dcast(OptError, q + p ~ variable,value.var="value",fun.aggregate=mean,na.rm=TRUE)
			OptMeanError <- OptMeanError[OptMeanError[,"p"]!=0,]
			ErrTypepq <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptMeanError[,"ErrType"],na.rm = TRUE),])
			ErrTypepq <- cbind(ErrTypepq,ErrTypepq[,1]+ErrTypepq[,2])
			colnames(ErrTypepq)<-c("q","p","ErrType","TotalPQ")
			eval(parse(text= paste(VarSeriesName,"BestOptimErrType <<- na.exclude(ErrTypepq[ErrTypepq[,\"TotalPQ\"]==max(ErrTypepq[,\"TotalPQ\"],na.rm = TRUE),])",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]", sep="")))
			PreResults <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreResults) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ResultsPre","<<- as.data.frame(PreResults)",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==1,]", sep="")))
			PreForecast <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreForecast) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ForecastsPre","<<- as.data.frame(PreForecast)",sep="")))
			eval(parse(text= paste(VarSeriesName,"ResultsPrep <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPrep,", VarSeriesName, "Forecasts)",sep="")))
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}else if (optoptimPar == "both" & (optfcast=="select" | optfcast=="comb")){
			eval(parse(text= paste("PreOptError","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]",sep="")))
			eval(parse(text= paste("PreOptError","<- PreOptError[,c(\"q\",\"p\",\"",opterrType,"\")]",sep="")))
			PreOptError <- as.data.frame(cbind(PreOptError[,c(1,2)],as.data.frame(as.numeric(as.character(PreOptError[,3])))))
			colnames(PreOptError) <- c("q","p", "ErrType") 
			OptError <- melt(PreOptError, varnames=c("q","p"),measure.vars=c("ErrType"),na.rm = FALSE)
			OptMeanError <- dcast(OptError, q + p ~ variable,value.var="value",fun.aggregate=mean,na.rm=TRUE)
			OptMeanError <- OptMeanError[OptMeanError[,"p"]!=0,]
			OptPstarError <- OptMeanError[OptMeanError[,"q"]==0,]
			OptPoneError <- OptMeanError[OptMeanError[,"q"]==1,]
			ErrTypePstar <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptPstarError[,"ErrType"],na.rm = TRUE),])
			ErrTypePstar <- cbind(ErrTypePstar,ErrTypePstar[,1]+ErrTypePstar[,2])
			colnames(ErrTypePstar)<-c("q","p","ErrType","TotalPQ")
			ErrTypePstar <- na.exclude(ErrTypePstar[ErrTypePstar[,"TotalPQ"]==max(ErrTypePstar[,"TotalPQ"],na.rm = TRUE),])
			ErrTypePone <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptPoneError[,"ErrType"],na.rm = TRUE),])
			ErrTypePone <- cbind(ErrTypePone,ErrTypePone[,1]+ErrTypePone[,2])
			colnames(ErrTypePone)<-c("q","p","ErrType","TotalPQ")
			ErrTypePone <- na.exclude(ErrTypePone[ErrTypePone[,"TotalPQ"]==max(ErrTypePone[,"TotalPQ"],na.rm = TRUE),])
			eval(parse(text= paste(VarSeriesName,"BestOptimErrType <<- rbind(ErrTypePstar, ErrTypePone)",sep="")))
			eval(parse(text= paste("rownames(",VarSeriesName,"BestOptimErrType) <<- c(\"pStar\", \"pOne\")",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]", sep="")))
			PreResults <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreResults) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ResultsPre","<<- as.data.frame(PreResults)",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==1,]", sep="")))
			PreForecast <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreForecast) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ForecastsPre","<<- as.data.frame(PreForecast)",sep="")))
			if (optfcast=="select"){
				eval(parse(text= paste("selectFcast <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pStar\",\"pOne\"),])",sep="")))
				Pselect <- selectFcast[selectFcast[,"ErrType"]==min(selectFcast[,"ErrType"]),]
				eval(parse(text= paste(VarSeriesName,"ResultsPre <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(Pselect[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(Pselect[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(Pselect[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(Pselect[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPre,", VarSeriesName, "Forecasts)",sep="")))
			}else if (optfcast=="comb"){
				eval(parse(text= paste("combResultStar <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pStar\"),])",sep="")))
				eval(parse(text= paste("combResultOne <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pOne\"),])",sep="")))
				eval(parse(text= paste(VarSeriesName,"ResultsPreStar <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(combResultStar[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(combResultStar[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"ResultsPreOne <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(combResultOne[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(combResultOne[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"ResultsComb","<<- rbind(", VarSeriesName, "ResultsPreStar,", VarSeriesName, "ResultsPreOne)",sep="")))
				eval(parse(text= paste("combFcast <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pStar\",\"pOne\"),])",sep="")))
				eval(parse(text= paste("CombFcastStar <- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(combFcast[1,2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(combFcast[1,1]),]",sep="")))
				eval(parse(text= paste("CombFcastOne <- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(combFcast[2,2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(combFcast[2,1]),]",sep="")))
				CombFcastPre <- as.data.frame(rowMeans(cbind(as.data.frame(as.numeric(as.character(CombFcastStar[,7]))),as.data.frame(as.numeric(as.character(CombFcastOne[,7]))))))
				PreCombFcast <- cbind(matrix("0&1",nrow(CombFcastStar),1),matrix(paste(combFcast[1,2],"&",combFcast[2,2],sep=""),nrow(CombFcastStar),1),as.data.frame(as.numeric(as.character(CombFcastStar[,6]))),CombFcastPre)
				Qcomb <- "0&1"
				Pcomb <- paste(combFcast[1,2],"&",combFcast[2,2],sep="")
				colnames(PreCombFcast, do.NULL = FALSE) 
				colnames(PreCombFcast) <- c("qComb", "pComb", "Actual", "ForecastComb") 
				ForecastsCombPre <<- as.data.frame(PreCombFcast)
				Err_Nrm <- (ForecastsCombPre[,4] - ForecastsCombPre[,3])
				if (opterrType=="MAE"){
					Err_ErrType <- abs((ForecastsCombPre[,4]) - ForecastsCombPre[,3])
				}else if (opterrType=="MAPE"){
					Err_ErrType <- abs((ForecastsCombPre[,4] - ForecastsCombPre[,3])/ForecastsCombPre[,3])*100
				}else if (opterrType=="sMAPE"){
					Err_ErrType <- (abs(ForecastsCombPre[,4] - ForecastsCombPre[,3])/(abs(ForecastsCombPre[,4]) + abs(ForecastsCombPre[,3])))*200
				}else if (opterrType=="MSE"){
					Err_ErrType <- (ForecastsCombPre[,4] - ForecastsCombPre[,3])^2
				}else {
				}
				ForecastsComb <- cbind(rep(1,times=(lenXh+1)),ForecastsCombPre[,1:2],rep(NA,times=(lenXh+1)),rep(NA,times=(lenXh+1)), ForecastsCombPre[,3:4],rep(NA,times=(lenXh+1)),rep(NA,times=(lenXh+1)),rep(NA,times=(lenXh+1)),Err_ErrType) 
				colnames(ForecastsComb, do.NULL = FALSE) 
				colnames(ForecastsComb) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType))
				eval(parse(text= paste(VarSeriesName,"ForecastsCP","<<- as.data.frame(ForecastsComb)",sep="")))
				eval(parse(text= paste("character_vars <- lapply(",VarSeriesName,"ForecastsCP, class) == \"numeric\"",sep="")))
				eval(parse(text= paste(VarSeriesName,"ForecastsCP[, character_vars] <<- lapply(",VarSeriesName,"ForecastsCP[, character_vars], as.factor)",sep="")))
				eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsComb,", VarSeriesName, "ForecastsCP)",sep="")))
			}else{
			}
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}else if (optoptimPar == "optPQ" & optfcast=="single"){
			eval(parse(text= paste("PreOptError","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]",sep="")))
			eval(parse(text= paste("PreOptError","<- PreOptError[,c(\"q\",\"p\",\"",opterrType,"\")]",sep="")))
			PreOptError <- as.data.frame(cbind(PreOptError[,c(1,2)],as.data.frame(as.numeric(as.character(PreOptError[,3])))))
			colnames(PreOptError) <- c("q","p", "ErrType") 
			OptError <- melt(PreOptError, varnames=c("q","p"),measure.vars=c("ErrType"),na.rm = FALSE)
			OptMeanError <- dcast(OptError, q + p ~ variable,value.var="value",fun.aggregate=mean,na.rm=FALSE)
			OptMeanError <- OptMeanError[OptMeanError[,"p"]!=0,]
			ErrTypepq <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptMeanError[,"ErrType"],na.rm = FALSE),])
			ErrTypepq <- cbind(ErrTypepq,ErrTypepq[,1]+ErrTypepq[,2])
			colnames(ErrTypepq)<-c("q","p","ErrType","TotalPQ")
			eval(parse(text= paste(VarSeriesName,"BestOptimErrType <<- na.exclude(ErrTypepq[ErrTypepq[,\"TotalPQ\"]==max(ErrTypepq[,\"TotalPQ\"],na.rm = FALSE),])",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]", sep="")))
			PreResults <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreResults) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ResultsPre","<<- as.data.frame(PreResults)",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==1,]", sep="")))
			PreForecast <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreForecast) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ForecastsPre","<<- as.data.frame(PreForecast)",sep="")))
			eval(parse(text= paste(VarSeriesName,"ResultsPrep <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPrep,", VarSeriesName, "Forecasts)",sep="")))
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}else if (optoptimPar == "pStar" & optfcast=="single"){
			eval(parse(text= paste("PreOptError","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]",sep="")))
			eval(parse(text= paste("PreOptError","<- PreOptError[,c(\"q\",\"p\",\"",opterrType,"\")]",sep="")))
			PreOptError <- as.data.frame(cbind(PreOptError[,c(1,2)],as.data.frame(as.numeric(as.character(PreOptError[,3])))))
			colnames(PreOptError) <- c("q","p", "ErrType") 
			OptError <- melt(PreOptError, varnames=c("q","p"),measure.vars=c("ErrType"),na.rm = FALSE)
			OptMeanError <- dcast(OptError, q + p ~ variable,value.var="value",fun.aggregate=mean,na.rm=TRUE)
			OptMeanError <- OptMeanError[OptMeanError[,"p"]!=0,]
			OptPstarError <- OptMeanError[OptMeanError[,"q"]==0,]
			ErrTypePstar <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptPstarError[,"ErrType"],na.rm = TRUE),])
			ErrTypePstar <- cbind(ErrTypePstar,ErrTypePstar[,1]+ErrTypePstar[,2])
			colnames(ErrTypePstar)<-c("q","p","ErrType","TotalPQ")
			eval(parse(text= paste(VarSeriesName,"BestOptimErrType <<- na.exclude(ErrTypePstar[ErrTypePstar[,\"TotalPQ\"]==max(ErrTypePstar[,\"TotalPQ\"],na.rm = TRUE),])", sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]", sep="")))
			PreResults <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreResults) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ResultsPre","<<- as.data.frame(PreResults)",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==1,]", sep="")))
			PreForecast <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreForecast) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ForecastsPre","<<- as.data.frame(PreForecast)",sep="")))
			eval(parse(text= paste(VarSeriesName,"ResultsPrep <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPrep,", VarSeriesName, "Forecasts)",sep="")))
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}else if (optoptimPar == "pOne" & optfcast=="single"){
			eval(parse(text= paste("PreOptError","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]",sep="")))
			eval(parse(text= paste("PreOptError","<- PreOptError[,c(\"q\",\"p\",\"",opterrType,"\")]",sep="")))
			PreOptError <- as.data.frame(cbind(PreOptError[,c(1,2)],as.data.frame(as.numeric(as.character(PreOptError[,3])))))
			colnames(PreOptError) <- c("q","p", "ErrType") 
			OptError <- melt(PreOptError, varnames=c("q","p"),measure.vars=c("ErrType"),na.rm = FALSE)
			OptMeanError <- dcast(OptError, q + p ~ variable,value.var="value",fun.aggregate=mean,na.rm=TRUE)
			OptMeanError <- OptMeanError[OptMeanError[,"p"]!=0,]
			OptPoneError <- OptMeanError[OptMeanError[,"q"]==1,]
			ErrTypePone <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptPoneError[,"ErrType"],na.rm = TRUE),])
			ErrTypePone <- cbind(ErrTypePone,ErrTypePone[,1]+ErrTypePone[,2])
			colnames(ErrTypePone)<-c("q","p","ErrType","TotalPQ")
			eval(parse(text= paste(VarSeriesName,"BestOptimErrType <<- na.exclude(ErrTypePone[ErrTypePone[,\"TotalPQ\"]==max(ErrTypePone[,\"TotalPQ\"],na.rm = TRUE),])", sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]", sep="")))
			PreResults <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreResults) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ResultsPre","<<- as.data.frame(PreResults)",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==1,]", sep="")))
			PreForecast <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreForecast) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ForecastsPre","<<- as.data.frame(PreForecast)",sep="")))
			eval(parse(text= paste(VarSeriesName,"ResultsPrep <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPrep,", VarSeriesName, "Forecasts)",sep="")))
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}else if (optoptimPar == "pStarQ" & optfcast=="single"){
			eval(parse(text= paste("PreOptError","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]",sep="")))
			eval(parse(text= paste("PreOptError","<- PreOptError[,c(\"q\",\"p\",\"",opterrType,"\")]",sep="")))
			PreOptError <- as.data.frame(cbind(PreOptError[,c(1,2)],as.data.frame(as.numeric(as.character(PreOptError[,3])))))
			colnames(PreOptError) <- c("q","p", "ErrType") 
			OptError <- melt(PreOptError, varnames=c("q","p"),measure.vars=c("ErrType"),na.rm = FALSE)
			OptMeanError <- dcast(OptError, q + p ~ variable,value.var="value",fun.aggregate=mean,na.rm=TRUE)
			OptMeanError <- OptMeanError[OptMeanError[,"p"]!=0,]
			OptPstarError <- OptMeanError[OptMeanError[,"q"]==0,]
			OptPoneError <- OptMeanError[OptMeanError[,"q"]==1,]
			ErrTypePstar <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptPstarError[,"ErrType"],na.rm = TRUE),])
			ErrTypePstar <- cbind(ErrTypePstar,ErrTypePstar[,1]+ErrTypePstar[,2])
			colnames(ErrTypePstar)<-c("q","p","ErrType","TotalPQ")
			ErrTypePstar <- na.exclude(ErrTypePstar[ErrTypePstar[,"TotalPQ"]==max(ErrTypePstar[,"TotalPQ"],na.rm = TRUE),])
			OptQPstarErrType <- OptMeanError[OptMeanError[,"p"]==ErrTypePstar[1,2],]
			ErrTypeQPstar <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptQPstarErrType[,"ErrType"],na.rm = TRUE),])
			ErrTypeQPstar <- cbind(ErrTypeQPstar,ErrTypeQPstar[,1]+ErrTypeQPstar[,2])
			colnames(ErrTypeQPstar)<-c("q","p","ErrType","TotalPQ")
			eval(parse(text= paste(VarSeriesName,"BestOptimErrType <<- na.exclude(ErrTypeQPstar[ErrTypeQPstar[,\"TotalPQ\"]==max(ErrTypeQPstar[,\"TotalPQ\"],na.rm = TRUE),])",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]", sep="")))
			PreResults <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreResults) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ResultsPre","<<- as.data.frame(PreResults)",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==1,]", sep="")))
			PreForecast <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreForecast) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ForecastsPre","<<- as.data.frame(PreForecast)",sep="")))
			eval(parse(text= paste(VarSeriesName,"ResultsPrep <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(",VarSeriesName,"BestOptimErrType[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(",VarSeriesName,"BestOptimErrType[1]),]",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPrep,", VarSeriesName, "Forecasts)",sep="")))
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}else{
			eval(parse(text= paste("PreOptError","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]",sep="")))
			eval(parse(text= paste("PreOptError","<- PreOptError[,c(\"q\",\"p\",\"",opterrType,"\")]",sep="")))
			PreOptError <- as.data.frame(cbind(PreOptError[,c(1,2)],as.data.frame(as.numeric(as.character(PreOptError[,3])))))
			colnames(PreOptError) <- c("q","p", "ErrType") 
			OptError <- melt(PreOptError, varnames=c("q","p"),measure.vars=c("ErrType"),na.rm = FALSE)
			OptMeanError <- dcast(OptError, q + p ~ variable,value.var="value",fun.aggregate=mean,na.rm=TRUE)
			OptMeanError <- OptMeanError[OptMeanError[,"p"]!=0,]
			OptPstarError <- OptMeanError[OptMeanError[,"q"]==0,]
			OptPoneError <- OptMeanError[OptMeanError[,"q"]==1,]
			
			ErrTypePstar <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptPstarError[,"ErrType"],na.rm = TRUE),])
			ErrTypePstar <- cbind(ErrTypePstar,ErrTypePstar[,1]+ErrTypePstar[,2])
			colnames(ErrTypePstar)<-c("q","p","ErrType","TotalPQ")
			ErrTypePstar <- na.exclude(ErrTypePstar[ErrTypePstar[,"TotalPQ"]==max(ErrTypePstar[,"TotalPQ"],na.rm = TRUE),])
			
			OptQPstarErrType <- OptMeanError[OptMeanError[,"p"]==ErrTypePstar[1,2],]
			ErrTypeQPstar <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptQPstarErrType[,"ErrType"],na.rm = TRUE),])
			ErrTypeQPstar <- cbind(ErrTypeQPstar,ErrTypeQPstar[,1]+ErrTypeQPstar[,2])
			colnames(ErrTypeQPstar)<-c("q","p","ErrType","TotalPQ")
			ErrTypeQPstar <- na.exclude(ErrTypeQPstar[ErrTypeQPstar[,"TotalPQ"]==max(ErrTypeQPstar[,"TotalPQ"],na.rm = TRUE),])
			
			ErrTypepq <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptMeanError[,"ErrType"],na.rm = TRUE),])
			ErrTypepq <- cbind(ErrTypepq,ErrTypepq[,1]+ErrTypepq[,2])
			colnames(ErrTypepq)<-c("q","p","ErrType","TotalPQ")
			ErrTypepq <- na.exclude(ErrTypepq[ErrTypepq[,"TotalPQ"]==max(ErrTypepq[,"TotalPQ"],na.rm = TRUE),])
			
			ErrTypePone <- na.exclude(OptMeanError[OptMeanError[,"ErrType"]==min(OptPoneError[,"ErrType"],na.rm = TRUE),])
			ErrTypePone <- cbind(ErrTypePone,ErrTypePone[,1]+ErrTypePone[,2])
			colnames(ErrTypePone)<-c("q","p","ErrType","TotalPQ")
			ErrTypePone <- na.exclude(ErrTypePone[ErrTypePone[,"TotalPQ"]==max(ErrTypePone[,"TotalPQ"],na.rm = TRUE),])
			
			BestErrTypeOpt <- rbind(ErrTypepq, ErrTypePstar, ErrTypePone, ErrTypeQPstar)
			rownames(BestErrTypeOpt) <- c("optPQ","pStar", "pOne", "pStarQ")
			eval(parse(text= paste(VarSeriesName,"BestOptimErrType <<- BestErrTypeOpt[,c(1,2,3,4)]",sep="")))
			if (optoptimPar=="optPQ" & optoptimPar=="both"){
				eval(parse(text= paste("minoptBest <- min(",VarSeriesName,"BestOptimErrType[,\"ErrType\"],na.rm = TRUE)",sep=""))) 
				eval(parse(text= paste("optBestPar <- na.exclude(",VarSeriesName,"BestOptimErrType[",VarSeriesName,"BestOptimErrType[,\"ErrType\"]==minoptBest,])",sep="")))
				optBestPar <- na.exclude(optBestPar[optBestPar[,"TotalPQ"]==max(optBestPar[,"TotalPQ"],na.rm = TRUE),])
				optBestPar <- optBestPar[rownames(unique(optBestPar[,c(3,4)],fromLast = TRUE)),]
			}else {
				eval(parse(text= paste("optBestPar <- na.exclude(",VarSeriesName,"BestOptimErrType[row.names(",VarSeriesName,"BestOptimErrType)==\"",optoptimPar,"\",])",sep="")))
				optBestPar <- na.exclude(optBestPar[optBestPar[,"TotalPQ"]==max(optBestPar[,"TotalPQ"],na.rm = TRUE),])
				optBestPar <- optBestPar[rownames(unique(optBestPar[,c(3,4)],fromLast = TRUE)),]
			}
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==0,]", sep="")))
			PreResults <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreResults) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ResultsPre","<<- as.data.frame(PreResults)",sep="")))
			eval(parse(text= paste("PreATAComponents","<-",varATA,"[",varATA,"[,\"Sample\"]==1,]", sep="")))
			PreForecast <- as.data.frame(cbind(PreATAComponents[,1:3],PreATAComponents[,5:12]))
			colnames(PreForecast) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType)) 
			eval(parse(text= paste(VarSeriesName,"ForecastsPre","<<- as.data.frame(PreForecast)",sep="")))
			if (optfcast=="select"){
				eval(parse(text= paste("selectFcast <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pStar\",\"pOne\"),])",sep="")))
				Pselect <- selectFcast[selectFcast[,"ErrType"]==min(selectFcast[,"ErrType"]),]
				eval(parse(text= paste(VarSeriesName,"ResultsPre <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(Pselect[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(Pselect[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(Pselect[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(Pselect[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPre,", VarSeriesName, "Forecasts)",sep="")))
			}else if (optfcast=="comb"){
				eval(parse(text= paste("combResultStar <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pStar\"),])",sep="")))
				eval(parse(text= paste("combResultOne <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pOne\"),])",sep="")))
				eval(parse(text= paste(VarSeriesName,"ResultsPreStar <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(combResultStar[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(combResultStar[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"ResultsPreOne <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(combResultOne[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(combResultOne[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"ResultsComb","<<- rbind(", VarSeriesName, "ResultsPreStar,", VarSeriesName, "ResultsPreOne)",sep="")))
				eval(parse(text= paste("combFcast <- as.data.frame(", VarSeriesName,"BestOptimErrType[c(\"pStar\",\"pOne\"),])",sep="")))
				eval(parse(text= paste("CombFcastStar <- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(combFcast[1,2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(combFcast[1,1]),]",sep="")))
				eval(parse(text= paste("CombFcastOne <- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(combFcast[2,2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(combFcast[2,1]),]",sep="")))
				CombFcastPre <- as.data.frame(rowMeans(cbind(as.data.frame(as.numeric(as.character(CombFcastStar[,7]))),as.data.frame(as.numeric(as.character(CombFcastOne[,7]))))))
				PreCombFcast <- cbind(matrix("0&1",nrow(CombFcastStar),1),matrix(paste(combFcast[1,2],"&",combFcast[2,2],sep=""),nrow(CombFcastStar),1),as.data.frame(as.numeric(as.character(CombFcastStar[,6]))),CombFcastPre)
				Qcomb <- "0&1"
				Pcomb <- paste(combFcast[1,2],"&",combFcast[2,2],sep="")
				colnames(PreCombFcast, do.NULL = FALSE) 
				colnames(PreCombFcast) <- c("qComb", "pComb", "Actual", "ForecastComb") 
				ForecastsCombPre <<- as.data.frame(PreCombFcast)
				Err_Nrm <- (ForecastsCombPre[,4] - ForecastsCombPre[,3])
				if (opterrType=="MAE"){
					Err_ErrType <- abs((ForecastsCombPre[,4]) - ForecastsCombPre[,3])
				}else if (opterrType=="MAPE"){
					Err_ErrType <- abs((ForecastsCombPre[,4] - ForecastsCombPre[,3])/ForecastsCombPre[,3])*100
				}else if (opterrType=="sMAPE"){
					Err_ErrType <- (abs(ForecastsCombPre[,4] - ForecastsCombPre[,3])/(abs(ForecastsCombPre[,4]) + abs(ForecastsCombPre[,3])))*200
				}else if (opterrType=="MSE"){
					Err_ErrType <- (ForecastsCombPre[,4] - ForecastsCombPre[,3])^2
				}else {
				}
				ForecastsComb <- cbind(rep(1,times=(lenXh+1)),ForecastsCombPre[,1:2],rep(NA,times=(lenXh+1)),rep(NA,times=(lenXh+1)), ForecastsCombPre[,3:4],rep(NA,times=(lenXh+1)),rep(NA,times=(lenXh+1)),rep(NA,times=(lenXh+1)),Err_ErrType) 
				colnames(ForecastsComb, do.NULL = FALSE) 
				colnames(ForecastsComb) <- c("Sample","q", "p", "coefq", "coefp", "Actual", "Fit", "SeasComp", "Level", "Trend", as.character(opterrType))
				eval(parse(text= paste(VarSeriesName,"ForecastsCP","<<- as.data.frame(ForecastsComb)",sep="")))
				eval(parse(text= paste("character_vars <- lapply(",VarSeriesName,"ForecastsCP, class) == \"numeric\"",sep="")))
				eval(parse(text= paste(VarSeriesName,"ForecastsCP[, character_vars] <<- lapply(",VarSeriesName,"ForecastsCP[, character_vars], as.factor)",sep="")))
				eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsComb,", VarSeriesName, "ForecastsCP)",sep="")))
				remove(BestErrTypeOpt)
			}else {
				eval(parse(text= paste(VarSeriesName,"ResultsPrep <<- ",VarSeriesName,"ResultsPre[",VarSeriesName,"ResultsPre[,\"p\"]==as.numeric(optBestPar[2])&",VarSeriesName,"ResultsPre[,\"q\"]==as.numeric(optBestPar[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"Forecasts <<- ",VarSeriesName,"ForecastsPre[",VarSeriesName,"ForecastsPre[,\"p\"]==as.numeric(optBestPar[2])&",VarSeriesName,"ForecastsPre[,\"q\"]==as.numeric(optBestPar[1]),]",sep="")))
				eval(parse(text= paste(VarSeriesName,"Results","<<- rbind(", VarSeriesName, "ResultsPrep,", VarSeriesName, "Forecasts)",sep="")))
			}	
			addinColNames <- c("Series", "N", "Nh", "Category", "Frequency", "Seasonality")
			eval(parse(text= paste("rColNames <- colnames(",VarSeriesName,"Results)",sep="")))
			resultColNames <- as.vector(cbind(t(addinColNames),t(rColNames)))
			eval(parse(text= paste(VarSeriesName,"tmpResults <<- cbind(rep(\"",VarSeriesName,"\",times=nrow(", VarSeriesName,"Results)),rep(\"",obsvNum,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									fHorizon,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",tscategory,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",obsPeriod,"\",times=nrow(", VarSeriesName,"Results)), rep(\"",
									is.season,"\",times=nrow(", VarSeriesName,"Results)),", VarSeriesName, "Results)",sep="")))
			eval(parse(text= paste("colnames(",VarSeriesName,"tmpResults) <<- resultColNames",sep="")))
			eval(parse(text= paste(VarSeriesName,"Results <<- ",VarSeriesName,"tmpResults",sep="")))
		}
		WD <- getwd()
		subDir <- "ATA Method"
		ifelse(!dir.exists(file.path(WD, subDir)), dir.create(file.path(WD, subDir)), FALSE)	  
		txtWrite <- paste("write.csv(",VarSeriesName,"Results, file=\"",WD,"/",subDir,"/",VarSeriesName,"Results.csv\")",sep="")	
		eval(parse(text= txtWrite))
		eval(parse(text= paste("in_sample <- sum(na.exclude(as.data.frame(as.numeric(as.character(",VarSeriesName,"Results[",VarSeriesName,"Results[,\"Sample\"]==0,\"",opterrType,"\"])))))/(obsvNum-1)",sep="")))
		eval(parse(text= paste("out_sample <- sum(na.exclude(as.data.frame(as.numeric(as.character(",VarSeriesName,"Results[",VarSeriesName,"Results[,\"Sample\"]==1,\"",opterrType,"\"])))))/fHorizon",sep="")))
		if (optfcast=="comb"){
			eval(parse(text= paste("PreM3Results <- as.data.frame(cbind(VarSeriesName,obsvNum,fHorizon,is.season,tscategory,obsPeriod,in_sample,out_sample,Qcomb,Pcomb,as.data.frame(t(",VarSeriesName,"Results[",VarSeriesName,"Results[,\"Sample\"]==1,\"Fit\"])),as.data.frame(matrix(NA,1,maxfHorizon-fHorizon))))",sep="")))
		}else {
			eval(parse(text= paste("PreM3Results <- as.data.frame(cbind(VarSeriesName,obsvNum,fHorizon,is.season,tscategory,obsPeriod,in_sample,out_sample,",VarSeriesName,"Results[1,8],",VarSeriesName,"Results[1,9],as.data.frame(t(",VarSeriesName,"Results[",VarSeriesName,"Results[,\"Sample\"]==1,\"Fit\"])),as.data.frame(matrix(NA,1,maxfHorizon-fHorizon))))",sep="")))
		}
		colnames(PreM3Results) <- c("Series", "N", "Nh", "Seasonality", "Category", "Frequency","in_sample","out_sample","q", "p","F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14","F15","F16","F17","F18","F19")
		M3Results <<- as.data.frame(rbind(M3Results,PreM3Results))
		
		if (optoptimPar=="optPQ" & lenQ == 1 & lenP == 1 & optParP[1,] != "opt" & optParQ[1,] != "opt"){
			remove(PreATAComponents,PreM3Results,tempCalc)
			eval(parse(text= paste("remove(SeasActual,",VarSeriesName,"ATAComponents,", VarSeriesName,"Results,", VarSeriesName,"tmpResults, envir=.GlobalEnv)",sep="")))						
		}else if (optoptimPar=="optPQ"){
			remove(PreATAComponents,PreM3Results,tempCalc,PreResults,PreForecast,ErrTypepq,PreOptError,OptError,OptMeanError)
			eval(parse(text= paste("remove(SeasActual,",VarSeriesName,"ATAComponents,", VarSeriesName,"ForecastsPre,", VarSeriesName, "ResultsPrep,", VarSeriesName, "ResultsPre,", VarSeriesName, "Forecasts,", VarSeriesName,"Results,", VarSeriesName,"tmpResults,",VarSeriesName,"BestOptimErrType, envir=.GlobalEnv)",sep="")))						
		}else if (optoptimPar=="pStar"){
			remove(PreATAComponents,PreM3Results,tempCalc,PreResults,PreForecast,ErrTypePstar,OptPstarError,PreOptError,OptError,OptMeanError)
			eval(parse(text= paste("remove(SeasActual,",VarSeriesName,"ATAComponents,", VarSeriesName,"ForecastsPre,", VarSeriesName, "ResultsPrep,", VarSeriesName, "ResultsPre,", VarSeriesName, "Forecasts,", VarSeriesName,"Results,", VarSeriesName,"tmpResults,",VarSeriesName,"BestOptimErrType, envir=.GlobalEnv)",sep="")))						
		}else if (optoptimPar=="pOne"){
			remove(PreATAComponents,PreM3Results,tempCalc,PreResults,PreForecast,ErrTypePone,OptPoneError,PreOptError,OptError,OptMeanError)
			eval(parse(text= paste("remove(SeasActual,",VarSeriesName,"ATAComponents,", VarSeriesName,"ForecastsPre,", VarSeriesName, "ResultsPrep,", VarSeriesName, "ResultsPre,", VarSeriesName, "Forecasts,", VarSeriesName,"Results,", VarSeriesName,"tmpResults,",VarSeriesName,"BestOptimErrType, envir=.GlobalEnv)",sep="")))						
		}else if (optoptimPar=="pStarQ"){
			remove(PreATAComponents,PreM3Results,tempCalc,PreResults,PreForecast,ErrTypeQPstar,PreOptError,OptError,OptMeanError,OptQPstarErrType)
			eval(parse(text= paste("remove(SeasActual,",VarSeriesName,"ATAComponents,", VarSeriesName,"ForecastsPre,", VarSeriesName, "ResultsPrep,", VarSeriesName, "ResultsPre,", VarSeriesName, "Forecasts,", VarSeriesName,"Results,", VarSeriesName,"tmpResults,",VarSeriesName,"BestOptimErrType, envir=.GlobalEnv)",sep="")))						
		}else if (optoptimPar=="both"){
			remove(PreATAComponents,PreM3Results,tempCalc,PreResults,PreForecast,ErrTypePstar,OptPstarError,ErrTypePone,OptPoneError,PreOptError,OptError,OptMeanError)
			if (optfcast=="select"){
				eval(parse(text= paste("remove(SeasActual,",VarSeriesName,"ATAComponents,", VarSeriesName,"ForecastsPre,", VarSeriesName, "ResultsPre,", VarSeriesName, "Forecasts,", VarSeriesName,"Results,", VarSeriesName,"tmpResults,",VarSeriesName,"BestOptimErrType, envir=.GlobalEnv)",sep="")))
			}
			if (optfcast=="comb"){
				eval(parse(text= paste("remove(SeasActual,ForecastsCombPre,",VarSeriesName,"ATAComponents,", VarSeriesName,"ForecastsPre,", VarSeriesName, "ResultsPre,", VarSeriesName,"Results,", VarSeriesName,"tmpResults,",VarSeriesName,"BestOptimErrType,",VarSeriesName, "ResultsPreStar,", VarSeriesName, "ResultsComb,", VarSeriesName, "ResultsPreOne,", VarSeriesName, "ForecastsCP, envir=.GlobalEnv)",sep="")))
			}
		}else{
		
		}	
		if (is.season ==1){
			remove(SeasIndex,envir=.GlobalEnv)
		}
		gc()
	}
} 