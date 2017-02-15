ataResultOutputCSV <- function(sRange,WD,subDir)
{
	if (sRange[1] == "full"){
		colSeries <- length(M3)
	}else{
		colsRange <- as.data.frame(sRange)
		colSeries <- nrow(colsRange)
	}	
	for (ll in 1:colSeries){
		if (sRange[1]=="full"){
			k <- ll
		}else{
			k <- as.numeric(colsRange[ll,])
		}
		VarSeriesName <- M3[[k]]$sn 
		if (ll==1){
			txtRead <- paste("SN_1 <- read.csv(file=\"",WD,"/",subDir,"/",VarSeriesName,"Results.csv\")",sep="")	
			eval(parse(text= txtRead))
			varResult <- SN_1
			eval(parse(text= paste("remove(SN_",ll,")",sep="")))
		}else{
			txtRead <- paste("SN_",ll," <- read.csv(file=\"",WD,"/",subDir,"/",VarSeriesName,"Results.csv\")",sep="")	
			eval(parse(text= txtRead))	
			eval(parse(text= paste("varResult <- rbind(varResult,SN_", ll,")", sep="")))
			eval(parse(text= paste("remove(SN_",ll,")",sep="")))
		}	
		file.remove(paste(WD,"/",subDir,"/",VarSeriesName,"Results.csv",sep=""))
	}
	M3ResultsNew <- unique(varResult)
	M3ForecastsNew <- M3ResultsNew[M3ResultsNew[,"Sample"]==1,]
	ifelse(!dir.exists(file.path(WD, subDir)), dir.create(file.path(WD, subDir)), FALSE)	  
	txtWrite <- paste("write.csv(M3ResultsNew, file=\"",WD,"/",subDir,"/","M3ResultsNew.csv\")",sep="")	
	eval(parse(text= txtWrite))
	txtWrite <- paste("write.csv(M3ForecastsNew, file=\"",WD,"/",subDir,"/","M3ForecastsNew.csv\")",sep="")	
	eval(parse(text= txtWrite))
}