ata.sub.func <- function(X, pk, qk, i, SerieName, mdlType, opterrType, seasIndex)
{
# X  			: nx1 matrix and ts object
# pk 			: Level parameter
# qk 			: Trend parameter
# i  			: i th serie in the Mcomp data set
# SerieName		: The name of i th serie 
# mdlType		: Additive or Multiplicative
# opterrType	: Selected Error Model
# seasIndex		: Seasonal component of the decomposed serie

	varATA <- paste(SerieName,"ATAComponents",sep="")
	if (i == 1){
		if (mdlType=="A"){
			S <- X[i]
			T <- 0
			if (opterrType=="MAE"){
				E_ErrType <- abs((S + T) - X[i+1])
			}else if (opterrType=="MAPE"){
				E_ErrType <- abs(((S + T) - X[i+1])/X[i+1])*100
			}else if (opterrType=="sMAPE"){
				E_ErrType <- (abs((S + T) - X[i+1])/(abs(S + T) + abs(X[i+1])))*200
			}else if (opterrType=="MSE"){
				E_ErrType <- ((S + T) - X[i+1])^2
			}else {
			}
			tempCalc <- matrix(c(0, qk, pk, i, NA, NA, X[i], (S + T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
			eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
			eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
		}else {
			if (mdlType=="M"){
				S <- X[i]
				T <- 1
				if (opterrType=="MAE"){
					E_ErrType <- abs((S * T) - X[i+1])
				}else if (opterrType=="MAPE"){
					E_ErrType <- abs(((S * T) - X[i+1])/X[i+1])*100
				}else if (opterrType=="sMAPE"){
					E_ErrType <- (abs((S * T) - X[i+1])/(abs(S * T) + abs(X[i+1])))*200
				}else if (opterrType=="MSE"){
					E_ErrType <- ((S * T) - X[i+1])^2
				}else {
				}
				tempCalc <- matrix(c(0, qk, pk, i, NA, NA, X[i], (S * T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
				eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
				eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
			}
		}
		my_list <- list("level" = S, "trend" = T)
		return(my_list)
	}else if (i<=pk & i<=qk & pk>=qk){
		if (mdlType=="A"){
			ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
			T_1 <- ATA$trend
			S_1 <- ATA$level	
			S <- X[i]
			T <- X[i] - X[i-1]
			if (opterrType=="MAE"){
				E_ErrType <- abs((S + T) - X[i+1])
			}else if (opterrType=="MAPE"){
				E_ErrType <- abs(((S + T) - X[i+1])/X[i+1])*100
			}else if (opterrType=="sMAPE"){
				E_ErrType <- (abs((S + T) - X[i+1])/(abs(S + T) + abs(X[i+1])))*200
			}else if (opterrType=="MSE"){
				E_ErrType <- ((S + T) - X[i+1])^2
			}else {
			}
			tempCalc <- matrix(c(0, qk, pk, i, NA, NA, X[i], (S + T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
			eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
			eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
		}else{
			if (mdlType=="M"){
				ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
				T_1 <- ATA$trend
				S_1 <- ATA$level
				S <- X[i]
				T <- X[i] / X[i-1]
				if (opterrType=="MAE"){
					E_ErrType <- abs((S * T) - X[i+1])
				}else if (opterrType=="MAPE"){
					E_ErrType <- abs(((S * T) - X[i+1])/X[i+1])*100
				}else if (opterrType=="sMAPE"){
					E_ErrType <- (abs((S * T) - X[i+1])/(abs(S * T) + abs(X[i+1])))*200
				}else if (opterrType=="MSE"){
					E_ErrType <- ((S * T) - X[i+1])^2
				}else {
				}
				tempCalc <- matrix(c(0, qk, pk, i, NA, NA, X[i], (S * T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
				eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
				eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
			}
		}
		my_list <- list("level" = S, "trend" = T)
		return(my_list)
	}else if (i<=pk & i>qk & pk>=qk){
		if (mdlType=="A"){
			coefqk <- abs(qk/i)
			ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
			T_1 <- ATA$trend
			S_1 <- ATA$level	
			S <- X[i]
			T <- coefqk * (S-S_1) + (1-coefqk) * (T_1)
			if (opterrType=="MAE"){
				E_ErrType <- abs((S + T) - X[i+1])
			}else if (opterrType=="MAPE"){
				E_ErrType <- abs(((S + T) - X[i+1])/X[i+1])*100
			}else if (opterrType=="sMAPE"){
				E_ErrType <- (abs((S + T) - X[i+1])/(abs(S + T) + abs(X[i+1])))*200
			}else if (opterrType=="MSE"){
				E_ErrType <- ((S + T) - X[i+1])^2
			}else {
			}
			tempCalc <- matrix(c(0, qk, pk, i, coefqk, NA, X[i], (S + T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
			eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
			eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
		}else {
			if (mdlType=="M"){
				coefqk <- abs(qk/i)
				ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
				T_1 <- ATA$trend
				S_1 <- ATA$level	
				S <- X[i]
				T <- coefqk * (S/S_1) + (1-coefqk) * (T_1)
				if (opterrType=="MAE"){
					E_ErrType <- abs((S * T) - X[i+1])
				}else if (opterrType=="MAPE"){
					E_ErrType <- abs(((S * T) - X[i+1])/X[i+1])*100
				}else if (opterrType=="sMAPE"){
					E_ErrType <- (abs((S * T) - X[i+1])/(abs(S * T) + abs(X[i+1])))*200
				}else if (opterrType=="MSE"){
					E_ErrType <- ((S * T) - X[i+1])^2
				}else {
				}
				tempCalc <- matrix(c(0, qk, pk, i, coefqk, NA, X[i], (S * T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
				eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
				eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))						
			}
		}
		my_list <- list("level" = S, "trend" = T)
		return(my_list)
	}else if (i>pk & i<=qk & pk>=qk){
		if (mdlType=="A"){
			coefpk <- abs(pk/i)
			ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
			T_1 <- ATA$trend
			S_1 <- ATA$level
			S <- coefpk * X[i] + (1-coefpk)*(S_1 + T_1)
			T <- X[i] - X[i-1]
			if (opterrType=="MAE"){
				E_ErrType <- abs((S + T) - X[i+1])
			}else if (opterrType=="MAPE"){
				E_ErrType <- abs(((S + T) - X[i+1])/X[i+1])*100
			}else if (opterrType=="sMAPE"){
				E_ErrType <- (abs((S + T) - X[i+1])/(abs(S + T) + abs(X[i+1])))*200
			}else if (opterrType=="MSE"){
				E_ErrType <- ((S + T) - X[i+1])^2
			}else {
			}
			tempCalc <- matrix(c(0, qk, pk, i, NA, coefpk, X[i], (S + T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
			eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
			eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))			
		}else {
			if (mdlType=="M"){
				coefpk <- abs(pk/i)
				ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
				T_1 <- ATA$trend
				S_1 <- ATA$level
				S <- coefpk * X[i] + (1-coefpk)*(S_1 * T_1)
				T <- X[i] / X[i-1]
				if (opterrType=="MAE"){
					E_ErrType <- abs((S * T) - X[i+1])
				}else if (opterrType=="MAPE"){
					E_ErrType <- abs(((S * T) - X[i+1])/X[i+1])*100
				}else if (opterrType=="sMAPE"){
					E_ErrType <- (abs((S * T) - X[i+1])/(abs(S * T) + abs(X[i+1])))*200
				}else if (opterrType=="MSE"){
					E_ErrType <- ((S * T) - X[i+1])^2
				}else {
				}
				tempCalc <- matrix(c(0, qk, pk, i, NA, coefpk, X[i], (S * T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE) 
				eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
				eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))	
			}
		}
		my_list <- list("level" = S, "trend" = T)
		return(my_list)
	}else if (i>pk & i>qk & pk>=qk){
		if (mdlType=="A"){
			coefpk <- abs(pk/i)
			coefqk <- abs(qk/i)
			ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
			T_1 <- ATA$trend
			S_1 <- ATA$level
			S <- coefpk * X[i] + (1-coefpk)*(S_1 + T_1)
			T <- coefqk * (S-S_1) + (1-coefqk) * (T_1)
			if (opterrType=="MAE"){
				E_ErrType <- abs((S + T) - X[i+1])
			}else if (opterrType=="MAPE"){
				E_ErrType <- abs(((S + T) - X[i+1])/X[i+1])*100
			}else if (opterrType=="sMAPE"){
				E_ErrType <- (abs((S + T) - X[i+1])/(abs(S + T) + abs(X[i+1])))*200
			}else if (opterrType=="MSE"){
				E_ErrType <- ((S + T) - X[i+1])^2
			}else {
			}
			tempCalc <- matrix(c(0, qk, pk, i, coefqk, coefpk, X[i], (S + T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
			eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
			eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))	
		}else {
			if (mdlType=="M"){
				coefpk <- abs(pk/i)
				coefqk <- abs(qk/i)
				ATA <- ata.sub.func(X, pk, qk, i-1, SerieName, mdlType, opterrType, SeasActual)
				T_1 <- ATA$trend
				S_1 <- ATA$level
				S <- coefpk * X[i] + (1-coefpk)*(S_1 * T_1)
				T <- coefqk * (S/S_1) + (1-coefqk) * (T_1)
				if (opterrType=="MAE"){
					E_ErrType <- abs((S * T) - X[i+1])
				}else if (opterrType=="MAPE"){
					E_ErrType <- abs(((S * T) - X[i+1])/X[i+1])*100
				}else if (opterrType=="sMAPE"){
					E_ErrType <- (abs((S * T) - X[i+1])/(abs(S * T) + abs(X[i+1])))*200
				}else if (opterrType=="MSE"){
					E_ErrType <- ((S * T) - X[i+1])^2
				}else {
				}
				tempCalc <- matrix(c(0, qk, pk, i, coefqk, coefpk, X[i], (S * T), seasIndex[i+1], S, T, E_ErrType), nrow=1, ncol=12, byrow = TRUE)
				eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
				eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
			}
		}
		my_list <- list("level" = S, "trend" = T)
		return(my_list)
	}else {
		tempCalc <- matrix(c(0, qk, pk, i, NA, NA, X[i], NA, seasIndex[i+1], NA, NA, NA), nrow=1, ncol=12, byrow = TRUE)
		eval(parse(text= paste("colnames(tempCalc) <- c(\"Sample\",\"q\", \"p\", \"t\", \"coefq\", \"coefp\", \"Actual\", \"Fit\", \"SeasComp\", \"Level\", \"Trend\",\"", as.character(opterrType),"\")", sep="")))
		eval(parse(text= paste(varATA, " <<- rbind(",varATA,",tempCalc)",sep="")))
		my_list <- list("level" = NA, "trend" = NA)
		return(my_list)
	}
}