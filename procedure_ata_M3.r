
require(forecast)
require(Mcomp)
require(tcltk2)
require(reshape2)

data(M3)
source("atamethod_M3.R")
source("atamethod_func_M3.R")
source("ata_sub_func_M3.R")
source("ataResultOutputCSV_M3.R")

atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "optPQ", errType = "sMAPE", fcast="single")
atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "optPQ", errType = "MAE", fcast="single")
atamethod(M3,sRange=1006:1007, parP="opt", parQ="opt", modelType = "A", optimPar = "optPQ", errType = "sMAPE", fcast="single")
atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "M", optimPar = "optPQ", errType = "sMAPE", fcast="single")

atamethod(M3,sRange=1:2, parP="opt", parQ=1:3, modelType = "A", optimPar = "optPQ", errType = "sMAPE", fcast="single")
atamethod(M3,sRange=1:2, parP=1:3, parQ="opt", modelType = "A", optimPar = "optPQ", errType = "sMAPE", fcast="single")
atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "pStar", errType = "sMAPE", fcast="single")
atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "pOne", errType = "sMAPE", fcast="single")
atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "pStarQ", errType = "sMAPE", fcast="single")

atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "both", errType = "sMAPE", fcast="single")
atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "optPQ", errType = "sMAPE", fcast="comb")
atamethod(M3,sRange=1:2, parP="opt", parQ="opt", modelType = "A", optimPar = "optPQ", errType = "sMAPE", fcast="select")


atamethod(M3,sRange=880, parP=63, parQ=1, modelType = "A", optimPar = "optPQ", errType = "sMAPE", fcast="single")