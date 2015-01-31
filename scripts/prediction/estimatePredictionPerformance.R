library(EARNMC)
library(parallel)
# cluster of size 2, each cluster makes cluster of size 3, maxes out 8 cores
cl = makeCluster(2)

runFunc = function(params){
    DF = params[[1]]
    tpts = params[[2]]
    result = WestAfricaAnalysisScript(DF, throwAwayTpts=tpts)
    save("result", file=paste("./DF",DF, "/pred_performance_dropping_", tpts, ".Rda.bz2", sep =""), compress="bzip2")
}

DFs = c(0,1,2,3)
tpts = c(20,15,10,5)
params = expand.grid(DFs, tpts)
paramsList = lapply(1:nrow(params), function(x){params[x,]})
parLapply(cl, paramsList, runFunc)
