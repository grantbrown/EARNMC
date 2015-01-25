library(EARNMC)

for (tpts in c(20,15,10,5)){
    result = WestAfricaAnalysisScript(3, throwAwayTpts=tpts)
    save("result", file=paste("./pred_performance_dropping_", result, ".Rda.bz2", sep =""), compress="bzip2")
}
