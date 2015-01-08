runWestAfricaAnalyses = function(){
    cat("This is going to take a *while*\n")
    result1 = WestAfricaAnalysisScript(0)
    cat("Zero degree of freedom analysis complete.\n")
    result2 = WestAfricaAnalysisScript(3)
    cat("Three degree of freedom analysis complete.\n")
    result3 = WestAfricaAnalysisScript(5)
    cat("Five degree of freedom analysis complete.\n")
    return(list(zeroDF=result1,
                threeDF=result2,
                fiveDF=result3))
}
