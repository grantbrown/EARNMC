uncumulate = function(x)
{
  out = c(x[2:length(x)]-x[1:(length(x)-1)])
  ifelse(out >= 0, out, 0)
}

processGraphDataWA = function(filename=NA){
  pred.days = 120
  targetDaysPerRecord = 7
  if (all(is.na(filename)))
  {
      print("Using included data")
      data(WestAfricaEbola)
  }
  else{
    WestAfricaEbola = read.csv(filename)
  }

  dat = WestAfricaEbola[,c(1, 3,4,5)]
  
  charDate = as.character(dat[,1])
  
  rptDate = as.Date(charDate, "%m/%d/%Y")
  numDays = max(rptDate) - min(rptDate)
  numDays.pred = numDays + pred.days
  
  original.rptDate = rptDate
  ascendingOrder = order(rptDate)
  rptDate = rptDate[ascendingOrder]
  original.rptDate = original.rptDate[ascendingOrder]
  
  
  cleanData = function(dataColumn, ascendingOrder)
  {
    # Remove commas
    charCol = as.character(dataColumn)[ascendingOrder]
    
    if (is.na(charCol[1]))
    {
      charCol[1] = "0"
    }
    charCol = as.numeric(charCol)
    
    for (i in 2:length(charCol))
    {
      if (is.na(charCol[i]))
      {
        charCol[i] = charCol[i-1]  
      }
    }
    # charCol
    # Correct for undercounts
    for (i in seq(length(charCol), 2))
    {
      if (charCol[i-1] > charCol[i])
      {
        charCol[i-1] = charCol[i]
      }
    }
    charCol
  }
  
  Guinea = cleanData(dat[,2], ascendingOrder)
  Liberia = cleanData(dat[,3], ascendingOrder)
  Sierra.Leone = cleanData(dat[,4], ascendingOrder)
  rawData = cbind(Guinea, Sierra.Leone, Liberia)
  rownames(rawData) = as.character(original.rptDate)
  colnames(rawData) = paste(paste("&nbsp;&nbsp;", c("Guinea", "Liberia", "Sierra Leone")), "&nbsp;&nbsp;")
  
  # The data needs to be aggregated: there's some error in the measurements, 
  # and we're not actually observing infection times. The data is therefore
  # recorded at an artifically high time scale.
  
  nDays = uncumulate(original.rptDate)
  
  thinIndices = function(minDays, weights)
  {
    keepIdx = c(length(weights))
    currentWeight = 0
    lastIdx = -1
    for (i in seq(length(weights)-1, 1))
    {
      currentWeight = currentWeight + weights[i]
      if (currentWeight >= minDays)
      {
        currentWeight = 0
        keepIdx = c(keepIdx, i)
        lastIdx = i
      }
    }
    if (currentWeight != 0)
    {
      keepIdx = c(keepIdx, lastIdx-1)
    }
    keepIdx
  }
  
  keepIdx = thinIndices(targetDaysPerRecord, c(nDays,1))
  keepIdx = keepIdx[order(keepIdx)]
  
  return(list(Guinea=Guinea,
              Sierra.Leone=Sierra.Leone,
              Liberia=Liberia,
              original.rptDate=original.rptDate,
              rptDate=rptDate,
              keepIdx=keepIdx,
              rawData=rawData))
}

westAfricaRawDataPlot = function(){
  
  processedData = processGraphDataWA()
  Guinea = processedData$Guinea[processedData$keepIdx]
  Sierra.Leone = processedData$Sierra.Leone[processedData$keepIdx]
  Liberia = processedData$Liberia[processedData$keepIdx]
  
  original.rptDate = processedData$original.rptDate[processedData$keepIdx]
  rptDate = original.rptDate[2:length(original.rptDate)]
  
  
  I_star = cbind(uncumulate(Guinea), 
                 uncumulate(Liberia), 
                 uncumulate(Sierra.Leone))
  I0 = c(Guinea[1], Liberia[1], Sierra.Leone[1])
  
  offsets = uncumulate(original.rptDate)
  
  
  ## Do a ggplot2 summary
  plotData = data.frame(country=rep(c("Guinea", "Liberia", "Sierra Leone"), each = length(uncumulate(Guinea))),
                   date=rep(rptDate, 3),cases=c(uncumulate(Guinea), uncumulate(Liberia), uncumulate(Sierra.Leone)),
                   d = c(uncumulate(original.rptDate)))
  
  cases.plot = qplot(x=date, y=cases, fill=country,
                     data=plotData, geom="bar", stat="identity",
                      position="dodge",xlab="Date", ylab ="Cases")
  cases.plot.greyscale = cases.plot + scale_fill_grey(start = 0, end = .9)
  cases.plot.greyscale = cases.plot.greyscale + theme_minimal()
  
  list(colorPlot=cases.plot,
       greyscalePlot=cases.plot.greyscale)
}

generateWestAfricaRawDataPlotPDF = function(greyscaleFileName="./InfectionsBarplotBW.pdf",
                                    colorFileName = "./InfectionsBarplotColor.pdf"){
  plots = westAfricaRawDataPlot()
  plot1 = plots[[1]];
  plot2 =  plots[[2]];
  pdf(file=colorFileName, width=8, height = 4)
    print(plot1)
  dev.off()
  
  pdf(file=greyscaleFileName, width=8, height = 4)
    print(plot2) 
  dev.off()
}


makeR0ComparisonPlot = function(outputFile, ylim=c(0,3)){
  ModelOutput = outputFile
  processedData = processGraphDataWA()
  
  original.rptDate = processedData$original.rptDate[processedData$keepIdx]
  rptDate = original.rptDate[2:length(original.rptDate)]
  
  nTpt = nrow(ModelOutput[[2]]$R0$empiricalR0$mean) -1
  EA_R = ModelOutput[[2]]$R0$empiricalR0
  R0 = ModelOutput[[2]]$R0$R0
  
  
  color1 = rgb(0,0,0,0.9)
  color2 = rgb(0,0,1,0.9)
  lwdMean = 3
  lwdCI = 1
  ltyCI=2
  lty1 = 1
  lty2 = 4
  pch1=1
  pch2=3
  ltyThreshold = 5
  pointCex = 0
  colThreshold = rgb(0.5,0.5,0.5,0.5)
  fig1Ylim=c(0,2)
  
  
  R0ComparisonPlot = function(EARList,R0List,colIdx,main, ...)
  {    
    EARMean = EARList$mean[2:nTpt,colIdx]
    EARUB = EARList$UB[2:nTpt,colIdx]
    EARLB = EARList$LB[2:nTpt,colIdx]
    
    R0Mean = R0List$mean[2:nTpt,colIdx]
    R0UB = R0List$UB[2:nTpt,colIdx]
    R0LB = R0List$LB[2:nTpt,colIdx]
    
    plot(rptDate[2:nTpt],EARMean,
         type = "n", ...)
    title(main=main, line=-2.5)
    abline(h = 1, lty = ltyThreshold, col = colThreshold)
    lines(rptDate[2:nTpt], EARMean, col = color1, lwd = lwdMean, lty = lty1,)
    
    points(rptDate[2:nTpt], EARMean, col = color1, pch=pch1, cex = pointCex)
    lines(rptDate[2:nTpt],EARUB, col = color1, lwd = lwdCI, lty = ltyCI)
    lines(rptDate[2:nTpt],EARLB, col = color1, lwd = lwdCI, lty = ltyCI)
    
    lines(rptDate[2:nTpt],R0Mean,  col = color2, lwd = lwdMean, lty = lty2)
    points(rptDate[2:nTpt], R0Mean, col = color2, pch = pch2, cex = pointCex)
    lines(rptDate[2:nTpt],R0UB, col = color2, lwd = lwdCI, lty = ltyCI)
    lines(rptDate[2:nTpt],R0LB, col = color2, lwd = lwdCI, lty = ltyCI)
    
  }
  
  figureFunc = function(EARList, R0List, titleVal, ylim)
  {
    par("oma" = c(1,1,1,1))
    layout(matrix(c(1,2,3,
                    4,4,4), 2, 3, byrow = TRUE),
           widths=c(3,3,3),
           heights = c(5,1))
    
    R0ComparisonPlot(
      EARList,
      R0List,
      1,
      main = "Guinea",
      xlab = "Date",
      ylab = "Reproductive Number",
      ylim=ylim)
    
    R0ComparisonPlot( EARList,
                      R0List,
                      2,
                      main = "Liberia",
                      xlab = "Date",
                      ylab = "Reproductive Number",
                      ylim=ylim)
    
    R0ComparisonPlot( EARList,
                      R0List,
                      3,
                      main = "Sierra Leone",
                      xlab = "Date",
                      ylab = "Reproductive Number",
                      ylim = ylim)
    par(xaxt="n")
    par(yaxt="n")
    par(xpd=TRUE)
    par(bty="n")
    plot(c(0,10),c(0,10), type = "n",xlab="",ylab="")
    legend(x = 5, 
           y = 5, 
           lty = c(lty1, lty2), 
           col = c(color1, color2), 
           legend=c(expression('R'^{(EA)}), expression('R'[0](t))),
           lwd = c(lwdMean, lwdMean), 
           pt.cex = c(pointCex,pointCex),
           pt.lwd= c(pointCex, pointCex),
           pch = c(pch1, pch2),
           cex=3,
           xjust=0.5,
           yjust=0.5, 
           horiz = TRUE,
           text.width=c(1))
    par(xaxt="s")
    par(yaxt="s")
    par(xpd=FALSE)
    par(bty="o")
    par(cex.main=1.5)
    title(titleVal, outer = TRUE, line = -2)
    par(cex.main=1)
  } 
  
  
  figureFunc(EA_R, R0, "", ylim)
}



makeR0ComparisonPlotPDFs = function(zeroDFOutputFile = NA,
                                    threeDFOutputFile = NA,
                                    fiveDFOutputFile = NA){
  if (is.na(zeroDFOutputFile)){
    print("Zero degree of freedom results file not provided, using included version.")
    data(WestAfrica0DFResults)
    zeroDFOutputFile = WestAfrica0DFResults
  }
  if (is.na(threeDFOutputFile)){
    print("Three degree of freedom results file not provided, using included version.")
    data(WestAfrica3DFResults)
    threeDFOutputFile = WestAfrica3DFResults
  }
  if (is.na(fiveDFOutputFile)){
    print("Five degree of freedom results file not provided, using included version.")
    data(WestAfrica5DFResults)
    fiveDFOutputFile = WestAfrica5DFResults
  }
  
  pdf(file="./0DF_Pairwise.pdf", width = 12, height =8)
    makeR0ComparisonPlot(zeroDFOutputFile, c(0, 1.8))
  dev.off()
  
  pdf(file="./3DF_Pairwise.pdf", width = 12, height =8)
    makeR0ComparisonPlot(threeDFOutputFile, c(0, 3.2))
  dev.off()
  
  pdf(file="./5DF_Pairwise.pdf", width = 12, height =8)
    makeR0ComparisonPlot(fiveDFOutputFile, c(0, 3.5))
  dev.off()
  
}

makeWestAfricaPredictionPlot = function(ModelResults=NA)
{
  if (is.na(ModelResults)){
    print("No model results file provided, using included data.")
    data(WestAfrica3DFResults)
    ModelResults = WestAfrica3DFResults
  }
  pred.days=120
  processedData = processGraphDataWA()
  Guinea = processedData$Guinea[processedData$keepIdx]
  Sierra.Leone = processedData$Sierra.Leone[processedData$keepIdx]
  Liberia = processedData$Liberia[processedData$keepIdx]
  
  original.rptDate = processedData$original.rptDate[processedData$keepIdx]
  rptDate = original.rptDate[2:length(original.rptDate)]
  
  
  I_star = cbind(uncumulate(Guinea), 
                 uncumulate(Liberia), 
                 uncumulate(Sierra.Leone))
  I0 = c(Guinea[1], Liberia[1], Sierra.Leone[1])
  
  offsets = uncumulate(original.rptDate)
  
  
  
  modelDF=3
  N = matrix(c(10057975, 4128572, 6190280), nrow = nrow(I_star),ncol = 3, 
             byrow=TRUE)
  X = diag(ncol(N))
  X.predict = X
  
  daysSinceJan = as.numeric(rptDate - as.Date("2014-01-01"))
  daysSinceJan.predict = c(max(daysSinceJan) + 1, max(daysSinceJan) 
                           + seq(2,pred.days-2,2))
  
  splineBasis = ns(daysSinceJan, df = modelDF)
  splineBasis.predict = predict(splineBasis, daysSinceJan.predict)
  Z = matrix(splineBasis, ncol = modelDF)
  Z.predict = splineBasis.predict
  
  
  # These co-variates are the same for each spatial location, 
  # so duplicate them row-wise. 
  Z = Z[rep(1:nrow(Z), nrow(X)),,drop=FALSE]
  Z.predict = Z.predict[rep(1:nrow(Z.predict), nrow(X)),,drop=FALSE]
  
  # For convenience, let's combine X and Z for prediction.
  X.pred = cbind(X.predict[rep(1:nrow(X.predict), 
                               each = nrow(Z.predict)/nrow(X)),], Z.predict)
  
  DM1 = matrix(c(0,1,0,
                 1,0,0,
                 0,0,0), nrow = 3, ncol = 3, byrow = TRUE)
  DM2 = matrix(c(0,0,1,
                 0,0,0,
                 1,0,0), nrow = 3, ncol = 3, byrow = TRUE)
  DM3 = matrix(c(0,0,0,
                 0,0,1,
                 0,1,0), nrow = 3, ncol = 3, byrow = TRUE)
  dmList = list(DM1, DM2, DM3)
  
  
  # Define prediction offsets. 
  offset.pred = uncumulate(c(1,seq(2,pred.days-2,2)))
  
  
  chain1 = ModelResults[[1]]$chainOutput
  chain2 = ModelResults[[2]]$chainOutput
  chain3 = ModelResults[[3]]$chainOutput
  
  nbeta = 6
  c1 = chain1[floor(nrow(chain1)/2):nrow(chain1),c(1:(nbeta+length(dmList)+3))]
  c2 = chain2[floor(nrow(chain2)/2):nrow(chain2),c(1:(nbeta+length(dmList)+3))]
  c3 = chain3[floor(nrow(chain3)/2):nrow(chain3),c(1:(nbeta+length(dmList)+3))]
  
  c1$p_ei = 1-exp(-c1$gamma_ei)
  c1$p_ir = 1-exp(-c1$gamma_ir)
  c2$p_ei = 1-exp(-c2$gamma_ei)
  c2$p_ir = 1-exp(-c2$gamma_ir)
  c3$p_ei = 1-exp(-c3$gamma_ei)
  c3$p_ir = 1-exp(-c3$gamma_ir)
  dayFunc = function(x){rgeom(rep(1, length(x)), x)}
  c1$exposedDays = dayFunc(c1$p_ei)
  c1$infectiousDays = dayFunc(c1$p_ir)
  c2$exposedDays = dayFunc(c2$p_ei)
  c2$infectiousDays = dayFunc(c2$p_ir)
  c3$exposedDays = dayFunc(c3$p_ei)
  c3$infectiousDays = dayFunc(c3$p_ir)
  
  colnames(c1) = c("Guinea Intercept", "Liberia Intercept", 
                   "Sierra Leone Intercept", 
                   paste("Time component ", (1:(nbeta-3)), sep = ""),
                   paste("Spatial Dependence Parameter", 1:length(dmList)),
                   "Overdispersion Precision",
                   "E to I Transition Parameter", 
                   "I to R Transition Parameter",
                   "E to I Transition Probability", 
                   "I to R Transition Probability",
                   "Days in Exposed Category",
                   "Days in Infectious Category")
  colnames(c2) = colnames(c1)
  colnames(c3) = colnames(c1)
  mcl = mcmc.list(as.mcmc(c1), 
                  as.mcmc(c2),
                  as.mcmc(c3))
  summaryTable = summary(mcl)
  summaryTableLatex = kable(round(summaryTable$quantiles,2), format = "latex")
  cat(summaryTableLatex, file = "./DF3QuantlesTable.partial.tex")
  
  getMeanAndCI = function(loc,tpt,baseStr="I_")
  {
    vec = chain1[[paste(baseStr, loc, "_", tpt, sep = "")]]
    vec = vec[floor(length(vec)/2):length(vec)]
    return(c(mean(vec), quantile(vec, probs = c(0.05, 0.95))))
  }
  
  Guinea.I.Est = sapply(0:(nrow(I_star)- 1), getMeanAndCI, loc=0)
  Liberia.I.Est = sapply(0:(nrow(I_star)- 1), getMeanAndCI, loc=1)
  SierraLeone.I.Est = sapply(0:(nrow(I_star)- 1), getMeanAndCI, loc=2)
  maxIdx = nrow(I_star)
  # Declare prediction functions
  predictEpidemic = function(beta.pred, 
                             X.pred,
                             gamma.ei,
                             gamma.ir,
                             S0,
                             E0,
                             I0,
                             R0,
                             rho,
                             offsets.pred)
  {
    N = (S0+E0+I0+R0)
    offsets.pred=c(offsets.pred,1)
    p_se_components = matrix(exp(X.pred %*% beta.pred), ncol=length(S0))
    p_se = matrix(0, ncol = length(S0), nrow = nrow(p_se_components))
    p_ei = 1-exp(-gamma.ei*offsets.pred)
    p_ir = 1-exp(-gamma.ir*offsets.pred)
    S_star = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    E_star = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    I_star = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    R_star = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    S = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    E = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    I = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    R = matrix(0, ncol=length(S0),nrow = nrow(p_se_components))
    S[1,] = S0
    E[1,] = E0
    I[1,] = I0
    R[1,] = R0
    S_star[1,] = rbinom(rep(1, length(S0)), R0, 0)
    p_se[1,] = I[1,]/N*p_se_components[1,]
    for (i in 1:length(dmList))
    {
      p_se[1,] = p_se[1,] + rho[i]*(dmList[[i]] %*% (I[1,]/N*p_se_components[1,]))
    }
    p_se[1,] = 1-exp(-offsets.pred[1]*(p_se[1,]))
    
    E_star[1,] = rbinom(rep(1, length(S0)), S0, p_se[1,])
    I_star[1,] = rbinom(rep(1, length(S0)), E0, p_ei[1])
    R_star[1,] = rbinom(rep(1, length(S0)), I0, p_ir[1])
    
    for (i in 2:nrow(S))
    {
      
      S[i,] = S[i-1,] + S_star[i-1,] - E_star[i-1,]
      E[i,] = E[i-1,] + E_star[i-1,] - I_star[i-1,]
      I[i,] = I[i-1,] + I_star[i-1,] - R_star[i-1,]
      R[i,] = R[i-1,] + R_star[i-1,] - S_star[i-1,]
      
      p_se[i,] = I[i,]/N*p_se_components[i,]
      for (j in 1:length(dmList))
      {
        p_se[i,] = p_se[i,] + rho[j]*(dmList[[j]] %*% (I[i,]/N*p_se_components[i,]))
      }
      p_se[i,] = 1-exp(-offsets.pred[i]*(p_se[i,]))
      
      
      S_star[i,] = rbinom(rep(1, length(S0)), R[i,], 0)
      E_star[i,] = rbinom(rep(1, length(S0)), S[i,], p_se[i,])
      I_star[i,] = rbinom(rep(1, length(S0)), E[i,], p_ei[i])
      R_star[i,] = rbinom(rep(1, length(S0)), I[i,], p_ir[i])
    }
    return(list(S=S,E=E,I=I,R=R,
                S_star=S_star,E_star=E_star,
                I_star=I_star,R_star=R_star,
                p_se=p_se,p_ei=p_ei,p_ir=p_ir))
  }
  
  predict.i = function(i)
  {
    dataRow = currentChain[i,]
    rho = rep(0, length(dmList))
    for (i in 1:length(dmList))
    {
      rho[i] = dataRow[[paste("rho_", i-1, sep = "")]]
    }
    beta = rep(0, modelDF+ncol(X))
    for (i in 0:(modelDF+ncol(X) -1))
    {
      beta[i+1] = dataRow[[paste("BetaP_SE_", i, sep = "")]]
    }
    
    S0 = c(dataRow[[paste("S_0_", maxIdx-1, sep = "")]],
           dataRow[[paste("S_1_", maxIdx-1, sep = "")]],
           dataRow[[paste("S_2_", maxIdx-1, sep = "")]],
           dataRow[[paste("S_3_", maxIdx-1, sep = "")]])
    E0 = c(dataRow[[paste("E_0_", maxIdx-1, sep = "")]],
           dataRow[[paste("E_1_", maxIdx-1, sep = "")]],
           dataRow[[paste("E_2_", maxIdx-1, sep = "")]],
           dataRow[[paste("E_3_", maxIdx-1, sep = "")]])
    I0 = c(dataRow[[paste("I_0_", maxIdx-1, sep = "")]],
           dataRow[[paste("I_1_", maxIdx-1, sep = "")]],
           dataRow[[paste("I_2_", maxIdx-1, sep = "")]],
           dataRow[[paste("I_3_", maxIdx-1, sep = "")]])
    R0 = c(dataRow[[paste("R_0_", maxIdx-1, sep = "")]],
           dataRow[[paste("R_1_", maxIdx-1, sep = "")]],
           dataRow[[paste("R_2_", maxIdx-1, sep = "")]],
           dataRow[[paste("R_3_", maxIdx-1, sep = "")]])
    
    
    return(predictEpidemic(beta,  
                           X.pred,
                           dataRow$gamma_ei,
                           dataRow$gamma_ir,
                           S0,
                           E0,
                           I0,
                           R0,
                           rho,
                           offset.pred
    ))
  }
  
  
  
  currentChain = chain1
  preds = lapply(rep(seq(floor(nrow(chain1)/4),
                         nrow(chain1), 5), each = 5), predict.i)
  
  pred.dates = c(rptDate[(which.max(rptDate))],
                 rptDate[(which.max(rptDate))] + seq(2,pred.days-2,2))
  pred.xlim = c(min(rptDate), max(pred.dates))
  lastIdx = nrow(I_star)
  Guinea.Pred = preds[[1]]$I[,1]
  Liberia.Pred = preds[[1]]$I[,2]
  SierraLeone.Pred = preds[[1]]$I[,3]
  
  
  breakpoint = mean(c(max(rptDate), min(pred.dates)))
  
  for (predIdx in 2:length(preds))
  {
    Guinea.Pred = rbind(Guinea.Pred, preds[[predIdx]]$I[,1])
    Liberia.Pred = rbind(Liberia.Pred, preds[[predIdx]]$I[,2])
    SierraLeone.Pred = rbind(SierraLeone.Pred, preds[[predIdx]]$I[,3])
  }
  
  Guinea.mean = apply(Guinea.Pred, 2, mean)
  Liberia.mean = apply(Liberia.Pred, 2, mean)
  SierraLeone.mean = apply(SierraLeone.Pred, 2, mean)
  
  
  Guinea.LB = apply(Guinea.Pred, 2, quantile, probs = c(0.025))
  Guinea.UB = apply(Guinea.Pred, 2, quantile, probs = c(0.975))
  
  Liberia.LB = apply(Liberia.Pred, 2, quantile, probs = c(0.025))
  Liberia.UB = apply(Liberia.Pred, 2, quantile, probs = c(0.975))
  
  SierraLeone.LB = apply(SierraLeone.Pred, 2, quantile, probs = c(0.025))
  SierraLeone.UB = apply(SierraLeone.Pred, 2, quantile, probs = c(0.975))
  
  maxI = max(c(max(c(Guinea.I.Est, Liberia.I.Est, SierraLeone.I.Est)), Guinea.UB, Liberia.UB, SierraLeone.UB))
  
  est.idx = seq(1, length(Guinea.I.Est[1,]), 2)
  pred.table1 = cbind(Guinea.I.Est[1,], 
                      Liberia.I.Est[1,],
                      SierraLeone.I.Est[1,]
  )[est.idx,]
  pred.table2 = cbind(Guinea.mean, 
                      Liberia.mean,
                      SierraLeone.mean)
  pred.table = rbind(pred.table1, pred.table2)
  rownames(pred.table) = paste("", c(as.character(rptDate)[est.idx], as.character(pred.dates)),
                               sep = "")
  rownames(pred.table) = paste(rownames(pred.table), "&nbsp;", sep = "")
  colnames(pred.table) = c("Guinea", 
                           "Liberia", 
                           "Sierra Leone")
  
  
  
  
  x = rptDate - min(rptDate)
  guinea.interp = approx(x,Guinea.I.Est[1,],xout = 0:max(x))
  liberia.interp = approx(x,Liberia.I.Est[1,],xout = 0:max(x))
  sierraleone.interp = approx(x,SierraLeone.I.Est[1,],xout = 0:max(x))
  interpDate = min(rptDate):max(rptDate)
  
  nobs =  length(c(sierraleone.interp$y, Guinea.mean))
  predPlotData = data.frame(InfectiousSize = c(guinea.interp$y, Guinea.mean,
                                               liberia.interp$y, Liberia.mean,
                                               sierraleone.interp$y, SierraLeone.mean),
                            Date = as.Date(rep(c(interpDate, pred.dates), 3), origin="1970-01-01"),
                            country = rep(c("Guinea", "Liberia", "Sierra Leone"), each = nobs))
  
  
  
  infections.plot = qplot(x=Date, y=InfectiousSize, color=country,
                          data=predPlotData, geom="line", stat="identity",
                          xlab="Date", ylab ="Number Infectious")
  infections.plot = infections.plot + geom_line(size=1)
  infections.plot = infections.plot + geom_hline(yintercept=0, size = 0.2)
  infections.plot = infections.plot + geom_vline(xintercept= as.numeric(min(pred.dates)), size = 0.2, lty=2)
  infections.plot = infections.plot + annotate("text", x = min(pred.dates), -90, label="Prediction", 
                                               color="blue", hjust=-0.2)
  infections.plot = infections.plot + annotate("text", x = min(pred.dates), -90, label="Estimation", 
                                               color="blue", hjust=1.2)
  
  infections.plot.greyscale = infections.plot + scale_color_grey(start = 0, end = .9)
  infections.plot.greyscale = infections.plot.greyscale  + theme_bw()

  return(list(colorPlot = infections.plot, 
              greyscalePlot = infections.plot.greyscale))
}
# Make a spline basis plot


makeSplineBasisPlotPDF = function(){
  pred.days = 120
  pdf(file="./splineBasisDF3.pdf", width = 8, height = 4);
  
    processedData = processGraphDataWA()
    Guinea = processedData$Guinea[processedData$keepIdx]
    Sierra.Leone = processedData$Sierra.Leone[processedData$keepIdx]
    Liberia = processedData$Liberia[processedData$keepIdx]
    
    original.rptDate = processedData$original.rptDate[processedData$keepIdx]
    rptDate = original.rptDate[2:length(original.rptDate)]
    
    
    I_star = cbind(uncumulate(Guinea), 
                   uncumulate(Liberia), 
                   uncumulate(Sierra.Leone))
    I0 = c(Guinea[1], Liberia[1], Sierra.Leone[1])
    
    offsets = uncumulate(original.rptDate)
  
    modelDF=3
    N = matrix(c(10057975, 4128572, 6190280), nrow = nrow(I_star),ncol = 3, 
               byrow=TRUE)
    X = diag(ncol(N))
    X.predict = X
    
    daysSinceJan = as.numeric(rptDate - as.Date("2014-01-01"))
    daysSinceJan.predict = c(max(daysSinceJan) + 1, max(daysSinceJan) 
                             + seq(2,pred.days-2,2))
    
    splineBasis = ns(daysSinceJan, df = modelDF)
    splineBasis.predict = predict(splineBasis, daysSinceJan.predict)
    Z = matrix(splineBasis, ncol = modelDF)
    Z.predict = splineBasis.predict
    
  
    combinedIdx = c(daysSinceJan, daysSinceJan.predict) - min(daysSinceJan)
    combinedSplines = rbind(splineBasis, splineBasis.predict)
    par(mfrow=c(1,3))
    plot(combinedIdx, combinedSplines[,1], xlab = "Epidemic Days", 
         ylab = "Basis Parameter Value",
         main = "Temporal Basis Component 1", type = "l",
         lwd=2)
    abline(v = min(daysSinceJan.predict- min(daysSinceJan)), col = "lightgrey", lty=2)
    plot(combinedIdx, combinedSplines[,2], xlab = "Epidemic Days", 
         ylab = "Basis Parameter Value",
         main = "Temporal Basis Component 2", type = "l",
         lwd=2)
    abline(v = min(daysSinceJan.predict - min(daysSinceJan)), col = "lightgrey", lty=2)
    plot(combinedIdx, combinedSplines[,3], xlab = "Epidemic Days", 
         ylab = "Basis Parameter Value",
         main = "Temporal Basis Component 3", type = "l",
         lwd=2)
    abline(v = min(daysSinceJan.predict - min(daysSinceJan)), col = "lightgrey", lty=2)
  dev.off()
}

makeWestAfricaPredictionPlotPDFs = function(ModelResults=NA){
  plots = makeWestAfricaPredictionPlot(ModelResults)
  pdf(file="WestAfricaPrediction.pdf", width = 8, height =4)
    print(plots[[1]])
  dev.off()
  
  pdf(file="WestAfricaPredictionGreyscale.pdf", width = 8, height =4)
    print(plots[[2]])
  dev.off()
}





