library(dplyr)
library(EARNMC)
library(parallel)
library(Matrix)
library(ggplot2)

trueEffectiveR0 = function(infectious, eventArray)
{
  EA_R0 = numeric(nrow(infectious))
  for (i in 1:nrow(infectious))
  {
    infected = which(infectious[i,] == 1)
    if (length(infected) != 0)
    {
      EA_R0[i] = sum(sapply(infected, function(x){sum(eventArray[x,,i:nrow(infectious)])}))
    }
  }
  Ipop = apply(infectious, 1, sum)
  ksmooth(1:length(EA_R0), EA_R0/ifelse(Ipop != 0, Ipop, 1), bandwidth = length(EA_R0)/10)$y
}


trueEffectiveR02 = function(infectious, eventArray)
{
  EA_R0 = numeric(nrow(infectious))
  for (i in 1:nrow(infectious))
  {
    infected = which(infectious[i,] == 1)
    if (length(infected) != 0)
    {
      EA_R0[i] = sum(sapply(infected, function(x){sum(eventArray[x,,i:nrow(infectious)])}))
    }
  }
  Ipop = apply(infectious, 1, sum)
  EA_R0/ifelse(Ipop != 0, Ipop, 1)
}


estimateR0 = function(dataName)
{
  load(dataName)
  nSamples = nrow(result$dat1)
  out = lapply((floor(nSamples/2)):nSamples, 
               function(x){
                 r0.i(x,result$dat1, result$simResults$X)
               })


  trueER = matrix(trueEffectiveR02(result$simResults$I, result$simResults$infectionEvents), 
                  nrow = nrow(result$simResults$E_star),
                  ncol=(nSamples - floor(nSamples/2)) + 1)
  EARMatrix = matrix(unlist(
                    lapply(out, function(x){x$EAR})), 
                    nrow = result$simResults$nTpt)
  
  R0tMatrix = matrix(unlist(
    lapply(out, function(x){x$R0t})), 
    nrow = result$simResults$nTpt)
  
  EffR0tMatrix = matrix(unlist(
    lapply(out, function(x){x$effR0t})), 
    nrow = result$simResults$nTpt)
  
  sm.fnc = function(x)
  {
    rbind(apply(x, 1, quantile, probs = c(0.05, 0.5, 0.95)),
          apply(x, 1, mean))
  }
  EARQuantiles = sm.fnc(EARMatrix)
  R0tQuantiles = sm.fnc(R0tMatrix)
  EffR0tQuantiles = sm.fnc(EffR0tMatrix)
  EARDifferenceQuantiles =  sm.fnc(EARMatrix-trueER)
  R0tDifferenceQuantiles =  sm.fnc(R0tMatrix-trueER)
  EffR0tDifferenceQuantiles = sm.fnc(EffR0tMatrix-trueER)
  
  
  r0.concordance =  mean((R0tMatrix >= 1 & trueER >= 1) | (R0tMatrix < 1 & trueER < 1))
  EAR.concordance =  mean((EARMatrix >= 1 & trueER >= 1) | (EARMatrix < 1 & trueER < 1))
  EffR0.concordance =  mean((EffR0tMatrix >= 1 & trueER >= 1) | (EffR0tMatrix < 1 & trueER < 1))
                     
  list(EARQuantiles=EARQuantiles,
       R0tQuantiles=R0tQuantiles,
       EffR0tQuantiles=EffR0tQuantiles,
       R0tDifferenceQuantiles=R0tDifferenceQuantiles,
       EARDifferenceQuantiles=EARDifferenceQuantiles,
       EffR0tDifferenceQuantiles=EffR0tDifferenceQuantiles,
       r0.concordance=r0.concordance,
       EAR.concordance=EAR.concordance,
       EffR0.concordance=EffR0.concordance,
       EARvsR0 = mean(abs(EARMatrix-trueER) < abs(R0tMatrix-trueER)),
       EARvsEffR0 = mean(abs(EARMatrix-trueER) < abs(EffR0tMatrix-trueER))
       )
}

r0.i = function(i, chainObject, X_se){
  startPt=0
  dataRow = chainObject[i,]
  S = as.numeric(select(dataRow, starts_with("S_0_")))
  E = as.numeric(select(dataRow, starts_with("E_0_"))) 
  I = as.numeric(select(dataRow, starts_with("I_0_"))) 
  R = as.numeric(select(dataRow, starts_with("R_0_"))) 
  N = S[1] + E[1] + I[1] + R[1]
  nTpt = nrow(X_se)/length(N)

  beta.SE = as.numeric(dataRow[grep("BetaP_SE_", 
                                    colnames(dataRow))])
  if (length(beta.SE) != ncol(X_se))
  {
    X_se = X_se[,1:length(beta.SE), drop=FALSE]
  }
  eta.SE = X_se %*% beta.SE
  gamma_ei = dataRow$gamma_ei
  gamma_ir = dataRow$gamma_ir
  p_se = 1-exp(-I/N*exp(eta.SE))
  p_ir = 1-exp(-gamma_ir)
  EA_R0 = rep(0, nTpt)
  R0_t = exp(eta.SE)/gamma_ir
  Eff_R0 = R0_t*(S/N)
  
  tmp = ifelse(I != 0, S/I*p_se, 0)
  for (t in 1:nTpt)
  {
    cpir = 1
    idx = 0
    if (I[t] != 0){
      while(cpir > 1e-16 && idx <= nTpt - t)
      {
        EA_R0[t] = EA_R0[t] + tmp[t+idx]*cpir
        cpir = cpir * (1-p_ir)
        idx = idx + 1
      }
      if (idx >= nTpt){
        while(cpir > 1e-8){
          EA_R0[t] = EA_R0[t] + tmp[nTpt]*cpir
          cpir = cpir * (1-p_ir[nTpt])
        }
      }
    }
  }
  list(EAR=matrix(EA_R0, ncol = 1),
       R0t=matrix(R0_t, ncol = 1), 
       effR0t = matrix(Eff_R0, ncol = 1))
}

#cl = makeCluster(8, outfile="./err.txt")
#parLapply(cl, 1:8, function(x){library(dplyr); library(EARNMC)})
#clusterExport(cl, c("estimateR0", "r0.i", "trueEffectiveR0",  "trueEffectiveR02"))
fnames = dir("./results")
fnames = fnames[fnames != "README.txt"]
fnames.uspec = fnames[grep("TRUE", fnames)][1:100]
fnames.spec = fnames[grep("FALSE", fnames)][1:100]
fnames = c(fnames.uspec, fnames.spec)
#R0Results = parLapplyLB(cl, fnames, function(x){cat(paste(x, "\n", sep = ""));
#                                                  estimateR0(paste("./results/", x, sep = ""))})


R0Results =lapply(fnames, function(x){cat(paste(x, "\n", sep = ""));
                                                         estimateR0(paste("./results/", x, sep = ""))})

underspec.idx = grepl("Uspec_TRUE", fnames)
spec.idx = grepl("Uspec_FALSE", fnames)

EARQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                       function(x){R0Results[[x]][[1]]}))/sum(underspec.idx)
EARQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                       function(x){R0Results[[x]][[1]]}))/sum(spec.idx)

R0tQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                           function(x){R0Results[[x]][[2]]}))/sum(underspec.idx)
R0tQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                      function(x){R0Results[[x]][[2]]}))/sum(spec.idx)

EffR0tQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                           function(x){R0Results[[x]][[3]]}))/sum(underspec.idx)
EffR0tQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                      function(x){R0Results[[x]][[3]]}))/sum(spec.idx)

R0tDifferenceQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                           function(x){R0Results[[x]][[4]]}))/sum(underspec.idx)
R0tDifferenceQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                      function(x){R0Results[[x]][[4]]}))/sum(spec.idx)

EARDifferenceQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                                      function(x){R0Results[[x]][[5]]}))/sum(underspec.idx)
EARDifferenceQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                                 function(x){R0Results[[x]][[5]]}))/sum(spec.idx)


EffR0tDifferenceQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                                      function(x){R0Results[[x]][[6]]}))/sum(underspec.idx)
EffR0tDifferenceQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                                 function(x){R0Results[[x]][[6]]}))/sum(spec.idx)


save(file="./quantiles.Robj", "EARQuantiles.underspec",
                                "R0tQuantiles.underspec",
                                "EffR0tQuantiles.underspec",
                                "R0tDifferenceQuantiles.underspec",
                                "EARDifferenceQuantiles.underspec",
                                "EARQuantiles.spec",
                                "R0tQuantiles.spec",
                                "EffR0tQuantiles.spec",
                                "R0tDifferenceQuantiles.spec",
                                "EARDifferenceQuantiles.spec",
                                "EffR0tDifferenceQuantiles.spec",
                                "EffR0tDifferenceQuantiles.underspec",
                                "R0Results")




require(tikzDevice)

tikz(file = "./ABMBiasPlot.tex", width = 12, height = 8, standAlone=TRUE)
  showCI = FALSE
  showIdx = 2:ncol(EARDifferenceQuantiles.underspec)
  layout(matrix(c(1,2,
                  3,3), nrow = 2, ncol = 2, byrow = TRUE), 
         widths = c(8,8), heights = c(8,2)) 
  plot(EARDifferenceQuantiles.underspec[2,showIdx], type = "l", lwd = 3, ylim = c(-1,1),
       xlab = "Epidemic Time", ylab = "Bias", main = "Median Bias Under Misspecification")
  abline(h = 0, col = "blue", lty = 2, lwd = 2)
  abline(h = seq(-5, 5, 0.5), lty = 2, col = rgb(0,0,0,0.2))
  lines(EARDifferenceQuantiles.underspec[2,showIdx])
  lines(R0tDifferenceQuantiles.underspec[2,showIdx], col = "red", lwd = 3, lty = 3)
  lines(EffR0tDifferenceQuantiles.underspec[2,showIdx], col = "grey", lwd = 3, lty = 4)
  if (showCI){
    lines(R0tDifferenceQuantiles.underspec[1,showIdx], col = "red", lty = 2)
    lines(R0tDifferenceQuantiles.underspec[3,showIdx], col = "red", lty = 2)
    lines(EARDifferenceQuantiles.underspec[1,showIdx], lty = 2)
    lines(EARDifferenceQuantiles.underspec[3,showIdx], lty = 2)
  }
  
  
  plot(EARDifferenceQuantiles.spec[2,showIdx], type = "l", lwd = 3, ylim = c(-1,1),
       xlab = "Epidemic Time", ylab = "Bias", main = "Median Bias Under Correct Specification")
  abline(h = 0, col = "blue", lty = 2, lwd = 2)
  abline(h = seq(-5, 5, 0.5), lty = 2, col = rgb(0,0,0,0.2))
  lines(EARDifferenceQuantiles.spec[2,showIdx])
  lines(R0tDifferenceQuantiles.spec[2,showIdx], col = "red", lwd = 3, lty = 3)
  lines(EffR0tDifferenceQuantiles.spec[2,showIdx], col = "grey", lwd = 3, lty = 4)
  if (showCI){
    lines(R0tDifferenceQuantiles.spec[1,showIdx], col = "red", lty = 2)
    lines(R0tDifferenceQuantiles.spec[3,showIdx], col = "red", lty = 2)
    lines(EARDifferenceQuantiles.spec[1,showIdx], lty = 2)
    lines(EARDifferenceQuantiles.spec[3,showIdx], lty = 2)
  }
  par(xaxt = "n")
  par(yaxt = "n")
  par(bty = "n")
  par(xpd = NA)
  plot(c(-4, 4), c(-1,1), type = "n", xlab = "", ylab = "", main = "", xlim = c(-4,4), ylim = c(-1,1))
  legend(x = 0, y = 0, lty = c(1,3,4), lwd = c(3,3,3), col = c("black","red", "grey"), 
         legend = c("$\\mathcal{R}^{(EA)}(t)$", "$\\mathcal{R}_0(t)$", "effective $\\mathcal{R}_0(t)$"),
         xjust = 0.5, yjust =0.5, horiz=TRUE, cex = 1.25)
  par(xaxt = "s")
  par(yaxt = "s")
  par(bty = "o")
  par(xpd = FALSE)
dev.off()



tikz(file = "./ABMBiasPlot2.tex", width = 12, height = 8, standAlone=TRUE)
showCI = FALSE
showIdx = 2:ncol(EARDifferenceQuantiles.underspec)
layout(matrix(c(1,2,
                3,3), nrow = 2, ncol = 2, byrow = TRUE), 
       widths = c(8,8), heights = c(8,2)) 
plot(EARDifferenceQuantiles.underspec[4,showIdx], type = "l", lwd = 3, ylim = c(-1,1),
     xlab = "Epidemic Time", ylab = "Bias", main = "Mean Bias Under Misspecification")
abline(h = 0, col = "blue", lty = 2, lwd = 2)
abline(h = seq(-5, 5, 0.5), lty = 2, col = rgb(0,0,0,0.2))
lines(EARDifferenceQuantiles.underspec[4,showIdx])
lines(R0tDifferenceQuantiles.underspec[4,showIdx], col = "red", lwd = 3, lty = 3)
lines(EffR0tDifferenceQuantiles.underspec[2,showIdx], col = "grey", lwd = 3, lty = 4)
if (showCI){
  lines(R0tDifferenceQuantiles.underspec[1,showIdx], col = "red", lty = 2)
  lines(R0tDifferenceQuantiles.underspec[3,showIdx], col = "red", lty = 2)
  lines(EARDifferenceQuantiles.underspec[1,showIdx], lty = 2)
  lines(EARDifferenceQuantiles.underspec[3,showIdx], lty = 2)
}


plot(EARDifferenceQuantiles.spec[4,showIdx], type = "l", lwd = 3, ylim = c(-1,1),
     xlab = "Epidemic Time", ylab = "Bias", main = "Mean Bias Under Correct Specification")
abline(h = 0, col = "blue", lty = 2, lwd = 2)
abline(h = seq(-5, 5, 0.5), lty = 2, col = rgb(0,0,0,0.2))
lines(EARDifferenceQuantiles.spec[4,showIdx])
lines(R0tDifferenceQuantiles.spec[4,showIdx], col = "red", lwd = 3, lty = 3)
lines(EffR0tDifferenceQuantiles.spec[2,showIdx], col = "grey", lwd = 3, lty = 4)
if (showCI){
  lines(R0tDifferenceQuantiles.spec[1,showIdx], col = "red", lty = 2)
  lines(R0tDifferenceQuantiles.spec[3,showIdx], col = "red", lty = 2)
  lines(EARDifferenceQuantiles.spec[1,showIdx], lty = 2)
  lines(EARDifferenceQuantiles.spec[3,showIdx], lty = 2)
}
par(xaxt = "n")
par(yaxt = "n")
par(bty = "n")
par(xpd = NA)
plot(c(-4, 4), c(-1,1), type = "n", xlab = "", ylab = "", main = "", xlim = c(-4,4), ylim = c(-1,1))
legend(x = 0, y = 0, lty = c(1,3,4), lwd = c(3,3,3), col = c("black","red", "grey"), 
       legend = c("$\\mathcal{R}^{(EA)}(t)$", "$\\mathcal{R}_0(t)$", "effective $\\mathcal{R}_0(t)$"),
       xjust = 0.5, yjust =0.5, horiz=TRUE, cex = 1.25)
par(xaxt = "s")
par(yaxt = "s")
par(bty = "o")
par(xpd = FALSE)
dev.off()


#plot(abs(EARDifferenceQuantiles.spec[2,]), type = "l", lwd = 2, ylim = c(0,2))
#lines(abs(R0tDifferenceQuantiles.spec[2,]), col = "red", lwd = 2)
mean(abs(EARDifferenceQuantiles.underspec[2,]) < abs(R0tDifferenceQuantiles.underspec[2,]))
mean(abs(EARDifferenceQuantiles.underspec[2,]) < abs(EffR0tDifferenceQuantiles.underspec[2,]))
mean(abs(EARDifferenceQuantiles.spec[2,]) < abs(R0tDifferenceQuantiles.spec[2,]))
mean(abs(EARDifferenceQuantiles.spec[2,]) < abs(EffR0tDifferenceQuantiles.spec[2,]))




concordance.uspec = Reduce("+", 
                     lapply(R0Results[grep("TRUE",fnames)], function(x){
                       matrix(c(x$r0.concordance, x$EAR.concordance, x$EffR0.concordance), nrow = 1, byrow = TRUE) 
                     }))/length(grep("TRUE",fnames))
colnames(concordance.uspec) = c("R0", "EAR", "EffR0")
concordance.spec = Reduce("+", 
                           lapply(R0Results[grep("FALSE",fnames)], function(x){
                             matrix(c(x$r0.concordance, x$EAR.concordance, x$EffR0.concordance), nrow = 1, byrow = TRUE) 
                           }))/length(grep("FALSE",fnames))
colnames(concordance.spec) = c("R0", "EAR", "EffR0")
concordance = rbind(concordance.uspec, concordance.spec)
rownames(concordance) = c("Spec", "Underspec")


bias.uspec = Reduce("+", 
                           lapply(R0Results[grep("TRUE",fnames)], function(x){
                             matrix(c(x$EARvsR0, x$EARvsEffR0), nrow = 1, byrow = TRUE) 
                           }))/length(grep("TRUE",fnames))
colnames(bias.uspec) = c("EARvsR0", "EARvsEffR0")
bias.spec = Reduce("+", 
                          lapply(R0Results[grep("FALSE",fnames)], function(x){
                            matrix(c(x$EARvsR0, x$EARvsEffR0), nrow = 1, byrow = TRUE) 
                          }))/length(grep("FALSE",fnames))
colnames(bias.spec) = c("EARvsR0", "EARvsEffR0")
bias = rbind(bias.spec, bias.uspec)
rownames(bias) = c("Spec", "Underspec")


write.csv(bias, file = "./bias.csv")
write.csv(concordance, file = "./concordance.csv")

