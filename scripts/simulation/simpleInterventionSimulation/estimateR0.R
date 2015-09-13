library(dplyr)
library(EARNMC)
library(parallel)
trueR0t = function(simResult)
{
  eta = simResult$X %*% simResult$beta_SE
  gamma_ir = simResult$gamma_ir
  return(exp(eta)/gamma_ir)
}


trueEARt = function(simResult)
{
  S = simResult$S
  I = simResult$I
  N = simResult$N
  
  nTpt = length(I)
  
  eta.SE = matrix(exp(simResult$X %*% simResult$beta_SE), ncol = ncol(simResult$S))
  gamma_ir = simResult$gamma_ir
  p_ir = 1-exp(-gamma_ir)
  
  # Calculate G
  eta.R0  = I/N*eta.SE
  
  G = 1-exp(-eta.R0)
  
  EAR.components = S*G
  EAR.components = ifelse(I == 0, 0, EAR.components/I)
    
  EA_R0 = 0*EAR.components
  for (t in 1:nTpt)
  {
    cpir = 1
    idx = 0
    
    while(cpir > 1e-16 && idx <= nTpt - t)
    {
      #EA_R0[t,] = (EA_R0[t,] + EAR.components[t+idx,]*cpir)*(EAR.components[t,] != 0)
      EA_R0[t,] = (EA_R0[t] + EAR.components[t+idx]*cpir)
      cpir = cpir * (1-p_ir)
      idx = idx + 1
    }
    if (idx >= nTpt){
      while(cpir > 1e-16){
        #EA_R0[t,] = EA_R0[t,] + EAR.components[nTpt,]*(EAR.components[t,] != 0)*cpir
        EA_R0[t,] = EA_R0[t] + EAR.components[nTpt]*cpir
        cpir = cpir * (1-p_ir)
      }
    }
  }
  EA_R0
}

estimateR0 = function(dataName)
{
  load(dataName)
  nSamples = nrow(result$dat1)
  out = lapply((floor(nSamples/2)):nSamples, 
               function(x){
                 r0.i(x,result$dat1, result$simResults$X)
               })

  
  trueR0 = matrix(trueR0t(result$simResults), nrow = nrow(result$simResults$S_star),
                  ncol=(nSamples - floor(nSamples/2)) + 1)
  
  trueEAR = matrix(trueEARt(result$simResults), nrow = nrow(result$simResults$S_star),
                   ncol=(nSamples - floor(nSamples/2)) + 1)
  
  EARMatrx = matrix(unlist(
                    lapply(out, function(x){x$EAR})), 
                    nrow = result$simResults$nTpt)
  
  R0tMatrix = matrix(unlist(
    lapply(out, function(x){x$R0t})), 
    nrow = result$simResults$nTpt)
  
  EffR0tMatrix = matrix(unlist(
    lapply(out, function(x){x$effR0t})), 
    nrow = result$simResults$nTpt)
  
  r = cor(apply(EARMatrx, 1, mean), apply(EffR0tMatrix, 1, mean))  
  
  sm.fnc = function(x)
  {
    rbind(apply(x, 1, quantile, probs = c(0.05, 0.5, 0.95)),
          apply(x, 1, mean))
  }
  EARQuantiles = sm.fnc(EARMatrx)
  R0tQuantiles = sm.fnc(R0tMatrix)
  EffR0tQuantiles = sm.fnc(EffR0tMatrix)
  EffR0tDifferenceQuantiles = sm.fnc(EffR0tMatrix-trueR0)
  EARDifferenceQuantiles =  sm.fnc(EARMatrx-trueR0)
  R0tDifferenceQuantiles =  sm.fnc(R0tMatrix-trueR0)
  EARConcordanceQuantiles = apply((EARMatrx >= 1 & trueR0 >= 1) | (EARMatrx < 1 & trueR0 < 1), 1, mean)
  R0tConcordanceQuantiles = apply((R0tMatrix >= 1 & trueR0 >= 1) | (R0tMatrix < 1 & trueR0 < 1), 1, mean)

  EffR0tDifferenceQuantiles.EAR = sm.fnc(EffR0tMatrix-trueEAR)
  EARDifferenceQuantiles.EAR =  sm.fnc(EARMatrx-trueEAR)
  R0tDifferenceQuantiles.EAR =  sm.fnc(R0tMatrix-trueEAR)
  EARConcordanceQuantiles.EAR = apply((EARMatrx >= 1 & trueEAR >= 1) | (EARMatrx < 1 & trueEAR < 1), 1, mean)
  R0tConcordanceQuantiles.EAR = apply((R0tMatrix >= 1 & trueEAR >= 1) | (R0tMatrix < 1 & trueEAR < 1), 1, mean)
  
  
  
  list(EARQuantiles=EARQuantiles,
       R0tQuantiles=R0tQuantiles,
       EffR0tQuantiles=EffR0tQuantiles,
       R0tDifferenceQuantiles=R0tDifferenceQuantiles,
       EARDifferenceQuantiles=EARDifferenceQuantiles,
       EffR0tDifferenceQuantiles=EffR0tDifferenceQuantiles,
       EARConcordanceQuantiles=EARConcordanceQuantiles,
       R0tConcordanceQuantiles=R0tConcordanceQuantiles,
       R0tDifferenceQuantiles.EAR=R0tDifferenceQuantiles.EAR,
       EARDifferenceQuantiles.EAR=EARDifferenceQuantiles.EAR,
       EffR0tDifferenceQuantiles.EAR=EffR0tDifferenceQuantiles.EAR,
       EARConcordanceQuantiles.EAR=EARConcordanceQuantiles.EAR,
       R0tConcordanceQuantiles.EAR=R0tConcordanceQuantiles.EAR,
       r=r
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
        while(cpir > 1e-16){
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

cl = makeCluster(6, type = "PSOCK")
parLapply(cl, 1:length(cl), function(x){library(dplyr); library(EARNMC)})
clusterExport(cl, c("estimateR0", "r0.i", "trueR0t", "trueEARt"))
tmp.fnames = dir("./results/")
tmp.fnames = tmp.fnames[tmp.fnames != "README.txt"]
R0Results = parLapplyLB(cl, tmp.fnames, function(x){cat(paste(x, "\n", sep = ""));
                                                  estimateR0(paste("./results/", x, sep = ""))})

underspec.idx = grepl("Uspec_TRUE", tmp.fnames)
spec.idx = grepl("Uspec_FALSE", tmp.fnames)

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

EARConcordanceQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                                       function(x){R0Results[[x]][[7]]}))/sum(underspec.idx)
R0ConcordanceQuantiles.underspec = Reduce("+", lapply(which(underspec.idx), 
                                                       function(x){R0Results[[x]][[8]]}))/sum(underspec.idx)

EARConcordanceQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                                       function(x){R0Results[[x]][[7]]}))/sum(spec.idx)
R0ConcordanceQuantiles.spec = Reduce("+", lapply(which(spec.idx), 
                                                      function(x){R0Results[[x]][[8]]}))/sum(spec.idx)




R0tDifferenceQuantiles.underspec.EAR = Reduce("+", lapply(which(underspec.idx), 
                                                      function(x){R0Results[[x]]$R0tDifferenceQuantiles.EAR}))/sum(underspec.idx)
R0tDifferenceQuantiles.spec.EAR = Reduce("+", lapply(which(spec.idx), 
                                                 function(x){R0Results[[x]]$R0tDifferenceQuantiles.EAR}))/sum(spec.idx)

EARDifferenceQuantiles.underspec.EAR = Reduce("+", lapply(which(underspec.idx), 
                                                      function(x){R0Results[[x]]$EARDifferenceQuantiles.EAR}))/sum(underspec.idx)
EARDifferenceQuantiles.spec.EAR = Reduce("+", lapply(which(spec.idx), 
                                                 function(x){R0Results[[x]]$EARDifferenceQuantiles.EAR}))/sum(spec.idx)

EffR0tDifferenceQuantiles.underspec.EAR = Reduce("+", lapply(which(underspec.idx), 
                                                         function(x){R0Results[[x]]$EffR0tDifferenceQuantiles.EAR}))/sum(underspec.idx)
EffR0tDifferenceQuantiles.spec.EAR = Reduce("+", lapply(which(spec.idx), 
                                                    function(x){R0Results[[x]]$EffR0tDifferenceQuantiles.EAR}))/sum(spec.idx)

EARConcordanceQuantiles.underspec.EAR = Reduce("+", lapply(which(underspec.idx), 
                                                       function(x){R0Results[[x]]$EARConcordanceQuantiles.EAR}))/sum(underspec.idx)
R0ConcordanceQuantiles.underspec.EAR = Reduce("+", lapply(which(underspec.idx), 
                                                      function(x){R0Results[[x]]$R0tConcordanceQuantiles.EAR}))/sum(underspec.idx)

EARConcordanceQuantiles.spec.EAR = Reduce("+", lapply(which(spec.idx), 
                                                  function(x){R0Results[[x]]$EARConcordanceQuantiles.EAR}))/sum(spec.idx)
R0ConcordanceQuantiles.spec.EAR = Reduce("+", lapply(which(spec.idx), 
                                                 function(x){R0Results[[x]]$R0tConcordanceQuantiles.EAR}))/sum(spec.idx)



mean(abs(EARDifferenceQuantiles.underspec[2,]) < abs(R0tDifferenceQuantiles.underspec[2,]))
mean(abs(EARDifferenceQuantiles.spec[2,]) < abs(R0tDifferenceQuantiles.spec[2,]))

mean(abs(EARDifferenceQuantiles.underspec.EAR[2,]) < abs(R0tDifferenceQuantiles.underspec.EAR[2,]))
mean(abs(EARDifferenceQuantiles.spec.EAR[2,]) < abs(R0tDifferenceQuantiles.spec.EAR[2,]))


mean(EARConcordanceQuantiles.spec)
mean(R0ConcordanceQuantiles.spec)
mean(EARConcordanceQuantiles.underspec)
mean(R0ConcordanceQuantiles.underspec)

mean(EARConcordanceQuantiles.spec.EAR)
mean(R0ConcordanceQuantiles.spec.EAR)
mean(EARConcordanceQuantiles.underspec.EAR)
mean(R0ConcordanceQuantiles.underspec.EAR)



r.uspec = Reduce("+", lapply(which(underspec.idx), 
                 function(x){R0Results[[x]]$r}))/sum(underspec.idx)
r.spec = Reduce("+", lapply(which(spec.idx), 
                             function(x){R0Results[[x]]$r}))/sum(spec.idx)

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
                                "R0Results"="R0Results")


require(tikzDevice)

tikz(file = "./SimpleInterventionBiasPlot.tex", width = 12, height = 8, standAlone=TRUE)
  layout(matrix(c(1,2,
                  3,3), nrow = 2, ncol = 2, byrow = TRUE), 
         widths = c(8,8), heights = c(8,2)) 
  plot(EARDifferenceQuantiles.underspec[2,], type = "l", lwd = 3, ylim = c(-2,2),
       xlab = "Epidemic Time", ylab = "Bias", main = "Median Bias Under Misspecification")
  abline(h = 0, col = "blue", lty = 2, lwd = 2)
  abline(h = seq(-5, 5, 0.5), lty = 2, col = rgb(0,0,0,0.2))
  lines(EARDifferenceQuantiles.underspec[2,])
  lines(R0tDifferenceQuantiles.underspec[2,], col = "red", lwd = 3, lty = 2)
  #lines(EffR0tDifferenceQuantiles.underspec[2,], col = "grey", lwd = 3, lty = 3)
  #lines(R0tDifferenceQuantiles.underspec[1,], col = "red", lty = 2)
  #lines(R0tDifferenceQuantiles.underspec[3,], col = "red", lty = 2)
  #lines(EARDifferenceQuantiles.underspec[1,], lty = 2)
  #lines(EARDifferenceQuantiles.underspec[3,], lty = 2)
  
  
  plot(EARDifferenceQuantiles.spec[2,], type = "l", lwd = 3, ylim = c(-2,2),
       xlab = "Epidemic Time", ylab = "Bias", main = "Median Bias Under Correct Specification")
  abline(h = 0, col = "blue", lty = 2, lwd = 2)
  abline(h = seq(-5, 5, 0.5), lty = 2, col = rgb(0,0,0,0.2))
  lines(EARDifferenceQuantiles.spec[2,])
  lines(R0tDifferenceQuantiles.spec[2,], col = "red", lwd = 3, lty = 2)
  #lines(EffR0tDifferenceQuantiles.spec[2,], col = "grey", lwd = 3, lty = 3)
  #lines(R0tDifferenceQuantiles.spec[1,], col = "red", lty = 2)
  #lines(R0tDifferenceQuantiles.spec[3,], col = "red", lty = 2)
  #lines(EARDifferenceQuantiles.spec[1,], lty = 2)
  #lines(EARDifferenceQuantiles.spec[3,], lty = 2)
  
  
  par(xaxt = "n")
  par(yaxt = "n")
  par(bty = "n")
  par(xpd = NA)
  plot(c(-4, 4), c(-1,1), type = "n", xlab = "", ylab = "", main = "", xlim = c(-4,4), ylim = c(-1,1))
  #legend(x = 0, y = 0, lty = c(1,2,3), lwd = c(2,2,2), col = c("black", "red", "grey"), 
  #      legend = c("Empirically Adjusted", expression('R'[0](t))),
  #      xjust = 0.5, yjust =0.5, horiz=TRUE, cex = 1.25)

  legend(x = 0, y = 0, lty = c(1,2), lwd = c(3,3), col = c("black","red"), 
       legend = c("$\\mathcal{R}^{(EA)}(t)$", "$\\mathcal{R}_0(t)$"),
       xjust = 0.5, yjust =0.5, horiz=TRUE, cex = 1.25)
  par(xaxt = "s")
  par(yaxt = "s")
  par(bty = "o")
  par(xpd = FALSE)
dev.off()



