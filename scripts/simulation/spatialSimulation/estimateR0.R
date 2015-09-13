library(dplyr)
library(EARNMC)
library(parallel)


trueR0t = function(simResult)
{
  eta = matrix(exp(simResult$X %*% simResult$beta_SE), ncol = ncol(simResult$S))
  gamma_ir = simResult$gamma_ir
  eta/gamma_ir
}

trueEAR = function(simResult)
{
  S = simResult$S
  I = simResult$I
  N = simResult$N
  nLoc = ncol(I)
  nTpt = nrow(I)
  
  rho = simResult$rho
  eta.SE = matrix(exp(simResult$X %*% simResult$beta_SE), ncol = ncol(simResult$S))
  gamma_ir = simResult$gamma_ir
  p_ir = 1-exp(-gamma_ir)
    
  # Calculate G
  eta.R0  = I/N*eta.SE
  
  G = array(0, dim=c(nLoc, nLoc, nTpt))
  dmList = simResult$distMatList
  for (i in 1:nTpt)
  {
    for (j in 1:nLoc)
    {
      for (l in 1:nLoc)
      {
        if (j != l)
        {
          for (z in 1:length(dmList))
          {
            G[j,l,i] = G[j,l,i] + rho[z]*dmList[[z]][j,l]*eta.R0[i,l]
          }
        }
        else
        {
          G[j,j,i] = eta.R0[i,l]
        }
      }
    }
  }
  G = 1-exp(-G)
  S_tmp = array(t(S)[rep(1:3, each = 3),], dim = c(3,3,150))
  I_tmp = array(t(I)[rep(1:3, each = 3),], dim = c(3,3,150))
  EAR.components = S_tmp*G
  EAR.components = ifelse(I_tmp == 0, 0, EAR.components/I_tmp)
  EAR.components = t(apply(EAR.components,c(1,3), sum))
  
  
  EA_R0 = 0*EAR.components
  for (t in 1:nTpt)
  {
    cpir = 1
    idx = 0
    
    while(cpir > 1e-16 && idx <= nTpt - t)
    {
      #EA_R0[t,] = (EA_R0[t,] + EAR.components[t+idx,]*cpir)*(EAR.components[t,] != 0)
      EA_R0[t,] = (EA_R0[t,] + EAR.components[t+idx,]*cpir)
      cpir = cpir * (1-p_ir)
      idx = idx + 1
    }
    if (idx >= nTpt){
      while(cpir > 1e-16){
        #EA_R0[t,] = EA_R0[t,] + EAR.components[nTpt,]*(EAR.components[t,] != 0)*cpir
        EA_R0[t,] = EA_R0[t,] + EAR.components[nTpt,]*cpir
        cpir = cpir * (1-p_ir)
      }
    }
  }
  EA_R0
}

estimateR0 = function(dataSeed)
{
  fname1 = paste("./results/simSp_results_Uspec_FALSE_", dataSeed, ".Rda.bz2", sep = "")
  fname2 = paste("./results/simSp_results_Uspec_TRUE_", dataSeed, ".Rda.bz2", sep = "")
  
  if (!(file.exists(fname1) && file.exists(fname2)))
  {
    cat(paste("Missing one or more files for data seed: ", dataSeed, ", skipping\n", sep = ""))
    return(FALSE)
  }
  
  r1 = estimateR0_(fname1)
  r2 = estimateR0_(fname2)
  
  calculateBias = function(rslt)
  {
    trueR0Array = array(rslt$trueR0, dim = c(dim(rslt$trueR0), dim(rslt$EARMatrix)[3]))
    trueEARArray = array(rslt$trueEAR, dim = c(dim(rslt$trueEAR), dim(rslt$EARMatrix)[3]))
    
    r0.bias = rslt$R0tMatrix - trueR0Array
    EAR.bias = rslt$EARMatrix - trueR0Array
    EffR0.bias = rslt$EffR0tMatrix - trueR0Array
    
    r0.bias.EAR = rslt$R0tMatrix - trueEARArray
    EAR.bias.EAR = rslt$EARMatrix - trueEARArray
    EffR0.bias.EAR = rslt$EffR0tMatrix - trueEARArray
    
    r0.concordance = (rslt$R0tMatrix >= 1 &  trueR0Array >= 1) | (rslt$R0tMatrix < 1 &  trueR0Array < 1)
    EAR.concordance = (rslt$EARMatrix >= 1 &  trueR0Array >= 1) | (rslt$EARMatrix < 1 &  trueR0Array < 1)
    EffR0.concordance = (rslt$EffR0tMatrix >= 1 &  trueR0Array >= 1) | (rslt$EffR0tMatrix < 1 &  trueR0Array < 1)
  
    r0.concordance.EAR = (rslt$R0tMatrix >= 1 &  trueEARArray >= 1) | (rslt$R0tMatrix < 1 &  trueEARArray < 1)
    EAR.concordance.EAR = (rslt$EARMatrix >= 1 &  trueEARArray >= 1) | (rslt$EARMatrix < 1 &  trueEARArray < 1)
    EffR0.concordance.EAR = (rslt$EffR0tMatrix >= 1 &  trueEARArray >= 1) | (rslt$EffR0tMatrix < 1 &  trueEARArray < 1)
    
    
    
    list(
      r0.bias.quantiles = apply(r0.bias, 1, quantile, probs = c(0.05, 0.5, 0.95)),
      EAR.bias.quantiles = apply(EAR.bias, 1, quantile, probs = c(0.05, 0.5, 0.95)),
      EffR0.bias.quantiles = apply(EffR0.bias, 1, quantile, probs = c(0.05, 0.5, 0.95)),
      EARvsR0 = mean(abs(EAR.bias) < abs(r0.bias)),
      EARvsEffR0 = mean(abs(EAR.bias) < abs(EffR0.bias)),
      r0.concordance=mean(r0.concordance),
      EAR.concordance=mean(EAR.concordance),
      EffR0.concordance=mean(EffR0.concordance),
      R0tMatrix = rslt$R0tMatrix,
      EARMatrix = rslt$EARMatrix,
      EffR0Matrix = rslt$EffR0tMatrix,
      r0.bias.EAR.quantiles = apply(r0.bias.EAR, 1, quantile, probs = c(0.05, 0.5, 0.95)),
      EAR.bias.EAR.quantiles = apply(EAR.bias.EAR, 1, quantile, probs = c(0.05, 0.5, 0.95)),
      EffR0.bias.EAR.quantiles = apply(EffR0.bias.EAR, 1, quantile, probs = c(0.05, 0.5, 0.95)),
      r0.concordance.EAR = mean(r0.concordance.EAR),
      EAR.concordance.EAR = mean(EAR.concordance.EAR),
      EffR0.concordance.EAR = mean(EffR0.concordance.EAR),
      EARvsR0.EAR = mean(abs(EAR.bias.EAR) < abs(r0.bias.EAR)),
      EARvsEffR0.EAR = mean(abs(EAR.bias.EAR) < abs(EffR0.bias.EAR)),
      trueR0Array=trueR0Array,
      trueEARArray=trueEARArray
    )
  }
  list(bias1 = calculateBias(r1),
       bias2 = calculateBias(r2))
}

estimateR0_ = function(dataName)
{

  load(dataName)
  nSamples = nrow(result$dat1)
  N = result$simResults$N[1,]
  out = lapply((floor(nSamples/2)):nSamples, 
               function(x){
                 r0.i(x,result$dat1, result$simResults)
               })

  underspec = length(grep("rho", colnames(result$dat1))) == 1
  nLoc = length(N)
  #nLoc = ifelse(underspec, 1,  result$simResults$nLoc)
  trueR0 = matrix(trueR0t(result$simResults), nrow = nrow(result$simResults$S_star))
  trueEAR = trueEAR(result$simResults)
  EARMatrix = array(unlist(
                    lapply(out, function(x){matrix(x$EAR, nrow = nrow(result$simResults$S))})), 
                    dim = c(result$simResults$nTpt, nLoc, length(out)))
  R0tMatrix = array(unlist(
    lapply(out, function(x){matrix(x$R0t, nrow = nrow(result$simResults$S))})), 
    dim = c(result$simResults$nTpt, nLoc, length(out)))
  
  EffR0tMatrix = array(unlist(
    lapply(out, function(x){matrix(x$effR0t, nrow = nrow(result$simResults$S))})), 
    dim = c(result$simResults$nTpt, nLoc, length(out)))
  
  list(trueR0=trueR0, 
       EARMatrix=EARMatrix, 
       R0tMatrix=R0tMatrix, 
       EffR0tMatrix=EffR0tMatrix,
       trueEAR=trueEAR,
       N=N)
}

r0.i = function(i, chainObject, simResults){
  startPt=0
  dataRow = chainObject[i,]
  
  rho = as.numeric(dataRow[grep("rho_", 
                                colnames(dataRow))])
  underspec = (length(rho) == 1)
  nLoc = ncol(simResults$S)
  nTpt = nrow(simResults$S)
  
  S = matrix(unlist(lapply(0:(nLoc - 1), 
                           function(x){
                             as.numeric(
                               select(dataRow, 
                                      starts_with(paste("S_", x, "_", sep = ""))))
                             }
                           )
                    ),
             nrow=nTpt, ncol = nLoc)
  E = matrix(unlist(lapply(0:(nLoc - 1), 
                           function(x){
                             as.numeric(
                               select(dataRow, 
                                      starts_with(paste("E_", x, "_", sep = ""))))
                           }
                           )
                    ),
             nrow=nTpt, ncol = nLoc)
  I = matrix(unlist(lapply(0:(nLoc - 1), 
                           function(x){
                             as.numeric(
                               select(dataRow, 
                                      starts_with(paste("I_", x, "_", sep = ""))))
                           }
                        )
                    ),
             nrow=nTpt, ncol = nLoc)
  I_star = matrix(unlist(lapply(0:(nLoc - 1), 
                           function(x){
                             as.numeric(
                               select(dataRow, 
                                      starts_with(paste("I_star_", x, "_", sep = ""))))
                           }
  )
  ),
  nrow=nTpt, ncol = nLoc)
  R = matrix(unlist(lapply(0:(nLoc - 1), 
                           function(x){
                             as.numeric(
                               select(dataRow, 
                                      starts_with(paste("R_", x, "_", sep = ""))))
                           }
                           )
                    ),
             nrow=nTpt, ncol = nLoc)
  
  N = S + E + I + R

  beta.SE = as.numeric(dataRow[grep("BetaP_SE_", 
                                    colnames(dataRow))])
  X_se = simResults$X
  if (length(beta.SE) != ncol(X_se))
  {
    stop("Incorrect number of columns")
    X_se = simResults$X_underspec
  }
  eta.SE = matrix(X_se %*% beta.SE, ncol=nLoc)
  gamma_ei = dataRow$gamma_ei
  gamma_ir = dataRow$gamma_ir
  
  
  # Calculate G
  
  eta.R0  = I/N*exp(eta.SE)
  
  G = array(0, dim=c(nLoc, nLoc, nTpt))
  dmList = ifelse(underspec, simResults$distMatList_uspec, simResults$distMatList)
  for (i in 1:nTpt)
  {
    for (j in 1:nLoc)
    {
      for (l in 1:nLoc)
      {
        if (j != l)
        {
          for (z in 1:length(dmList))
          {
            G[j,l,i] = G[j,l,i] + rho[z]*dmList[[z]][j,l]*eta.R0[i,l]
          }
        }
        else
        {
          G[j,j,i] = eta.R0[i,l]
        }
      }
    }
  }
  G = 1-exp(-G)
  S_tmp = array(t(S)[rep(1:3, each = 3),], dim = c(3,3,150))
  I_tmp = array(t(I)[rep(1:3, each = 3),], dim = c(3,3,150))
  EAR.components = S_tmp*G
  EAR.components = ifelse(I_tmp == 0, 0, EAR.components/I_tmp)
  EAR.components = t(apply(EAR.components,c(1,3), sum))
  
  #EAR.components = S*t(apply(G, c(1,3), sum))
  #EAR.components = ifelse(I == 0, 0, EAR.components/I)
  
  p_ir = 1-exp(-gamma_ir)
  EA_R0 = 0*EAR.components
  ##
  #eta = matrix(exp(eta.SE), ncol = ncol(S))
  #dm = diag(3) + rho[1]*simResults$distMatList[[1]] + rho[2]*simResults$distMatList[[2]]
  #gamma_ir = simResult$gamma_ir
  #t(apply(eta, 1, function(x){
  #  xm = matrix(x, nrow = 3, ncol =3, byrow = TRUE)*dm/gamma_ir
  #  apply(xm, 1, sum)
  #}))
  ##
  R0_t = matrix(exp(eta.SE)/gamma_ir, ncol = ncol(S))
  Eff_R0 = R0_t*(S/N)
  
  
  for (t in 1:nTpt)
  {
    cpir = 1
    idx = 0
    while(cpir > 1e-16 && idx <= nTpt - t)
    {
      EA_R0[t,] = EA_R0[t,] + EAR.components[t+idx,]*cpir
      cpir = cpir * (1-p_ir)
      idx = idx + 1
    }
    if (idx >= nTpt){
      while(cpir > 1e-16){
        EA_R0[t,] = EA_R0[t,] + EAR.components[nTpt,]*cpir
        cpir = cpir * (1-p_ir)
      }
    }
  }
  list(EAR = EA_R0,
       R0t = R0_t, 
       effR0t = Eff_R0,
       S=S,
       E=E,
       I=I,
       R=R,
       I_star = I_star)
}

tmp.fnames = dir("./results/")
tmp.fnames = tmp.fnames[tmp.fnames != "README.txt"]

dataSeeds = unique(gsub(".Rda.bz2", "", 
                        gsub("simSp_results_Uspec_TRUE_", "", 
                             gsub("simSp_results_Uspec_FALSE_", "",
                                  tmp.fnames))))

#cl = makeCluster(6, outfile = "./err.txt")
#parLapply(cl, 1:length(cl), function(x){library(dplyr); library(EARNMC)})
#clusterExport(cl, list("r0.i", "estimateR0_", "trueR0t", "trueR0t"))

f=function(x){
  cat(paste(x, "\n", sep = ""));
  fname = paste(x, "_processed.Robj", sep = "")
  if (file.exists(fname))
  {
    cat("...Already processed\n")
    return(TRUE)
  }
  else
  {
    result = estimateR0(x)
    save("result", file = fname)
    return(TRUE)
  }
}
read.res = function(x)
{
  load(paste(x, "_processed.Robj", sep = ""))
  result
}
R0Results = lapply(dataSeeds, f)
R0Results = lapply(dataSeeds, read.res)
#R0Results = parLapplyLB(cl, dataSeeds, estimateR0)
R0Results.cleaned = R0Results[which(as.logical(lapply(R0Results, function(x){length(x) != 1})))]
save(file="./results.Robj", "R0Results.cleaned")

concordance = Reduce("+", 
       lapply(R0Results.cleaned, function(x){
  b1 = x$bias1
  b2 = x$bias2
  matrix(c(b1$r0.concordance, b1$EAR.concordance, b1$EffR0.concordance, 
  b2$r0.concordance, b2$EAR.concordance, b2$EffR0.concordance), nrow = 2, byrow = TRUE) 
}))/length(R0Results.cleaned)
colnames(concordance) = c("R0", "EAR", "EffR0")
rownames(concordance) = c("Spec", "Underspec")


concordance.EAR = Reduce("+", 
                     lapply(R0Results.cleaned, function(x){
                       b1 = x$bias1
                       b2 = x$bias2
                       matrix(c(b1$r0.concordance.EAR, b1$EAR.concordance.EAR, b1$EffR0.concordance.EAR, 
                                b2$r0.concordance.EAR, b2$EAR.concordance.EAR, b2$EffR0.concordance.EAR), nrow = 2, byrow = TRUE) 
                     }))/length(R0Results.cleaned)
colnames(concordance.EAR) = c("R0", "EAR", "EffR0")
rownames(concordance.EAR) = c("Spec", "Underspec")



bias_comparison = Reduce("+", 
                     lapply(R0Results.cleaned, function(x){
                       b1 = x$bias1
                       b2 = x$bias2
                       matrix(c(b1$EARvsR0, b1$EARvsEffR0, 
                                b2$EARvsR0, b2$EARvsEffR0), nrow = 2, byrow = TRUE) 
                     }))/length(R0Results.cleaned)
colnames(bias_comparison) = c("R0", "EffR0")
rownames(bias_comparison) = c("Spec", "Underspec")


bias_comparison.EAR = Reduce("+", 
                         lapply(R0Results.cleaned, function(x){
                           b1 = x$bias1
                           b2 = x$bias2
                           matrix(c(b1$EARvsR0.EAR, b1$EARvsEffR0.EAR, 
                                    b2$EARvsR0.EAR, b2$EARvsEffR0.EAR), nrow = 2, byrow = TRUE) 
                         }))/length(R0Results.cleaned)
colnames(bias_comparison.EAR) = c("R0", "EffR0")
rownames(bias_comparison.EAR) = c("Spec", "Underspec")

write.csv(concordance, "./concordance.csv")
write.csv(bias_comparison, "./bias.csv")
write.csv(concordance.EAR, "./concordance_EAR.csv")
write.csv(bias_comparison.EAR, "./bias_EAR.csv")

EARArray.1 = matrix(Reduce("+", lapply(R0Results.cleaned, 
			function(x){
				apply(x$bias1$EARMatrix, 1:2, mean)
			}))/length(R0Results.cleaned), nrow = 150, ncol = 3)
EARArray.2 = matrix(Reduce("+", lapply(R0Results.cleaned, 
			function(x){
				apply(x$bias2$EARMatrix, 1:2, mean)
			}))/length(R0Results.cleaned), nrow = 150, ncol = 3)

R0Array.1 = matrix(Reduce("+", lapply(R0Results.cleaned, 
			function(x){
				apply(x$bias1$R0tMatrix, 1:2, mean)
			}))/length(R0Results.cleaned), nrow = 150, ncol = 3)
R0Array.2 = matrix(Reduce("+", lapply(R0Results.cleaned, 
			function(x){
				apply(x$bias2$R0tMatrix, 1:2, mean)
			}))/length(R0Results.cleaned), nrow = 150, ncol = 3)

EffR0Array.1 = matrix(Reduce("+", lapply(R0Results.cleaned, 
			function(x){
				apply(x$bias1$EffR0Matrix, 1:2, mean)
			}))/length(R0Results.cleaned), nrow = 150, ncol = 3)
EffR0Array.2 = matrix(Reduce("+", lapply(R0Results.cleaned, 
			function(x){
				apply(x$bias2$EffR0Matrix, 1:2, mean)
			}))/length(R0Results.cleaned), nrow = 150, ncol = 3)

trueR0Matrix = R0Results.cleaned[[1]]$bias1$trueR0Array[,,1]
trueEARMatrix = R0Results.cleaned[[1]]$bias1$trueEARArray[,,1]

plot(R0Array.1[,1], type = "l", lty = 1, ylim = c(0, 4))
lines(R0Array.1[,2], lty = 2)
lines(R0Array.1[,3], lty = 3)

lines(EARArray.1[,1], col = "red", lty = 1)
lines(EARArray.1[,2], col = "red", lty = 2)
lines(EARArray.1[,3], col = "red", lty = 3)

#lines(EffR0Array.1[,1], col = "blue", lty = 1)
#lines(EffR0Array.1[,2], col = "blue", lty = 2)
#lines(EffR0Array.1[,3], col = "blue", lty = 3)

lines(trueR0Matrix[,1], col = "green", lty = 1)
lines(trueR0Matrix[,2], col = "green", lty = 2)
lines(trueR0Matrix[,3], col = "green", lty = 3)


lines(trueEARMatrix[,1], col = "blue", lty = 1)
lines(trueEARMatrix[,2], col = "blue", lty = 2)
lines(trueEARMatrix[,3], col = "blue", lty = 3)



plot(R0Array.2[,1], type = "l", lty = 1, ylim = c(0,4))
lines(R0Array.2[,2], lty = 2)
lines(R0Array.2[,3], lty = 3)

lines(EARArray.2[,1], col = "red", lty = 1)
lines(EARArray.2[,2], col = "red", lty = 2)
lines(EARArray.2[,3], col = "red", lty = 3)


lines(trueR0Matrix[,1], col = "green", lty = 1)
lines(trueR0Matrix[,2], col = "green", lty = 2)
lines(trueR0Matrix[,3], col = "green", lty = 3)


lines(trueEARMatrix[,1], col = "blue", lty = 1)
lines(trueEARMatrix[,2], col = "blue", lty = 2)
lines(trueEARMatrix[,3], col = "blue", lty = 3)






