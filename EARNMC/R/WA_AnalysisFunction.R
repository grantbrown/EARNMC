WestAfricaAnalysisScript = function(modelDF, fileName = NA, modelMode=c("DF", "knots"), pred.days=120, 
                                    throwAwayTpts = 0){  
  ## set number of samples/batches
  numBurnInBatches =  1000
  convergenceCriterion =  1.02
  convergenceSampleSize = 100000
  convergenceBatchSize = 100000
  minimumSamples =  500000
  extraR0Iterations = 500
  extraR0BatchSize = 1000
  iterationStride = 1000
  targetDaysPerRecord = 7
  totalSamples = 0
  
  ## Read in the data
  if (all(is.na(fileName))){
    print("Using included data")
    data(WestAfricaEbola)
  }
  else{
    westAfricaEbola = read.csv(fileName)  
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
  uncumulate = function(x)
  {
    out = c(x[2:length(x)]-x[1:(length(x)-1)])
    ifelse(out >= 0, out, 0)
  }
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
  
  Guinea = Guinea[keepIdx]
  Sierra.Leone = Sierra.Leone[keepIdx]
  Liberia = Liberia[keepIdx]
  
  original.rptDate = original.rptDate[keepIdx]
  rptDate = original.rptDate[2:length(original.rptDate)]
  
  I_star = cbind(uncumulate(Guinea), 
                 uncumulate(Liberia), 
                 uncumulate(Sierra.Leone))
  I0 = c(Guinea[1], Liberia[1], Sierra.Leone[1])
  
  # If requested, get rid of some of the most recent measurements
  # in order to assess prediction performance. 
  I_star = I_star[1:(nrow(I_star) - throwAwayTpts),]
  rptDate = rptDate[1:(nrow(I_star))]
  
  # Define the temporal offset vector to be the number of days reflected in each 
  # aggregated record (time between reports).
  offsets = uncumulate(original.rptDate)[1:nrow(I_star)]
  if (any(offsets <= 0))
  {
    cat("Invalid Date Information. The data source has likely changed.\n")
    stop(-1)
  }
  
  maxIdx = nrow(I_star)
  # Guinea, Liberia, Sierra Leone, Nigeria
  N = matrix(c(10057975, 4128572, 6190280), nrow = nrow(I_star),ncol = 3, 
             byrow=TRUE)
  X = diag(ncol(N))
  X.predict = X
  
  daysSinceJan = as.numeric(rptDate - as.Date("2014-01-01"))
  daysSinceJan.predict = c(max(daysSinceJan) + 1, max(daysSinceJan) 
                           + seq(2,pred.days-2,2))
{
    if (modelDF != 0 && modelMode == "DF")
    {
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
    }
    else if (modelDF != 0 && modelMode == "Knots")
    {
      splineBasis = ns(daysSinceJan, knots = floor(quantile(daysSinceJan, 
                                                            seq(0,1, length = modelDF)))[1:(modelDF-1)])
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
    }
    else
    {
      splineBasis = c()
      splineBasis.predict = c()
      Z = NA
      Z.predict = NA
      X.pred = cbind(X.predict[rep(1:nrow(X.predict), each = length(daysSinceJan.predict)),])
    }
  }
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


# There's no reinfection process for Ebola, but we still need to provide dummy
# values for the reinfection terms. This will be changed (along with most of 
# the R level API) Dummy covariate matrix:
X_p_rs = matrix(0)

# Dummy value for reinfection params
beta_p_rs = rep(0, ncol(X_p_rs))
# Dummy value for reinfection params prior precision
betaPrsPriorPrecision = 0.5



# Declare prior parameters for the E to I and I to R probabilities. 

# Latent period between 1 and 21 days
prior_gamma_ei = 1/5
prior_gamma_ir = 1/7
transitionEffectiveSampleSizes = 10000

# Declare prior parameters for the overdispersion precision
priorAlpha_phi = 10000
priorBeta_phi = 1500

# Declare prior precision for exposure model paramters
betaPriorPrecision = 1

# Declare a function which can come up with several different starting values 
# for the model parameters. This will allow us to assess convergence. 
proposeParameters = function(seedVal, chainNumber)
{
  set.seed(seedVal) 
  
  # 2 to 21 day incubation period according to who
  p_ei = 0.25 + rnorm(1, 0, 0.02) 
  # Up to 7 weeks even after recovery
  p_ir = 0.14 + rnorm(1, 0, 0.01) 
  gamma_ei=-log(1-p_ei)
  gamma_ir=-log(1-p_ir)
  
  # Starting value for exposure regression parameters
  beta = rep(0, ncol(X) + ifelse(all(is.na(Z)), 0, ncol(Z)))
  beta[1:(ncol(X))] = rnorm(ncol(X), -3, 3)
  #beta[1] = 2.5 + rnorm(1,0,0.5)
  
  phi = 0.01 # Overdispersion precision
  
  outFileName = paste("./chain_output_ebola_", chainNumber ,".txt", sep = "")
  
  # Make a crude guess as to the true compartments:
  # S_star, E_star, R_star, and thus S,E,I and R
  DataModel = buildDataModel(I_star, type = "overdispersion", 
                             params=c(priorAlpha_phi,priorBeta_phi))
  ExposureModel = buildExposureModel(X, Z, beta, betaPriorPrecision, 0, offsets, nTpt = nrow(I_star))
  ReinfectionModel = buildReinfectionModel("SEIR")
  SamplingControl = buildSamplingControl(iterationStride=iterationStride,
                                         sliceWidths=c(1, # S_star
                                                       1, # E_star
                                                       1,  # R_star
                                                       1,  # S_0
                                                       1,  # I_0
                                                       0.05,  # beta
                                                       0.0,  # beta_p_rs, fixed in this case
                                                       0.01, # rho
                                                       0.01, # gamma_ei
                                                       0.01,  # gamma_ir
                                                       0.01)) # phi)
  
  InitContainer = buildInitialValueContainer(I_star, N, 
                                             S0 = N[1,]-I_star[1,] - 2*I0,
                                             E0 = I0,
                                             I0 = I0)
  DistanceModel = buildDistanceModel(dmList, priorAlpha = 1, priorBeta = 500)
  TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-prior_gamma_ei), 1-exp(-prior_gamma_ir),
                                                            transitionEffectiveSampleSizes, transitionEffectiveSampleSizes)
  return(list(DataModel=DataModel,
              ExposureModel=ExposureModel,
              ReinfectionModel=ReinfectionModel,
              SamplingControl=SamplingControl,
              InitContainer=InitContainer,
              DistanceModel=DistanceModel,
              TransitionPriors=TransitionPriors,
              outFileName=outFileName))
}


paramsList = list(list("estimateR0"=FALSE, "traceCompartments"=TRUE, "seedVal"=133,"chainNumber"=4),
                  list("estimateR0"=TRUE, "traceCompartments"=FALSE, "seedVal"=1224,"chainNumber"=5),
                  list("estimateR0"=FALSE,"traceCompartments"=FALSE, "seedVal"=12325,"chainNumber"=6))

buildAndBurnInModel = function(params)
{
  # Need this call because this function executes in a node
  library(spatialSEIR)
  # save proposal and params to node workspace
  proposal <<- proposeParameters(params[["seedVal"]], params[["chainNumber"]])
  params <<- params
  SEIRmodel =  buildSEIRModel(proposal$outFileName,
                              proposal$DataModel,
                              proposal$ExposureModel,
                              proposal$ReinfectionModel,
                              proposal$DistanceModel,
                              proposal$TransitionPriors,
                              proposal$InitContainer,
                              proposal$SamplingControl)
  
  SEIRmodel$setRandomSeed(params[["seedVal"]])
  # Save model object to node workspace
  localModelObject <<- SEIRmodel
  
  # Do we need to keep track of compartment values for prediction? 
  # No sense doing this for all of the chains.
  if (params[["traceCompartments"]])
  {
    SEIRmodel$setTrace(0) #Guinea 
    SEIRmodel$setTrace(1) #Liberia
    SEIRmodel$setTrace(2) #Sierra Leone
  }
  
  # Make a helper function to run each chain, as well as update the metropolis 
  # tuning parameters. 
  runSimulation = function(modelObject,
                           numBatches=500, 
                           batchSize=20, 
                           targetAcceptanceRatio=0.2,
                           tolerance=0.05,
                           proportionChange = 0.1
  )
  {
    for (batch in 1:numBatches)
    {
      modelObject$simulate(batchSize)
      modelObject$updateSamplingParameters(targetAcceptanceRatio, 
                                           tolerance, 
                                           proportionChange)
    }
  }
  
  # Burn in tuning parameters
  runSimulation(SEIRmodel, numBatches = numBurnInBatches)
  runSimulation(SEIRmodel, batchSize = 100, numBatches = numBurnInBatches)
  
  SEIRmodel$compartmentSamplingMode = 17
  SEIRmodel$performHybridStep = 100
  if (modelDF > 0)
  {
    SEIRmodel$useDecorrelation = 10
  }
}


finishSimulation = function(iterationNumber)
{
  dat = read.csv(proposal$outFileName)
  
  ## Do we need to estimate R0 for this chain?
  if (params[["estimateR0"]])
  {  
    R0 = array(0, dim = c(nrow(I_star), ncol(I_star), extraR0Iterations))
    effectiveR0 = array(0, dim = c(nrow(I_star), ncol(I_star), extraR0Iterations))
    empiricalR0 = array(0, dim = c(nrow(I_star), ncol(I_star), extraR0Iterations))
    for (i in 1:extraR0Iterations)
    {
      localModelObject$simulate(extraR0BatchSize)
      for (j in 0:(nrow(I_star)-1))
      {
        R0[j,,i] = localModelObject$estimateR0(j)
        effectiveR0[j,,i] = localModelObject$estimateEffectiveR0(j)
        empiricalR0[j,,i] = apply(localModelObject$getIntegratedGenerationMatrix(j), 1, sum)
      }
    }
    
    R0Mean = apply(R0, 1:2, mean)
    R0LB = apply(R0, 1:2, quantile, probs = 0.05)
    R0UB = apply(R0, 1:2, quantile, probs = 0.95)
    effectiveR0Mean = apply(effectiveR0, 1:2, mean)
    effectiveR0LB = apply(effectiveR0, 1:2, quantile, probs = 0.05)
    effectiveR0UB = apply(effectiveR0, 1:2, quantile, probs = 0.95)
    empiricalR0Mean = apply(empiricalR0, 1:2, mean)
    empiricalR0LB = apply(empiricalR0, 1:2, quantile, probs = 0.05)
    empiricalR0UB = apply(empiricalR0, 1:2, quantile, probs = 0.95)
    orig.R0 = R0
    R0 = list("R0" = list("mean"=R0Mean, "LB" = R0LB, "UB" = R0UB),
              "effectiveR0" = list("mean"=effectiveR0Mean, "LB" = effectiveR0LB, 
                                   "UB" = effectiveR0UB),
              "empiricalR0" = list("mean"=empiricalR0Mean, "LB" = empiricalR0LB, 
                                   "UB" = empiricalR0UB))
  } else
  {
    R0 = NULL
    orig.R0 = NULL
  }  
  
  return(list("chainOutput" = dat, "R0" = R0, "rawSamples" = orig.R0))
}


cl = makeCluster(3, outfile = "err.txt")
clusterExport(cl, c( "offsets",
                     "X",
                     "Z",
                     "I0",
                     "X_p_rs",
                     "prior_gamma_ei",
                     "prior_gamma_ir",
                     "transitionEffectiveSampleSizes",
                     "priorAlpha_phi",
                     "priorBeta_phi",
                     "betaPriorPrecision",
                     "beta_p_rs",
                     "betaPrsPriorPrecision",
                     "N",
                     "dmList",
                     "iterationStride",
                     "proposeParameters",
                     "generateCompartmentProposal",
                     "extraR0Iterations",
                     "extraR0BatchSize",
                     "I_star",
                     "modelDF",
                     "numBurnInBatches"), envir=environment())

additionalIterations = function(params)
{
  
  N = params[[1]]
  batchSize = params[[2]]
  targetRatio = params[[3]]
  targetWidth=params[[4]]
  proportionChange = params[[5]]
  updateParameters = params[[6]]
  for (i in 1:(N/batchSize))
  {
    localModelObject$simulate(batchSize)
    if (updateParameters)
    {
      localModelObject$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
    }
  }
}

iterationParams = list(convergenceSampleSize, convergenceBatchSize,
                       targetAcceptanceRatio=0.2,   
                       tolerance=0.05,
                       proportionChange = 0.1,
                       updateSamplingParams = FALSE)
iterationParams = list(iterationParams, iterationParams, iterationParams)
fileNames = paste(paste("chain_output_ebola_", c(paramsList[[1]]$chainNumber,
                                                 paramsList[[2]]$chainNumber,
                                                 paramsList[[3]]$chainNumber), sep = ""), ".txt", sep = "")                      
chains = parLapply(cl, paramsList, buildAndBurnInModel)

conv = FALSE
while (!conv)
{
  cat(paste("Not converged, adding iterations. Total so far: ", totalSamples, 
            "\n", sep =""))
  parLapply(cl, iterationParams, additionalIterations)
  totalSamples = totalSamples + convergenceSampleSize
  conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3], maxVal = convergenceCriterion)
  conv = (conv && (minimumSamples < totalSamples))
}
cat("Chains converged, finishing up...\n")

cleanUpParamsList = list(1,2,3)
chains = parLapply(cl, cleanUpParamsList, finishSimulation)
stopCluster(cl)


chain1 = chains[[1]]$chainOutput 
chain2 = chains[[2]]$chainOutput 
chain3 = chains[[3]]$chainOutput 

nbeta = ncol(X) + ifelse(class(Z) == "matrix", ncol(Z), 0)
c1 = chain1[floor(nrow(chain1)/2):nrow(chain1),c(1:(nbeta + length(dmList)),
                                                 (nbeta+length(dmList)+1):
                                                   (nbeta+length(dmList)+3))]
c2 = chain2[floor(nrow(chain2)/2):nrow(chain2),c(1:(nbeta + length(dmList)),
                                                 (nbeta+length(dmList)+1):
                                                   (nbeta+length(dmList)+3))]
c3 = chain3[floor(nrow(chain3)/2):nrow(chain3),c(1:(nbeta + length(dmList)),
                                                 (nbeta+length(dmList)+1):
                                                   (nbeta+length(dmList)+3))]

c1$gamma_ei = 1-exp(-c1$gamma_ei)
c1$gamma_ir = 1-exp(-c1$gamma_ir)
c2$gamma_ei = 1-exp(-c2$gamma_ei)
c2$gamma_ir = 1-exp(-c2$gamma_ir)
c3$gamma_ei = 1-exp(-c3$gamma_ei)
c3$gamma_ir = 1-exp(-c3$gamma_ir)
{
  if (modelDF != 0)
  {
    colnames(c1) = c("Guinea Intercept", "Liberia Intercept", 
                     "Sierra Leone Intercept", 
                     paste("Time component ", (1:(nbeta-ncol(X))), sep = ""),
                     "Overdispersion Precision",
                     paste("Spatial Dependence Parameter", 1:length(dmList)),
                     "E to I probability", 
                     "I to R probability")
  }
  else
  {
    colnames(c1) = c("Guinea Intercept", "Liberia Intercept", 
                     "Sierra Leone Intercept",
                     "Overdispersion Precision",
                     paste("Spatial Dependence Parameter", 1:length(dmList)),
                     "E to I probability", 
                     "I to R probability")
  }
}
colnames(c2) = colnames(c1)
colnames(c3) = colnames(c1)

mcl = mcmc.list(as.mcmc(c1), 
                as.mcmc(c2),
                as.mcmc(c3))

R0_list = chains[[2]]$R0


getMeanAndCI = function(loc,tpt,baseStr="I_")
{
  vec = chain1[[paste(baseStr, loc, "_", tpt, sep = "")]]
  vec = vec[floor(length(vec)/2):length(vec)]
  return(c(mean(vec), quantile(vec, probs = c(0.05, 0.95))))
}

Guinea.I.Est = sapply(0:(nrow(I_star)- 1), getMeanAndCI, loc=0)
Liberia.I.Est = sapply(0:(nrow(I_star)- 1), getMeanAndCI, loc=1)
SierraLeone.I.Est = sapply(0:(nrow(I_star)- 1), getMeanAndCI, loc=2)

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
  dataRow = chain1[i,]
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

preds = lapply((nrow(chain1) - floor(nrow(chain1)/2)):
                 nrow(chain1), predict.i)


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


Guinea.LB = apply(Guinea.Pred, 2, quantile, probs = c(0.05))
Guinea.UB = apply(Guinea.Pred, 2, quantile, probs = c(0.95))

Liberia.LB = apply(Liberia.Pred, 2, quantile, probs = c(0.05))
Liberia.UB = apply(Liberia.Pred, 2, quantile, probs = c(0.95))

SierraLeone.LB = apply(SierraLeone.Pred, 2, quantile, probs = c(0.05))
SierraLeone.UB = apply(SierraLeone.Pred, 2, quantile, probs = c(0.95))

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

predlist=list(Guinea.mean=Guinea.mean,
              Liberia.mean=Liberia.mean,
              SierraLeone.mean=SierraLeone.mean,
              Guinea.LB=Guinea.LB,
              Liberia.LB=Liberia.LB,
              SierraLeone.LB=SierraLeone.LB,
              Guinea.UB=Guinea.UB,
              Liberia.UB=Liberia.UB,
              SierraLeone.UB=SierraLeone.UB)

save("chains", "predlist", file="./chainOutput.Robj")

return(list(chainOutput=chains,
            predictions=predlist))
}
