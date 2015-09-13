generateABSEIRData = function(seed, 
                              nTpt,
                              S0,
                              E0,
                              I0,
                              R0,
                              timeIndex,
                              beta_SE,
                              X,
                              gamma_ei,
                              gamma_ir
                              )
{
    set.seed(seed)
    nLoc = length(S0)   
    if (nrow(X) != length(S0)*nTpt)
    {
        stop(paste("Arguments are inconsistent with ", nTpt, 
                   " time points.\n", sep = ""))
    }
    # Naiive ABM is too slow for spatial models
    if (length(S0) != 1)
    {
        stop(paste("Spatial simulations not currently supported for agent", 
                   " based approach.", sep = ""))
    }

    bigX_SE = X

    # Calculate temporal offsets. 
    uncumulate = function(x)
    {
        out = c(x[2:length(x)]-x[1:(length(x)-1)])
        ifelse(out >= 0, out, 0)
    }
    offsets = uncumulate(timeIndex)[1:nTpt]

    # Calculate linear predictor for the intensity process. 
    eta_SE = exp(matrix((bigX_SE %*% beta_SE), ncol = nLoc))

    # Calculate p_EI and p_IR
    p_EI = 1-exp(-gamma_ei*offsets) 
    p_IR = 1-exp(-gamma_ir*offsets)


    # Declare N
    N = matrix(S0+E0+I0+R0, ncol = nLoc, nrow = length(offsets), byrow=TRUE)

    # Allocate compartments
    E_star = I_star = R_star = S = E = I = R = matrix(0, nrow = nTpt, 
                                                               ncol = N[1])
    # Row infects column
    infectionEvents = array(0, dim = c(N[1], N[1], nTpt))

    # Run Simulation 
    offsetMatrix = matrix(offsets, nrow = length(offsets), ncol = nLoc)
    idx = sample(1:ncol(E_star))
    
    S[1, idx[1:(S0[1])]] = 1
    E[1, idx[(S0[1] + 1):(S0[1] + E0[1])]] = 1
    I[1, idx[(S0[1] + E0[1] + 1):(S0[1] + E0[1] +I0[1])]] = 1
    # R[1,] = 0

    p_SE = eta_SE
    for (i in 1:(nTpt))
    {
        cat(paste(i, "\n"));
        # Includes contact and infection events
        p_SE[i,] = 1-exp(-offsetMatrix[i,]*p_SE[i,])
        n= N[1,1]
        cmat = Matrix(0, nrow = n, ncol = n)    
        n.tri = n*(n-1)/2 
        cmat[upper.tri(cmat)] = rbinom(rep(1, n.tri), rep(1, n.tri), rep(p_SE[i,1], n.tri))
        for (k in 1:nrow(cmat))
        {
            for (l in which(cmat[k,]==1))
            {
 
                e.1 = I[i,l] && S[i,k] && (!E_star[i,k])
                e.2 = S[i,l] && I[i,k] && (!E_star[i,l])
                infectionEvents[l,k,i] = e.1
                infectionEvents[k,l,i] = e.2
                E_star[i,k] = (e.1 || E_star[i,k])
                E_star[i,l] = (e.2 || E_star[i,l])
            }
        }
        I_star[i,] = rbinom(rep(1, ncol(S)), E[i,], rep(p_EI[i], ncol(S)))
        R_star[i,] = rbinom(rep(1, ncol(S)), I[i,], rep(p_IR[i], ncol(S)))

        if (i < nTpt)
        {
            S[i+1,] = S[i,] - E_star[i,] 
            E[i+1,] = E[i,] + E_star[i,] - I_star[i,]
            I[i+1,] = I[i,] + I_star[i,] - R_star[i,]
            R[i+1,] = R[i,] + R_star[i,]
        }
    }

    return(list(E_star=E_star,
                I_star=I_star,
                R_star=R_star,
                S=S,
                I=I,
                E=E,
                R=R,
                S0=S0,
                E0=E0,
                I0=I0,
                R0=R0,
                infectionEvents=infectionEvents,
                N =N, 
                X=bigX_SE,
                beta_SE=beta_SE,
                gamma_ei=gamma_ei,
                gamma_ir=gamma_ir,
                p_SE=p_SE,
                p_EI=p_EI,
                p_IR=p_IR,
                nTpt=nTpt
                ))
}


abmInterventionSimulationKernel = function(cl, genSeed, fitSeeds, 
                                              underspecified, beta_SE)
{
  #TODO: Vary starting linear predictor parameters on each iteration 
  
  nTpt = 150
  enableJIT(3)
  fc = cmpfun(generateABSEIRData, options = list(optimize = 3))
  X = cbind(1, cumsum(1:nTpt > floor(nTpt/2)))
  X[,2] = X[,2]/max(X[,2])
  simResults = fc(genSeed, nTpt, 500, 10, 10, 0, 0:nTpt, beta_SE, X, 1/5, 1/7)
  
  fileNames = paste(paste(paste("sim3_", genSeed, "_", sep = ""), 1:3, sep = ""), ".txt", sep = "")

  paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], 
                         simResults, underspecified),
                    list(seed=fitSeeds[2], outFileName = fileNames[2], 
                         simResults, underspecified),
                    list(seed=fitSeeds[3], outFileName = fileNames[3], 
                         simResults, underspecified))
  
  trueVals = parLapply(cl, paramsList, 
                       buildAbmInterventionSimulationKernel)
  
  iterationParams = list(list(200000, 1000, 0.2, 0.05, 0.1),
                         list(200000, 1000, 0.2, 0.05, 0.1),
                         list(200000, 1000, 0.2, 0.05, 0.1))
  
  additionalIterations = function(params)
  {
    
    N = params[[1]]
    batchSize = params[[2]]
    targetRatio = params[[3]]
    targetWidth=params[[4]]
    proportionChange = params[[5]]
    for (i in 1:(N/batchSize))
    {
      localModelObject$simulate(batchSize)
      localModelObject$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
    }
  }
  
  
  conv = FALSE
  while (!conv)
  {
    cat("Not converged, adding iterations...\n")
    parLapply(cl, iterationParams, additionalIterations)
    conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3], 
                            maxVal=1.02)
  }
  
  dat1=read.csv(fileNames[1])
  dat2=read.csv(fileNames[2])
  dat3=read.csv(fileNames[3])
  list(dat1=dat1,dat2=dat2,dat3=dat3,simResults=simResults)
}

buildAbmInterventionSimulationKernel = function(params) 
{
  library(spatialSEIR)
  seed = params[[1]]
  outFileName = params[[2]]
  simResults = params[[3]]
  underspecified = params[[4]]
  
  set.seed(seed)
  
  simResults$I_star = matrix(apply(simResults$I_star, 1, sum), ncol = 1)
  #DataModel = buildDataModel(simResults$I_star, type = "overdispersion", 
  #                           phi = 1)
  DataModel = buildDataModel(simResults$I_star, type = "identity")

  
  priorBetaIntercept = log(mean(-log(1-(simResults$I_star/(simResults$N))))) 
  if (underspecified)
  {
    ExposureModel = buildExposureModel(simResults$X[,1,drop=FALSE], nTpt = simResults$nTpt,
                                       nLoc = 1,
                                       beta = c(priorBetaIntercept),
                                       betaPriorPrecision = 0.1)
  }
  else
  {
    ExposureModel = buildExposureModel(simResults$X, nTpt = simResults$nTpt,
                                       nLoc = 1,
                                       beta = c(priorBetaIntercept, 
                                                rep(0,((length(simResults$beta_SE))-1))), 
                                       betaPriorPrecision = 0.1)
  }
  
  ReinfectionModel = buildReinfectionModel("SEIR")
  SamplingControl = buildSamplingControl(iterationStride=1000,
                                         sliceWidths = c(0.26,  # S_star
                                                         0.1,  # E_star
                                                         0.15, # I_star
                                                         0.22, # S0
                                                         0.24, # I0
                                                         0.8, # beta
                                                         0.2, # betaPrs
                                                         0.015, # rho
                                                         0.01, # gamma_ei
                                                         0.01, # gamma_ir
                                                         0.01 # phi
                                         ))
  DistanceModel = buildDistanceModel(list(matrix(0)))
  TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-simResults$gamma_ei), 
                                                            1-exp(-simResults$gamma_ir), 
                                                            1000,
                                                            1000) 
  
  I0 = max(simResults$I_star[1:10], simResults$I_star[2], 2)
  E0 = I0
  S0 = simResults$N[1] - I0 - E0
  InitContainer = buildInitialValueContainer(simResults$I_star, simResults$N, 
                                             S0 = S0, I0 = I0, E0 = E0, 
                                             reinfection=FALSE, 
                                             dataType = "I_star")
  res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,
                       DistanceModel,TransitionPriors,
                       InitContainer,SamplingControl)
  
  res$setRandomSeed(seed)
  res$setTrace(0)
  # Burn in tuning parameters
  for (i in 1:(500))
  {
    res$simulate(10)
    res$updateSamplingParameters(0.2, 0.05, 0.01)
  }
  for (i in 1:(50))
  {
    res$simulate(100)
    res$updateSamplingParameters(0.2, 0.05, 0.01)
  }
  res$compartmentSamplingMode = 17
  res$useDecorrelation = 11
  res$performHybridStep = 13
  
  # Store the model object in the global namespace of the node 
  # - can't pass these between sessions
  localModelObject <<- res
  return(list("model"=res,
              "fileName"=outFileName))    
}



runAbmInterventionSimulation = function(cellIterations = 50, 
                                           genSeed=123123, 
                                           fitSeeds=c(123123,234234,345345),
                                           beta_SE=c(-8, -0.8))
{                     
  seeds = genSeed + 10*seq(1, cellIterations)
  params = lapply(seeds, function(x){list(seed = x, 
                                          beta_SE=beta_SE)})
  main.cluster = makeCluster(4)
  clusterExport(main.cluster, c("fitSeeds", "beta_SE"), envir = environment())
  outer.loop = function(paramVal){
    library(EARNMC)    
    genSeed = paramVal[["seed"]]
    
    cl = makeCluster(3)
    print("Cluster Created")
    clusterExport(cl, c("buildAbmInterventionSimulationKernel"))
    print("Variables Exported.") 
    for (underspec in c(TRUE, FALSE))
    {
      
      f = function(genSeed)
      {
        abmInterventionSimulationKernel(cl, genSeed, 
                                        fitSeeds + genSeed,
                                        underspec, beta_SE)
      }
      

      fname = paste("./sim3_results_Uspec_", underspec, "_", genSeed, 
                    ".Rda.bz2", sep="")
      if (!file.exists(fname))
      {
        result = f(genSeed)
        save(result, file=fname, 
             compress="bzip2")
      }
      else
      {
        cat(paste("Exists: ", fname, "\n"))
      }
    }     
    print("Results obtained")
    stopCluster(cl)
  }
  parLapplyLB(main.cluster, params, outer.loop)	
  stopCluster(main.cluster)
}




