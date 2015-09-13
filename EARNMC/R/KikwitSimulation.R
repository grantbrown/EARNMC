generateSEIRData = function(seed, 
                            nTpt,
                            S0,
                            E0,
                            I0,
                            R0,
                            timeIndex,
                            beta_SE,
                            beta_RS,
                            distMatList,
                            rho,
                            X,
                            X_RS,
                            gamma_ei,
                            gamma_ir,
                            effectiveTransitionSampleSize
                            )
{
    set.seed(seed)
    nLoc = length(S0)   
    if (nrow(X) != length(S0)*nTpt)
    {
        stop(paste("Arguments are inconsistent with ", nTpt, 
                   " time points.\n", sep = ""))
    }


    bigX_SE = X

    # Calculate temporal offsets. 
    uncumulate = function(x)
    {
        out = c(x[2:length(x)]-x[1:(length(x)-1)])
        ifelse(out >= 0, out, 0)
    }
    if (length(timeIndex) - nrow(X_RS) != 1)
    {
        stop(paste("ERROR: timeIndex should be one longer than the number of",
                   "time points in the study."))
    }
    offsets = uncumulate(timeIndex)[1:nTpt]

    # Calculate linear predictor for the intensity process. 
    eta_SE = exp(matrix((bigX_SE %*% beta_SE), ncol = nLoc))

    # Calculate linear predictor for reinfection process.
    X_RS = X_RS[1:nTpt,,drop=FALSE]
    eta_RS = exp(X_RS %*% beta_RS)

    # Calculate p_EI and p_IR
    p_EI = 1-exp(-gamma_ei*offsets) 
    p_IR = 1-exp(-gamma_ir*offsets)


    # p_RS
    p_RS = 1-exp(-offsets*(eta_RS))

    # Allocate compartments
    S_star = E_star = I_star = R_star = S = E = I = R = matrix(0, nrow = nTpt, 
                                                               ncol = nLoc)

    # Declare N
    N = matrix(S0+E0+I0+R0, ncol = nLoc, nrow = length(offsets), byrow=TRUE)

    # Run Simulation 
    offsetMatrix = matrix(offsets, nrow = length(offsets), ncol = nLoc)
    S[1,] = S0
    E[1,] = E0
    I[1,] = I0
    R[1,] = R0

    p_SE = eta_SE
    for (i in 1:(nTpt))
    {
        # Calculate p_SE
        if (nLoc > 1)
        {
            # Case 1: spatial
            p_SE[i,] = p_SE[i,]*(I[i,]/N[i,]) 
            for (j in 1:length(distMatList))
            {
               p_SE[i,] = p_SE[i,] + rho[[j]]*distMatList[[j]] %*% 
                          (eta_SE[i,] * (I[i,]/N[i,]))
            }
            p_SE[i,] = 1-exp(-offsetMatrix[i,]*p_SE[i,])
        }
        else
        {
            # Case 2: single location
            p_SE[i,] = (p_SE[i,]*(I[i,]/N[i,]))
            p_SE[i,] = 1-exp(-offsetMatrix[i,]*p_SE[i,])
        }

        S_star[i,] = rbinom(rep(1, ncol(S)), R[i,], rep(p_RS[i], ncol(S)))
        E_star[i,] = rbinom(rep(1, ncol(S)), S[i,], p_SE[i,])
        I_star[i,] = rbinom(rep(1, ncol(S)), E[i,], rep(p_EI[i], ncol(S)))
        R_star[i,] = rbinom(rep(1, ncol(S)), I[i,], rep(p_IR[i], ncol(S)))

        if (i < nTpt)
        {
            S[i+1,] = S[i,] + S_star[i,] - E_star[i,] 
            E[i+1,] = E[i,] + E_star[i,] - I_star[i,]
            I[i+1,] = I[i,] + I_star[i,] - R_star[i,]
            R[i+1,] = R[i,] + R_star[i,] - S_star[i,]
        }
    }

    return(list(S_star=S_star,
                E_star=E_star,
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
                N =N, 
                X=bigX_SE,
                X_RS=X_RS,
                beta_SE=beta_SE,
                beta_RS=beta_RS,
                gamma_ei=gamma_ei,
                gamma_ir=gamma_ir,
                distMatList=distMatList,
                p_SE=p_SE,
                p_EI=p_EI,
                p_IR=p_IR,
                p_RS=p_RS,
                effectiveTransitionSampleSize=effectiveTransitionSampleSize,
        		rho=rho,
                nTpt=nTpt,
                nLoc=nLoc
                ))
}

generateSimpleInterventionData = function(seed, 
                                  population=5363500, 
                                  ThrowAwayTpt=0,
                                  beta_SE=c(-1.4,-0.055),
                                  beta_RS=c(-100000))
{    
    set.seed(seed)
    data(Kikwit95) 
    E0 = 0  
    I0 = 1 
    R0 = 0
    S0 = population - R0 - E0 - I0

    maxTpt = 150
    nTpt = maxTpt - ThrowAwayTpt
    


    interventionDate = as.Date("05-09-1995", "%m-%d-%Y")
    hasIntervention = (Kikwit95$Date > interventionDate)
    X = cbind(1, hasIntervention*(Kikwit95$Date - interventionDate))[1:nTpt,]

    X_RS = matrix(1, nrow = nTpt)



    distMatList = -1   
    rho = -1
    gamma_ei = 1/5
    gamma_ir = 1/7
    
    effectiveTransitionSampleSize = 1000
    timeIndex = 0:maxTpt

    out = generateSEIRData(seed + 1, 
                       nTpt,
                       S0,
                       E0,
                       I0,
                       R0,
                       timeIndex,
                       beta_SE,
                       beta_RS,
                       distMatList,
                       rho,
                       X,
                       X_RS,
                       gamma_ei,
                       gamma_ir,
                       effectiveTransitionSampleSize
                       )

    if (sum(out$I == 0) > 100 || max(out$I) < 15)
    {
        cat("Rejecting simulated data - epidemic died out too soon.\n")
        return(generateSimpleInterventionData(seed + rpois(1, 1000), 
                                              ThrowAwayTpt=ThrowAwayTpt, 
                                              beta_SE=beta_SE))
    }
    out
}


simpleInterventionSimulationKernel = function(cl, genSeed, fitSeeds, 
                                              underspecified, beta_SE)
{
    #TODO: Vary starting linear predictor parameters on each iteration 
    simResults = generateSimpleInterventionData(genSeed, 
                                                ThrowAwayTpt=0, 
                                                beta_SE=beta_SE)


    fileNames = c("sim2_1.txt", "sim2_2.txt", "sim2_3.txt")
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], 
                           simResults, underspecified),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], 
                           simResults, underspecified),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], 
                           simResults, underspecified))

    trueVals = parLapply(cl, paramsList, 
                         buildSimpleInterventionSimulationKernel)

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
        #localModelObject$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
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

buildSimpleInterventionSimulationKernel = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]
    underspecified = params[[4]]

    set.seed(seed)


    DataModel = buildDataModel(simResults$I_star, type = "overdispersion", 
                               phi = 1)

    priorBetaIntercept = log(mean(-log(1-(simResults$I_star/(simResults$N))))) 
    if (underspecified)
    {
        ExposureModel = buildExposureModel(simResults$X[,1,drop=FALSE], nTpt = simResults$nTpt,
                                           nLoc = simResults$nLoc,
                                           beta = c(priorBetaIntercept),
                                           betaPriorPrecision = 1)
    }
    else
    {
        ExposureModel = buildExposureModel(simResults$X, nTpt = simResults$nTpt,
                                           nLoc = simResults$nLoc,
                                           beta = c(priorBetaIntercept, 
                                     rep(0,((length(simResults$beta_SE))-1))), 
                                     betaPriorPrecision = 1)
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
                                                              simResults$effectiveTransitionSampleSize,
                                                              simResults$effectiveTransitionSampleSize) 

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
    res$useDecorrelation = 50
    res$performHybridStep = 50

    # Store the model object in the global namespace of the node 
    # - can't pass these between sessions
    localModelObject <<- res
    return(list("model"=res,
                "fileName"=outFileName))    
}



runSimpleInterventionSimulation = function(cellIterations = 50, 
                                           genSeed=123123, 
                                           fitSeeds=c(812123,12301,5923),
                                           beta_SE=c(-1.4,-0.055))
{                     
    cl = makeCluster(3, outfile = "err.txt")
    print("Cluster Created")
    clusterExport(cl, c("buildSimpleInterventionSimulationKernel"))
    print("Variables Exported.") 
 

    for (underspec in c(TRUE, FALSE))
    {

        f = function(genSeed)
        {
            simpleInterventionSimulationKernel(cl, genSeed, 
                                               fitSeeds + genSeed,
                                               underspec, beta_SE)
        }

        itrSeeds = genSeed + seq(1, cellIterations)
        i = 1
        for (itrSeed in itrSeeds){
            fname = paste("./sim2_results_Uspec_", underspec, "_", i, 
                                    ".Rda.bz2", sep="")
            if (!file.exists(fname))
            {
                result = f(itrSeed)
                save(result, file=fname, 
                     compress="bzip2")
                i=i+1
            }
            else
            {
                cat(paste("Exists: ", fname, "\n"))
                i=i+1
            }
        }    
    }

    print("Results obtained")
    stopCluster(cl)
    TRUE
}
