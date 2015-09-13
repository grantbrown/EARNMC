generateSpatialMisspecificationData = function(seed, 
                                  population=floor(c(0.4, 
                                                     0.2, 
                                                     0.4)*5363500), 
                                  ThrowAwayTpt=0,
                                  beta_SE=c(-1.1,-0.055),
                                  beta_RS=c(-100000))
{    
    set.seed(seed)
    data(Kikwit95) 
    # Spatial Locations: s1, s2, s3, s4, s5, s6, s7, s8
    # Clusters: (s1, s2, s3), (s4, s5), (s6, s7, s8)
    E0 = rep(0, 3)  
    I0 = c(0,0,2) 
    R0 = rep(0, 3)
    S0 = population - R0 - E0 - I0

    maxTpt = 150
    nTpt = maxTpt - ThrowAwayTpt
    nLoc = length(I0)
    nLoc_uspec = 3
    


    interventionDate = as.Date("05-09-1995", "%m-%d-%Y")
    hasIntervention = (Kikwit95$Date > interventionDate)
    X = cbind(1, hasIntervention*(Kikwit95$Date - interventionDate))  
    X =  X[rep(1:nTpt, nLoc),]
    X_underspec = X
    
    X_RS = matrix(1, nrow = maxTpt)

    # s1 <-> s2
    DM1 = diag(nLoc)*0 
    DM1[1,2] = 1
    DM1[2,1] = 1
    
    # s2 <-> s3
    DM2 = diag(nLoc)*0
    DM2[2,3] = 1
    DM2[3,2] = 1
    
    #DM1_uspec = matrix(c(0,1,0,
    #                     1,0,0,
    #                     0,0,0), nrow = 3, byrow = TRUE)  
    #DM2_uspec = matrix(c(0,0,1,
    #                     0,0,0,
    #                     1,0,0), nrow = 3, byrow = TRUE)
    #DM3_uspec = matrix(c(0,0,0,
    #                     0,0,1,
    #                     0,1,0), nrow = 3, byrow = TRUE)
    DM_uspec = (DM1 + DM2)/2
    
    #standardizeDM = function(DM)
    #{
    #  sm = apply(DM, 1, sum)
    #  1.0*DM/(matrix(ifelse(sm == 0, 1, sm), nrow = nLoc, ncol = nLoc))
    #}
    #distMatList = lapply(list(DM1, DM2), standardizeDM)
    distMatList = list(DM1, DM2)
    distMatList_uspec = list(DM_uspec)
    #distMatList_uspec = list(DM1_uspec, 
    #                         DM2_uspec,
    #                         DM3_uspec)
    rho = c(0.05,0.1)
    gamma_ei = 1/5
    gamma_ir = 1/7
    
    effectiveTransitionSampleSize = 10000
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

    if (any(apply(out$I, 2, sum) < 5))
    {
        cat("Rejecting simulated data - epidemic died out too soon.\n")
        return(generateSpatialMisspecificationData(seed + rpois(1, 1000), 
                                              ThrowAwayTpt=ThrowAwayTpt, 
                                              beta_SE=beta_SE))
    }
    #out$I_star_uspec = cbind(out$I_star[,1] + out$I_star[,2] + out$I_star[,3])
    out$I_star_uspec = out$I_star
    out$X_underspec = X_underspec
    out$distMatList_uspec = distMatList_uspec
    #out$N_uspec = cbind(out$N[,1] + out$N[,2] + out$N[,3])
    out$N_uspec = out$N
    out
}


spatialMisspecificationSimulationKernel = function(cl, genSeed, fitSeeds, 
                                              underspecified,
                                              scratchLoc)
{
    simResults = generateSpatialMisspecificationData(genSeed)

    fileNames = paste(scratchLoc, 
                      paste(paste("sim2_", genSeed, sep = ""), "_", 1:3, ".txt", sep = ""),
                      sep = "")
    paramsList = list(list(seed=fitSeeds[1], outFileName = fileNames[1], 
                           simResults, underspecified),
                      list(seed=fitSeeds[2], outFileName = fileNames[2], 
                           simResults, underspecified),
                      list(seed=fitSeeds[3], outFileName = fileNames[3], 
                           simResults, underspecified))

    trueVals = parLapply(cl, paramsList, 
                         buildSpatialMisspecificationSimulationKernel)

    iterationParams = list(list(20000, 0.1, 0.01, 0.05),
                           list(20000, 0.1, 0.01, 0.05),
                           list(20000, 0.1, 0.01, 0.05))
    additionalIterations = function(params)
    {
      N = params[[1]]
      targetRatio = params[[2]]
      targetWidth=params[[3]]
      proportionChange = params[[4]]
      localModelObject$simulate(N)
      #localModelObject$updateSamplingParameters(targetRatio, targetWidth, proportionChange)
    }
    
    
    
    conv = FALSE
    itrs = 0
    while (!conv || itrs < 400000)
    {
        cat("Not converged, adding iterations...\n")
        parLapply(cl, iterationParams, additionalIterations)
        conv = checkConvergence(fileNames[1], fileNames[2], fileNames[3], 
                                maxVal=1.02)
        itrs = itrs + 20000
    }

    dat1=read.csv(fileNames[1])
    dat2=read.csv(fileNames[2])
    dat3=read.csv(fileNames[3])
    list(dat1=dat1,dat2=dat2,dat3=dat3,simResults=simResults)
}

buildSpatialMisspecificationSimulationKernel = function(params) 
{
    library(spatialSEIR)
    seed = params[[1]]
    outFileName = params[[2]]
    simResults = params[[3]]
    underspecified = params[[4]]

    set.seed(seed)
    
    if (underspecified)
    {
      DataModel = buildDataModel(simResults$I_star_uspec, type = "overdispersion", 
                                 phi = 1)
      
      priorBetaIntercept = -3
      ExposureModel = buildExposureModel(simResults$X_underspec, nTpt = simResults$nTpt,
                                           nLoc = ncol(simResults$I_star_uspec),
                                           beta = c(priorBetaIntercept, 
                                                    rep(0, (ncol(simResults$X_underspec) - 1))),
                                           betaPriorPrecision = 0.1)
      I0 = I0 = rep(2, ncol(simResults$I_star))
      E0 = I0
      S0 = simResults$N_uspec[1,] - I0 - E0
      InitContainer = buildInitialValueContainer(simResults$I_star_uspec, simResults$N_uspec, 
                                                 S0 = S0, I0 = I0, E0 = E0, 
                                                 reinfection=FALSE, 
                                                 dataType = "I_star")
      DistanceModel = buildDistanceModel(simResults$distMatList_uspec,
                                         priorAlpha = 1, priorBeta = 20)
    }
    else
    {
        DataModel = buildDataModel(simResults$I_star, type = "overdispersion", 
                                   phi = 1)
        priorBetaIntercept = -3
        ExposureModel = buildExposureModel(simResults$X, nTpt = simResults$nTpt,
                                           nLoc = simResults$nLoc,
                                           beta = c(priorBetaIntercept, 
                                     rep(0,((length(simResults$beta_SE))-1))), 
                                     betaPriorPrecision = 0.1)
        I0 = rep(2, ncol(simResults$I_star))
        E0 = I0
        S0 = simResults$N[1,] - I0 - E0
        InitContainer = buildInitialValueContainer(simResults$I_star, simResults$N, 
                                                   S0 = S0, I0 = I0, E0 = E0, 
                                                   reinfection=FALSE, 
                                                   dataType = "I_star")
        DistanceModel = buildDistanceModel(simResults$distMatList,
                                           priorAlpha = 1, priorBeta = 20)
    }

    ReinfectionModel = buildReinfectionModel("SEIR")
    SamplingControl = buildSamplingControl(iterationStride=1000,
                                           sliceWidths = c(1,  # S_star
                                                           1,  # E_star
                                                           1, # I_star
                                                           1, # S0
                                                           1, # I0
                                                           0.05, # beta
                                                           0.05, # betaPrs
                                                           0.01, # rho
                                                           0.01, # gamma_ei
                                                           0.01, # gamma_ir
                                                           0.01 # phi
                                                          ))
    
    TransitionPriors = buildTransitionPriorsFromProbabilities(1-exp(-simResults$gamma_ei), 
                                                              1-exp(-simResults$gamma_ir), 
                                                              1000,
                                                              1000) 

    
    
    res = buildSEIRModel(outFileName,DataModel,ExposureModel,ReinfectionModel,
                         DistanceModel,TransitionPriors,
                         InitContainer,SamplingControl)

    res$setRandomSeed(seed)
    for (i in 1:ncol(res$I_star))
    {
      res$setTrace(i-1)
    }
   
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
    res$performHybridStep = 111

    # Store the model object in the global namespace of the node 
    # - can't pass these between sessions
    localModelObject <<- res
    return(list("model"=res,
                "fileName"=outFileName))    
}


runSpatialMisspecificationSimulation = function(cellIterations = 50, 
                                           genSeed=123123, 
                                           fitSeeds=c(812123,12301,5923),
                                           scratchLoc="/localscratch/Users/gown/")
{                     
    itrSeeds = genSeed + seq(1, cellIterations)
    main.cluster = makeCluster(2,outfile = paste("err", genSeed, ".txt", sep = ""))

    parLapply(main.cluster, 1:length(main.cluster), function(x){library(EARNMC)})

    mainLoop = function(genSeed)
    {
        #cl = makeCluster(3, outfile = paste("innerErr", genSeed, ".txt", sep = ""))
        cl = makeCluster(3)
        print("Cluster Created")
        clusterExport(cl, c("buildSpatialMisspecificationSimulationKernel"))
        print("Variables Exported.") 
     
        for (underspec in c(FALSE, TRUE))
        {
            f = function(genSeed)
            {
                spatialMisspecificationSimulationKernel(cl, genSeed, 
                                                   fitSeeds + genSeed,
                                                   underspec,
                                                   scratchLoc)
            }
            fname = paste("./simSp_results_Uspec_", underspec, "_", genSeed, 
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
        stopCluster(cl)
    }

    parLapplyLB(main.cluster, itrSeeds, mainLoop)
    stopCluster(main.cluster)
    TRUE
}
