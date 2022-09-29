# Import MCMC function

jointModel_IncrCumInc_Misc = function(fitlme,fitCoxCause1,fitCoxCause2,prior,startingValues,nknotsCause1 = 3,nknotsCause2 = 3,
                                      scaleVarSurv = 1,iterMaxSurv_1 = 1,iterMaxSurv_2 = 1,DoubleSamplVar = "Robs",TrueStatusVar = "status",
                                      full.knotsCause1 = NULL,full.knotsCause2 = NULL,timeVar = "times",alpha1 = 1,alpha2 = 1,
                                      ndraw = 100,thin = 1,store.b = F,nburn = 0,useGauleg = T,Gauleg_points = 30,tdf = 10)
{
  #########################################################################
  ### The Longitudinal model should be specified through an lme object  ###
  ### The survival models, for each competing risk, should be specified ###
  ### by coxph objects ####################################################
  #########################################################################
  
  #################
  ### Functions ###
  #################
  
  # Function specifying the boundness constraints: For each individual the sum of the CIFs
  # should be less than 1 at the observed survival time
  hin = function(paramSurv)
  {
    paramCause1 = paramSurv[1:lpCause1]
    paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
    cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
    
    1 - cincCause1 - cincCause2
  }
  
  hin.jac = function(paramSurv)
  {
    paramCause1 = paramSurv[1:lpCause1]
    paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    a1s1b =  1 + alpha1*int1
    a2s2b =  1 + alpha2*int2 
    
    # Compute the score vector
    qq1Mat = xCause1CombLong*qq1
    qq2Mat = xCause2CombLong*qq2
    
    aa1 = rep( a1s1b^(-(alpha1+1)/alpha1), each = length(pps) )
    aa2 = rep( a2s2b^(-(alpha2+1)/alpha2), each = length(pps) )
    
    mat = cbind(qq1Mat*aa1,qq2Mat*aa2)
    
    for (i in 1:G)
    {
      jacMat[i,] = -.colSums(mat[ids[[i]],],length(pps),lpCause1 + lpCause2,F) 
    }
    jacMat
  }
  
  f = function(paramSurv,computeDeriv = F)
  {
    paramCause1 = paramSurv[1:lpCause1]
    paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
    cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = -((alpha1+1)/alpha1)*log(1 + alpha1*int1) + linpredCause1
    DercincCause2 = -((alpha2+1)/alpha2)*log(1 + alpha2*int2) + linpredCause2
    
    logl = sum(DercincCause1[indDelta1]) + sum(DercincCause2[indDelta2]) + sum(log(1-cincCause1-cincCause2)*(1-delta))
    return(-logl)
    
  }
  
  fDer = function(paramSurv)
  {
    paramCause1 = paramSurv[1:lpCause1]
    paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    a1s1b =  1 + alpha1*int1
    a2s2b =  1 + alpha2*int2 
    cincCause1 = 1 - a1s1b^(-1/alpha1)
    cincCause2 = 1 - a2s2b^(-1/alpha2)
    st = 1-cincCause1-cincCause2
    
    # Compute the score vector
    qq1Mat = xCause1CombLong*qq1
    qq2Mat = xCause2CombLong*qq2
    
    tt1 = rep( a1s1b^(-(alpha1+1)/alpha1)/st, each = length(pps) )
    tt2 = rep( a2s2b^(-(alpha2+1)/alpha2)/st, each = length(pps) )
    aa1 = rep( 1/a1s1b, each = length(pps) )
    aa2 = rep( 1/a2s2b, each = length(pps) )
    
    score1 = .colSums(deltaCause1*xCause1Comb,G,lpCause1,F) - (1+alpha1)*.colSums( deltaCause1Long*(qq1Mat*aa1),nLong,lpCause1,F) 
    score1 = score1 - .colSums( (1 - deltaLong)*(qq1Mat*tt1),nLong,lpCause1,F)
    
    score2 = .colSums(deltaCause2*xCause2Comb,G,lpCause2,F) - (1+alpha2)*.colSums( deltaCause2Long*(qq2Mat*aa2),nLong,lpCause2,F) 
    score2 = score2 - .colSums( (1 - deltaLong)*(qq2Mat*tt2),nLong,lpCause2,F)
    
    # Score vector
    gradient = c(score1,score2)
    
    return(-gradient)
  }
  
  fDer1 = function(paramCause1)
  {
    # paramCause1 = paramSurv[1:lpCause1]
    # paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    a1s1b =  1 + alpha1*int1
    a2s2b =  1 + alpha2*int2 
    cincCause1 = 1 - a1s1b^(-1/alpha1)
    cincCause2 = 1 - a2s2b^(-1/alpha2)
    st = 1-cincCause1-cincCause2
    
    # Compute the score vector
    qq1Mat = xCause1CombLong*qq1
    #qq2Mat = xCause2CombLong*qq2
    
    tt1 = rep( a1s1b^(-(alpha1+1)/alpha1)/st, each = length(pps) )
    #tt2 = rep( a2s2b^(-(alpha2+1)/alpha2)/st, each = length(pps) )
    aa1 = rep( 1/a1s1b, each = length(pps) )
    #aa2 = rep( 1/a2s2b, each = length(pps) )
    
    score1 = .colSums(deltaCause1*xCause1Comb,G,lpCause1,F) - (1+alpha1)*.colSums( deltaCause1Long*(qq1Mat*aa1),nLong,lpCause1,F) 
    score1 = score1 - .colSums( (1 - deltaLong)*(qq1Mat*tt1),nLong,lpCause1,F)
    
    # score2 = .colSums(deltaCause2*xCause2Comb,G,lpCause2,F) - (1+alpha2)*.colSums( deltaCause2Long*(qq2Mat*aa2),nLong,lpCause2,F) 
    # score2 = score2 - .colSums( (1 - deltaLong)*(qq2Mat*tt2),nLong,lpCause2,F)
    
    # Score vector
    #gradient = c(score1,score2)
    
    return(-score1)
  }
  
  fDer2 = function(paramCause2)
  {
    # paramCause1 = paramSurv[1:lpCause1]
    # paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    a1s1b =  1 + alpha1*int1
    a2s2b =  1 + alpha2*int2 
    cincCause1 = 1 - a1s1b^(-1/alpha1)
    cincCause2 = 1 - a2s2b^(-1/alpha2)
    st = 1-cincCause1-cincCause2
    
    # Compute the score vector
    #qq1Mat = xCause1CombLong*qq1
    qq2Mat = xCause2CombLong*qq2
    
    #tt1 = rep( a1s1b^(-(alpha1+1)/alpha1)/st, each = length(pps) )
    tt2 = rep( a2s2b^(-(alpha2+1)/alpha2)/st, each = length(pps) )
    #aa1 = rep( 1/a1s1b, each = length(pps) )
    aa2 = rep( 1/a2s2b, each = length(pps) )
    
    # score1 = .colSums(deltaCause1*xCause1Comb,G,lpCause1,F) - (1+alpha1)*.colSums( deltaCause1Long*(qq1Mat*aa1),nLong,lpCause1,F) 
    # score1 = score1 - .colSums( (1 - deltaLong)*(qq1Mat*tt1),nLong,lpCause1,F)
    
    score2 = .colSums(deltaCause2*xCause2Comb,G,lpCause2,F) - (1+alpha2)*.colSums( deltaCause2Long*(qq2Mat*aa2),nLong,lpCause2,F) 
    score2 = score2 - .colSums( (1 - deltaLong)*(qq2Mat*tt2),nLong,lpCause2,F)
    
    # Score vector
    #gradient = c(score1,score2)
    
    return(-score2)
  }
  
  fsurvEvalRdouble = function(paramSurv,computeDeriv = F)
  {
    paramCause1 = paramSurv[1:lpCause1]
    paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
    cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = -((alpha1+1)/alpha1)*log(1 + alpha1*int1) + linpredCause1
    DercincCause2 = -((alpha2+1)/alpha2)*log(1 + alpha2*int2) + linpredCause2
    
    out = list()
    out$DercincCause1 = DercincCause1
    out$DercincCause2 = DercincCause2
    return(out)
    
  }
  
  fsurvEval = function(paramSurv)
  {
    paramCause1 = paramSurv[1:lpCause1]
    paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
    cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = -((alpha1+1)/alpha1)*log(1 + alpha1*int1) + linpredCause1
    DercincCause2 = -((alpha2+1)/alpha2)*log(1 + alpha2*int2) + linpredCause2
    
    logl = sum(DercincCause1[indDelta1]) + sum(DercincCause2[indDelta2]) + sum(log(1-cincCause1-cincCause2)*(1-delta))
    return(logl)
    
  }
  
  fsurvEval1 = function(paramCause1)
  {
    # paramCause1 = paramSurv[1:lpCause1]
    # paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
    cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = -((alpha1+1)/alpha1)*log(1 + alpha1*int1) + linpredCause1
    DercincCause2 = -((alpha2+1)/alpha2)*log(1 + alpha2*int2) + linpredCause2
    
    logl = sum(DercincCause1[indDelta1]) + sum(DercincCause2[indDelta2]) + sum(log(1-cincCause1-cincCause2)*(1-delta))
    return(logl)
    
  }
  
  fsurvEval2 = function(paramCause2)
  {
    # paramCause1 = paramSurv[1:lpCause1]
    # paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
    cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = -((alpha1+1)/alpha1)*log(1 + alpha1*int1) + linpredCause1
    DercincCause2 = -((alpha2+1)/alpha2)*log(1 + alpha2*int2) + linpredCause2
    
    logl = sum(DercincCause1[indDelta1]) + sum(DercincCause2[indDelta2]) + sum(log(1-cincCause1-cincCause2)*(1-delta))
    return(logl)
    
  }
  
  fsurvDer = function(paramSurv,iter = 5)
  {
    paramCause1 = paramSurv[1:lpCause1]
    paramCause2 = paramSurv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    xbetasCause1 = c(xCause1Long%*%betasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = c(xCause1 %*% betasCause1 + splBasisCause1 %*% gammasCause1)
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    xbetasCause2 = c(xCause2Long%*%betasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = c(xCause2 %*% betasCause2 + splBasisCause2 %*% gammasCause2)
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    a1s1b =  1 + alpha1*int1
    a2s2b =  1 + alpha2*int2 
    cincCause1 = 1 - a1s1b^(-1/alpha1)
    cincCause2 = 1 - a2s2b^(-1/alpha2)
    st = 1-cincCause1-cincCause2
    
    # Compute the score vector
    qq1Mat = xCause1CombLong*qq1
    qq2Mat = xCause2CombLong*qq2
    
    for (j in 1:lpCause1)
    {
      vec = qq1Mat[,j]
      s1der[,j] = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(vec[x]))))
    }
    
    for (j in 1:lpCause2)
    {
      vec = qq2Mat[,j]
      s2der[,j] = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(vec[x]))))
    }
    
    gg1 = a1s1b^(-(alpha1+1)/alpha1)/st
    gg2 = a2s2b^(-(alpha2+1)/alpha2)/st
    
    tt1 = rep( gg1, each = length(pps) )
    tt2 = rep( gg2, each = length(pps) )
    aa1 = rep( 1/a1s1b, each = length(pps) )
    aa2 = rep( 1/a2s2b, each = length(pps) )
    
    score1 = .colSums( deltaCause1*(xCause1Comb-(1+alpha1)*s1der/a1s1b) - ((1-delta)*gg1)*s1der ,G,lpCause1,F) 
    score2 = .colSums( deltaCause2*(xCause2Comb-(1+alpha2)*s2der/a2s2b) - ((1-delta)*gg2)*s2der ,G,lpCause2,F)
    
    hess1 = (1+alpha1)*( crossprod(xCause1CombLong*(sqrt(qq1*aa1)*deltaCause1Long)) - alpha1*crossprod((s1der/a1s1b)*deltaCause1) )
    hess1 = hess1 + crossprod( s1der*(gg1*(1-delta)) ) - (alpha1+1)*crossprod(s1der*(sqrt(gg1/a1s1b)*(1-delta)) )
    hess1 = hess1 + crossprod( xCause1CombLong*(sqrt(qq1*tt1)*(1-deltaLong))  )
    
    hess2 = (1+alpha2)*( crossprod(xCause2CombLong*(sqrt(qq2*aa2)*deltaCause2Long)) - alpha2*crossprod((s2der/a2s2b)*deltaCause2) )
    hess2 = hess2 + crossprod( s2der*(gg2*(1-delta)) ) - (alpha2+1)*crossprod(s2der*(sqrt(gg2/a2s2b)*(1-delta)) )
    hess2 = hess2 + crossprod( xCause2CombLong*(sqrt(qq2*tt2)*(1-deltaLong))  )
    
    # lpcause1*lpcause2
    hess12 = crossprod(s1der*(gg1*(1-delta)),s2der*gg2)
    
    # Score vector
    gradient = c(score1,score2)
    
    hessmat[1:lpCause1,1:lpCause1] = hess1
    hessmat[-c(1:lpCause1),-c(1:lpCause1)] = hess2
    hessmat[1:lpCause1,-c(1:lpCause1)] = hess12
    hessmat[-c(1:lpCause1),1:lpCause1] = t(hess12)
    
    out = list()
    out$gradient = gradient
    out$Information = hessmat
    
    return(out)
  }
  
  fbetaEval = function(Lbeta)
  {
    # Log posterior from the Longitudinal model
    logl1 = - 0.5*sum( (LinfBetaLong %*% Lbeta)^2) + sum( Lbeta*right )
    
    # Predictions of the marker value at survival time
    mtpred = (xTime %*% Lbeta) + zTimeb
    
    # Predictions of the marker value at additional points (required to integrate the subdistribution hazard)
    mtpredLong = (xTimeLong %*% Lbeta) + zTimeLongb
    
    
    # Design matrices of the survival submodel
    if (lbetasCause1==1)
    {
      xCause1 = mtpred
      xCause1Long = mtpredLong
    } else{
      xCause1[,lbetasCause1] = mtpred
      xCause1Long[,lbetasCause1] = mtpredLong
    }
    
    if (lbetasCause2==1)
    {
      xCause2 = mtpred
      xCause2Long = mtpredLong
    } else{
      xCause2[,lbetasCause1] = mtpred
      xCause2Long[,lbetasCause1] = mtpredLong
    }
    
    xbetasCause1 = mtpredLong*assoc1
    linpredCause1Long = xbetasCause1 + uCause1Long
    explinpredCause1Long = exp(linpredCause1Long)
    linpredCause1 = mtpred %*% assoc1
    
    xbetasCause2 = mtpredLong*assoc2
    linpredCause2Long = xbetasCause2 + uCause2Long
    explinpredCause2Long = exp(linpredCause2Long)
    linpredCause2 = mtpred %*% assoc2
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq1 = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq1[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq2 = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq2[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
    cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = -((alpha1+1)/alpha1)*log(1 + alpha1*int1) + linpredCause1
    DercincCause2 = -((alpha2+1)/alpha2)*log(1 + alpha2*int2) + linpredCause2
    
    # Overall survival 
    st = 1 - cincCause1 - cincCause2
    
    logl2 = sum(DercincCause1[indDelta1]) + sum(DercincCause2[indDelta2]) + sum(log(st)*(1-delta)) 
    
    return(logl1+logl2)
    
  }
  
  # Checking arguments
  if ( !( class(fitlme) == "lme" & "coxph" %in% class(fitCoxCause1) & "coxph" %in% class(fitCoxCause2) ) )
  {
    print("Not a lme of coxph objects")
    break
  }
  
  if ( (!is.null(nknotsCause1) &  !is.null(full.knotsCause1)) |  (!is.null(nknotsCause2) &  !is.null(full.knotsCause2)) )
  {
    print("Number of knots and full knots both missing")
    break
  }
  
  if ( ncol(fitCoxCause1$y) != 2 | ncol(fitCoxCause2$y) != 2  )
  {
    print("Entry times not allowed")
    break
  }
  
  if ( sum(fitCoxCause1$y[,1] - fitCoxCause2$y[,1]) != 0 )
  {
    print("Different survival times in the coxph models")
    break
  }
  
  #############################
  # Prepare longitudinal data #
  #############################
  
  # Response variable
  y = as.vector(getResponse(fitlme))
  n = length(y)
  
  # Get the data from the LME object
  dataLme = fitlme$data
  
  # Design matrix of the fixed effects
  Lx = model.matrix(formula(fitlme),dataLme)
  Lxtx = crossprod(Lx)
  
  # Dimension of the fixed effects
  Lp = ncol(Lx)
  
  Lx = as.matrix(Lx[1:n,1:Lp])
  colnames(Lx) = NULL
  rownames(Lx) = NULL
  
  # Design matrix of the random effects
  Lz = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLme)
  
  # Number of random effects
  q = ncol(Lz)
  
  Lz = as.matrix(Lz[1:n,1:q])
  colnames(Lz) = NULL
  rownames(Lz) = NULL
  Lzmatrix = Lz
  
  # Number of observations by subject
  id = as.numeric(as.character(fitlme$groups[,1]))
  ni = as.numeric(table(id))
  
  # Number of subjects
  G = length(unique(id))
  
  # z: design matrix for random effects
  # Split z matrix by id
  Lz = data.frame(Lz)
  names(Lz) = NULL
  
  # Splitting the fixed and the random effects matrices
  zlist = lapply(split(Lz,id),function(x) {row.names(x)=NULL;x=as.matrix(x)})
  xlist = lapply(split(data.frame(Lx),id),function(x) {row.names(x)=NULL;x=as.matrix(x)})
  zzlist = lapply(zlist,FUN=function(x) crossprod(x) ) # cross products
  ylist = split(y,id)
  
  # Block diagonal z matrix dim = n*(q*G)
  ztbdiag = t(bdiag(zlist))
  
  # Prior paramters
  Lc0 = prior$Lc0
  Lmu0 = prior$Lmu0
  Smu0 = prior$Smu0
  Sc0 = prior$Sc0
  df = prior$df
  A = prior$A
  lambda1 = prior$lambda1
  lambba2 = prior$lambba2
  a1 = prior$a1
  a2 = prior$a2
  b1 = prior$b1
  b2 = prior$b2
  
  # Some manipulations
  Linvc0 = solve(Lc0)
  Linvc0mu0 = Linvc0 %*% Lmu0
  Sinvc0 = solve(Sc0)
  
  #####################################################
  ### Initial values for the longitudinal sub-model ###
  #####################################################
  Lbeta = fixed.effects(fitlme)
  names.Lbeta = names(Lbeta)
  names(Lbeta) = NULL
  Lbeta_lmm = Lbeta
  
  # Within subject precision
  omega = startingValues$omega
  
  # Variance matrix of random effects
  D = startingValues$D
  rownames(D) = NULL
  colnames(D) = NULL
  Dinv = solve(D)
  
  # Prediction of the random effects
  b = as.matrix(ranef(fitlme))
  colnames(b) = NULL
  rownames(b) = NULL
  b_lmm = b
  
  # For cholesky decomposition
  pivot = FALSE
  tol = -1
  
  #############################
  ### Prepare survival data ###
  #############################
  
  # Get the survival times and the failure indicators for cause 1
  if (ncol(fitCoxCause1$y) == 2)
  {
    times = fitCoxCause1$y[,1]
    deltaCause1 = fitCoxCause1$y[,2]
    
    names(times)=NULL;names(deltaCause1)=NULL
  }
  
  # Get the failure indicators for cause 1
  if (ncol(fitCoxCause2$y) == 2)
  {
    #times = fitCoxCause1$y[,1]
    deltaCause2 = fitCoxCause2$y[,2]
    
    names(times)=NULL;names(deltaCause2)=NULL
  }
  
  # Failure indicator
  # delta = deltaCause1 + deltaCause2
  # indcens = which(delta == 0)
  # indDelta1 = which(deltaCause1==1)
  # indDelta2 = which(deltaCause2==1)
  
  # Failure indicators based on the observed failure causes
  deltaCause1Obs = deltaCause1
  deltaCause2Obs = deltaCause2
  
  # data frame based for the cox models
  dataid = get(as.character(fitCoxCause1$call$data),envir = .GlobalEnv)
  
  # Indicator of doubly sampled individuals
  Rdouble = dataid[,DoubleSamplVar]
  statusTrue = dataid[,TrueStatusVar]
  deltaCause1True = 1*(statusTrue == 1)
  deltaCause2True = 1*(statusTrue == 2)
  Rind = which(Rdouble == 1)
  Rmis = which(Rdouble == 0)
  nRmis = length(Rmis)
  
  deltaCause1[Rind] <- deltaCause1True[Rind]
  deltaCause2[Rind] <- deltaCause2True[Rind]
  
  deltaCause1Start <- deltaCause1
  deltaCause2Start <- deltaCause2
  
  if ( !is.null(nknotsCause1) )
  {
    # nknots is the number of the internal knots
    # knots: internal knots
    centiles = seq(0,1,length.out = nknotsCause1 + 2)[-c(1,nknotsCause1 + 2)]
    knotsCause1 = quantile(times[deltaCause1==1],probs = centiles)
    mCause1  = length(knotsCause1)
    
    # Boundary knots
    bknotsCause1 = range(times)
  } else{
    
    knotsCause1 = full.knotsCause1[-c(1,length(full.knotsCause1))]
    bknotsCause1 = full.knotsCause1[c(1,length(full.knotsCause1))]
    mCause1  = length(knotsCause1)
  }
  
  if ( !is.null(nknotsCause2) )
  {
    # nknots is the number of the internal knots
    # knots: internal knots
    centiles = seq(0,1,length.out = nknotsCause2 + 2)[-c(1,nknotsCause2 + 2)]
    knotsCause2 = quantile(times[deltaCause2==1],probs = centiles)
    mCause2  = length(knotsCause2)
    
    # Boundary knots
    bknotsCause2 = range(times)
  } else{
    
    knotsCause2 = full.knotsCause2[-c(1,length(full.knotsCause2))]
    bknotsCause2 = full.knotsCause2[c(1,length(full.knotsCause2))]
    mCause2  = length(knotsCause2)
  }
  
  # Get the predicted values of the marker at the survival time for each subject
  dataLmeid = dataLme[!duplicated(id),]
  dataLmeid[,timeVar] = times
  predsLme = predict(fitlme,dataLmeid)
  
  xTime = model.matrix(formula(fitlme),data = dataLmeid)
  colnames(xTime) = NULL
  rownames(xTime) = NULL
  
  zTime = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeid)
  colnames(zTime) = NULL
  rownames(zTime) = NULL
  
  # Baseline covariate for the survival submodels
  xCause1 = fitCoxCause1$x
  xCause2 = fitCoxCause2$x
  
  fitCox1_temp = coxph(Surv(times,deltaCause1) ~ cbind(xCause1,predsLme),x = T)
  fitCox2_temp = coxph(Surv(times,deltaCause2) ~ cbind(xCause2,predsLme),x = T)
  
  # Initial values from the fixed and random effects
  for (i in 1:G)
  {
    Linf = Dinv + zzlist[[i]]*omega
    b[i,] = rmvnorm(1,mean = b[i,],sigma = solve(Linf))
  }
  Lbeta = startingValues$Lbeta
  
  # Level 0 residuals
  Lres = y - c(Lx%*%Lbeta)
  
  # z_{i}^{T}(y_{i}-x_{i}\beta)*omega
  ztres = as.vector(ztbdiag %*% Lres)*omega
  ztres = .Internal( matrix(ztres, nrow = G, ncol = q, byrow = TRUE, dimnames = NULL,FALSE, FALSE ) )
  
  # Starting values for the survival model
  fitStart = survModel_IncrCumInc_addConstr(fitCox1_temp,fitCox2_temp,fitlme,Lbeta = Lbeta,reffects = b,timeVar = timeVar,
                                            nknotsCause1 = nknotsCause1,nknotsCause2 = nknotsCause2,alpha1 = alpha1,alpha2 = alpha2,
                                            useGauleg = useGauleg,Gauleg_points = Gauleg_points)
  paramSurv = fitStart$param
  names(paramSurv) = NULL
  
  # Number of parameters for cause 1
  lpCause1 = fitStart$lpcause1
  lbetasCause1 = length(coef(fitCox1_temp))
  lgammasCause1 = lpCause1 - lbetasCause1
  
  # Number of parameters for cause 2
  lpCause2 = fitStart$lpcause2
  lbetasCause2 = length(coef(fitCox2_temp))
  lgammasCause2 = lpCause2 - lbetasCause2
  
  Smu0_1 = Smu0[1:lpCause1]
  Smu0_2 = Smu0[-c(1:lpCause1)]
  Sc0_1 = Sc0[1:lpCause1,1:lpCause1] 
  Sc0_2 = Sc0[-c(1:lpCause1),-c(1:lpCause1)]
  Sinvc0_1 = solve(Sc0_1)
  Sinvc0_2 = solve(Sc0_2)
  
  paramCause1 = paramSurv[1:lpCause1]
  paramCause2 = paramSurv[-c(1:lpCause1)]
  
  # Parameters for cause 1
  betasCause1 = paramCause1[1:lbetasCause1]
  gammasCause1 = paramCause1[-c(1:lbetasCause1)]
  
  # Parameters for cause 2
  betasCause2 = paramCause2[1:lbetasCause2]
  gammasCause2 = paramCause2[-c(1:lbetasCause2)]
  #rm(paramCause1,paramCause2)
  assoc1 = betasCause1[lbetasCause1]
  assoc2 = betasCause2[lbetasCause2]
  
  # 7-15 Gauss-Kronrod quadradure points and weights
  gaussKronrod <-
    function (k = 15) {
      sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
              0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
              -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
              0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
      wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
      wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975, 
               0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
      if (k == 7) 
        list(sk = sk[1:7], wk = wk7)
      else
        list(sk = sk, wk = wk15)
    }
  
  gauleg <- function(n,a=-1,b=1) {# Gauss-Legendre: returns x,w so that
    #\int_a^b f(x) dx \doteq \sum w_i f(x_i)
    EPS <- 3.e-14
    m <- trunc((n+1)/2)
    xm <- 0.5*(b+a)
    xl <- 0.5*(b-a)
    x <- w <- rep(-1,n)
    for (i in 1:m) {
      z <- cos(pi*(i-.25)/(n+.5))
      z1 <- z+1
      while (abs(z-z1) > EPS) {
        p1 <- 1
        p2 <- 0
        for (j in 1:n) {# recursively evaluates pn(x)
          p3 <- p2
          p2 <- p1
          p1 <- ((2*j-1)*z*p2-(j-1)*p3)/j
        }
        pp <- n*(z*p1-p2)/(z*z-1)
        z1 <- z
        z <- z1-p1/pp #Newton iteration
      }
      x[i] <- xm-xl*z
      x[n+1-i] <- xm+xl*z
      w[i] <- 2*xl/((1-z*z)*pp*pp)
      w[n+1-i] <- w[i]
    }
    list(x=x,w=w)
  }
  
  if (useGauleg == F)
  {
    gk = gaussKronrod()
    pps = gk$sk
    wws = gk$wk
  } else {
    gk = gauleg(Gauleg_points)
    pps = gk$x
    wws = gk$w
  }
  
  # Transformation of time needed for gaussKronrod integration
  t1 = 0.5*times
  t1Long = rep(t1,each = length(pps))
  ppsLong = rep(pps,G)
  wwsLong = rep(wws,G)
  ind = rep(1:G,each = length(pps))
  
  # Additional time points used to integrate the subditribution hazard
  timesForInt = t1Long*ppsLong + t1Long
  
  # Construct the spline basis for cause 1
  splBasisCause1Long = bSpline(timesForInt,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T)
  splBasisCause1 = bSpline(times,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T)
  
  # Construct the spline basis for cause 2
  splBasisCause2Long = bSpline(timesForInt,knots = knotsCause2,Boundary.knots = bknotsCause2,intercept = T)
  splBasisCause2 = bSpline(times,knots = knotsCause2,Boundary.knots = bknotsCause2,intercept = T)
  
  #uCause1Long = c(splBasisCause1Long %*% gammasCause1)
  #uCause2Long = c(splBasisCause2Long %*% gammasCause2)
  
  # Make predictions for the marker values at additional time points (used to integrate the subdistribution hazard) 
  dataLmeidLong = dataLmeid[ind,]
  dataLmeidLong[,timeVar] = timesForInt
  
  xTimeLong = model.matrix(formula(fitlme),data = dataLmeidLong)
  rownames(xTimeLong) = NULL
  colnames(xTimeLong) = NULL
  
  zTimeLong = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidLong)
  rownames(zTimeLong) = NULL
  colnames(zTimeLong) = NULL
  
  # Indices of id's
  ids = list()
  for (i in 1:G){ids[[i]] = which(ind==i)}
  
  # Based on the current value of the fixed and the random effects,
  # we create the appropriate design matrices for the survival models
  # Design matrices of the baseline covariates of the survival model
  xCause1 = fitCox1_temp$x
  colnames(xCause1) = NULL
  xCause1Long = as.matrix(xCause1[ind,])
  
  xCause2 = fitCox2_temp$x
  colnames(xCause2) = NULL
  xCause2Long = as.matrix(xCause2[ind,])
  
  timesLbeta = xTime %*% Lbeta
  timesLbetaLong = xTimeLong %*% Lbeta
  
  # Predictions of the marker value at the survival time of each subject
  mtpred = timesLbeta + rowSums(zTime*b)
  
  # Predictions of the marker value at additional points (required to integrate the subdistribution hazard)
  mtpredLong = timesLbetaLong + rowSums(zTimeLong*b[ind,])
  
  ####################################################################################
  ### Define starting point in the BFGS algorithm used for the survival parameters ###
  ####################################################################################
  if (lbetasCause1==1)
  {
    xCause1 = xTime %*% Lbeta_lmm + rowSums(zTime*b_lmm)
    xCause1Long = xTimeLong %*% Lbeta_lmm + rowSums(zTimeLong*b_lmm[ind,])
  } else{
    xCause1[,lbetasCause1] = xTime %*% Lbeta_lmm + rowSums(zTime*b_lmm)
    xCause1Long[,lbetasCause1] = xTimeLong %*% Lbeta_lmm + rowSums(zTimeLong*b_lmm[ind,])
  }
  
  if (lbetasCause2==1)
  {
    xCause2 = xTime %*% Lbeta_lmm + rowSums(zTime*b_lmm)
    xCause2Long = xTimeLong %*% Lbeta_lmm + rowSums(zTimeLong*b_lmm[ind,])
  } else{
    xCause2[,lbetasCause1] = xTime %*% Lbeta_lmm + rowSums(zTime*b_lmm)
    xCause2Long[,lbetasCause1] = xTimeLong %*% Lbeta_lmm + rowSums(zTimeLong*b_lmm[ind,])
  }
  
  jacMat = matrix(nr = G, nc = lpCause1 + lpCause2)
  s1der = matrix(NA,nr = G,nc = lpCause1)
  s2der = matrix(NA,nr = G,nc = lpCause2)
  hessmat = matrix(NA,nr = lpCause1 + lpCause2,nc = lpCause1 + lpCause2)
  
  pi11 = rbeta(1,shape1 = a1 + sum(deltaCause1*deltaCause1Obs),shape2 = b1 + sum(deltaCause1*deltaCause2Obs) )
  pi22 = rbeta(1,shape1 = a2 + sum(deltaCause2*deltaCause2Obs),shape2 = b2 + sum(deltaCause2*deltaCause1Obs) )
  
  # To get the starting point for the BFGS algorithm (not to be used as a starting point of the MCMC)
  for (j in 1:5)
  {
    fitCox1_temp = coxph(Surv(times,deltaCause1) ~ cbind(xCause1[,-ncol(xCause1)],predsLme),x = T)
    fitCox2_temp = coxph(Surv(times,deltaCause2) ~ cbind(xCause2[,-ncol(xCause2)],predsLme),x = T)
    
    # Calculate MLE of the survival model parameters based on the estimates from the LMM
    fitStart_approx = survModel_IncrCumInc_addConstr(fitCox1_temp,fitCox2_temp,fitlme,Lbeta = Lbeta_lmm,reffects = b_lmm,timeVar = timeVar,
                                                     nknotsCause1 = NULL,nknotsCause2 = NULL,full.knotsCause1 = c(bknotsCause1[1],knotsCause1,bknotsCause1[2]),
                                                     full.knotsCause2 = c(bknotsCause2[1],knotsCause2,bknotsCause2[2]),alpha1 = alpha1,alpha2 = alpha2,
                                                     useGauleg = useGauleg,Gauleg_points = Gauleg_points)
    
    # Approximated information matrix for the survival model parameters
    hessStart = fitStart_approx$hess + Sinvc0
    hessStart = (hessStart + t(hessStart))/2
    hessStart1 = hessStart[1:lpCause1,1:lpCause1]
    hessStart2 = hessStart[-c(1:lpCause1),-c(1:lpCause1)]
    
    Sigmastart = solve(hessStart/scaleVarSurv)
    Sigmastart1 = solve(hessStart[1:lpCause1,1:lpCause1]/scaleVarSurv)
    Sigmastart2 = solve(hessStart[-c(1:lpCause1),-c(1:lpCause1)]/scaleVarSurv)
    
    # Lstart = chol(hessStart/scaleVarSurv)
    # Sigmastart = solve(hessStart/scaleVarSurv)
    se_start = sqrt(diag(solve(hessStart)))
    
    # Starting value of the BFGS update for the survival 
    paramSurvStart = fitStart_approx$param
    names(paramSurvStart) = NULL
    
    # For the non-doubly sampled individuals, simulate the ``true'' failure cause conditional on all the observed information and 
    # the current parameter values
    temp = fsurvEvalRdouble(paramSurvStart)
    logf12 = temp$DercincCause1 - temp$DercincCause2
    
    for (i in 1:100)
    {
      miscProb = log(pi11/(1 - pi22))*deltaCause1Obs + log((1 - pi11)/pi22)*deltaCause2Obs
      miscProb = miscProb[Rmis]
      
      prob1 = exp(logf12[Rmis] + miscProb)
      prob1 = prob1/(1+prob1)
      
      deltaCause1[Rmis] = rbinom(nRmis,size = 1,prob = prob1)
      deltaCause2[Rmis] = 1 - deltaCause1[Rmis]
      
      deltaCause1Long = rep( deltaCause1, each = length(pps) )
      deltaCause2Long = rep( deltaCause2, each = length(pps) )
      deltaLong = deltaCause1Long + deltaCause2Long
      
      # Failure indicators
      delta = deltaCause1 + deltaCause2
      indcens = which(delta == 0)
      indDelta1 = which(deltaCause1==1)
      indDelta2 = which(deltaCause2==1)
      
      # Update the misclassification parameters
      # pi11 = rbeta(1,shape1 = a1 + sum(deltaCause1*deltaCause1Obs),shape2 = b1 + sum(deltaCause1*deltaCause2Obs) )
      # pi22 = rbeta(1,shape1 = a2 + sum(deltaCause2*deltaCause2Obs),shape2 = b2 + sum(deltaCause2*deltaCause1Obs) )
      pi11 <- mean(deltaCause1Obs[deltaCause1==1])
      pi22 <- mean(deltaCause2Obs[deltaCause2==1])
      #print(c(pi11,pi22))
    }
    
  }
  
  # We get the actual initial values
  deltaCause1 = deltaCause1Start
  deltaCause2 = deltaCause2Start
  
  deltaCause1Long = rep( deltaCause1, each = length(pps) )
  deltaCause2Long = rep( deltaCause2, each = length(pps) )
  deltaLong = deltaCause1Long + deltaCause2Long
  
  # Failure indicators
  delta = deltaCause1 + deltaCause2
  indcens = which(delta == 0)
  indDelta1 = which(deltaCause1==1)
  indDelta2 = which(deltaCause2==1)
  
  # Design matrices of the survival submodel
  if (lbetasCause1==1)
  {
    xCause1 = mtpred
    xCause1Long = mtpredLong
  } else{
    xCause1[,lbetasCause1] = mtpred
    xCause1Long[,lbetasCause1] = mtpredLong
  }
  
  if (lbetasCause2==1)
  {
    xCause2 = mtpred
    xCause2Long = mtpredLong
  } else{
    xCause2[,lbetasCause1] = mtpred
    xCause2Long[,lbetasCause1] = mtpredLong
  }
  
  # Combined design matrices of the survival submodel (including the B-splines basis)
  xCause1Comb = cbind(xCause1,splBasisCause1)
  xCause1CombLong = cbind(xCause1Long,splBasisCause1Long)
  
  xCause2Comb = cbind(xCause2,splBasisCause2)
  xCause2CombLong = cbind(xCause2Long,splBasisCause2Long)
  
  if (lbetasCause1==1)
  {
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
  } else{
    uCause1Long = c(splBasisCause1Long %*% gammasCause1 + as.matrix(xCause1Long[,-lbetasCause1]) %*% betasCause1[-lbetasCause1] )
  }
  
  if (lbetasCause2==1)
  {
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
  } else{
    uCause2Long = c(splBasisCause2Long %*% gammasCause2 + as.matrix(xCause2Long[,-lbetasCause2]) %*% betasCause2[-lbetasCause2] )
  }
  
  # Failure indicators for the `extended` design matrices of the survival submodels
  nLong = G*length(pps)
  deltaCause1Long = rep( deltaCause1, each = length(pps) )
  deltaCause2Long = rep( deltaCause2, each = length(pps) )
  deltaLong = deltaCause1Long + deltaCause2Long
  
  # Starting values for the misclassification parameters
  #pi11 = sum(deltaCause1*deltaCause1Obs)/sum(deltaCause1)
  #pi22 = sum(deltaCause2*deltaCause2Obs)/sum(deltaCause2)
  pi11 = (a1 + sum(deltaCause1*deltaCause1Obs))/(a1 + sum(deltaCause1*deltaCause1Obs)+b1 + sum(deltaCause1*deltaCause2Obs))
  pi22 = (a2 + sum(deltaCause2*deltaCause2Obs))/(a2 + sum(deltaCause2*deltaCause2Obs)+b2 + sum(deltaCause2*deltaCause1Obs))
  
  
  # For the non-doubly sampled individuals, simulate the ``true'' failure cause conditional on all the observed information and 
  # the current parameter values
  temp = fsurvEvalRdouble(paramSurv)
  logf12 = temp$DercincCause1 - temp$DercincCause2
  
  miscProb = log(pi11/(1 - pi22))*deltaCause1Obs + log((1 - pi11)/pi22)*deltaCause2Obs
  miscProb = miscProb[Rmis]
  
  prob1 = exp(logf12[Rmis] + miscProb)
  prob1 = prob1/(1+prob1)
  
  deltaCause1[Rmis] = rbinom(nRmis,size = 1,prob = prob1)
  deltaCause2[Rmis] = 1 - deltaCause1[Rmis]
  
  deltaCause1Long = rep( deltaCause1, each = length(pps) )
  deltaCause2Long = rep( deltaCause2, each = length(pps) )
  deltaLong = deltaCause1Long + deltaCause2Long
  
  # Failure indicators
  delta = deltaCause1 + deltaCause2
  indcens = which(delta == 0)
  indDelta1 = which(deltaCause1==1)
  indDelta2 = which(deltaCause2==1)
  
  # Update the misclassification parameters
  pi11 = rbeta(1,shape1 = a1 + sum(deltaCause1*deltaCause1Obs),shape2 = b1 + sum(deltaCause1*deltaCause2Obs) )
  pi22 = rbeta(1,shape1 = a2 + sum(deltaCause2*deltaCause2Obs),shape2 = b2 + sum(deltaCause2*deltaCause1Obs) )
  
  Sdraws = matrix(nrow = ndraw*thin,ncol = lpCause1+lpCause2)
  Ldraws = matrix(nrow = ndraw*thin,ncol = Lp + 1)
  Mdraws = matrix(nrow = ndraw*thin,ncol = 2)
  drawsD = array(dim = c(q,q,ndraw*thin))
  accept_re = rep(0,G)
  accept_surv1 = 0
  accept_surv2 = 0
  accept_Lbeta = 0
  
  if (store.b == T)
  {
    draws.b = array(dim = c(G,q,ndraw*thin))
  }
  
  # Iteration index
  it = - nburn
  jj = 0
  
  # MCMC begins
  while(it<ndraw*thin)
  {
    it = it + 1
    
    # Update the values of the random effects
    for (i in 1:G)
    {
      # Matrices/vectors for the ith subject
      zz = zzlist[[i]]
      zi = zlist[[i]]
      ziLong = zTimeLong[ids[[i]],];if (q==1){ziLong = as.matrix(ziLong)}
      LxibetasLong = timesLbetaLong[ids[[i]]]
      ui1Long = uCause1Long[ids[[i]]]
      ui2Long = uCause2Long[ids[[i]]]
      Linf = Dinv + zz*omega
      bi = b[i,]
      
      # Cholesky factorization of Linf (upper triangular)
      Llong = .Internal(La_chol(Linf, pivot, tol))
      
      # The posterior mode using the longitudinal model only
      bb = ztres[i,]
      mu1 = .Internal( backsolve(Llong, .Internal(backsolve(Llong,bb,q,TRUE,TRUE)) ,q,TRUE,FALSE )  )
      
      
      for (j in 1:1)
      {
        # Score vector (gradient at starting value)
        score  = - Linf %*% mu1 + bb
        
        if (deltaCause1[i]==1)
        {
          # \int_{0}^{Ti}\lambda_{01}(s)\exp\{ a1*m_{i}(t) \}
          ints1 = exp( ui1Long + ( LxibetasLong + c(ziLong %*% mu1) )*assoc1 )
          s1b = 0.5*times[i]*sum( wws*ints1 )
          a1s1b = 1 + alpha1*s1b
          
          s1derb = 0.5*times[i]*assoc1*colSums( ziLong*(ints1*wws) )
          s1derbbt = (0.5*times[i]*assoc1^2)*crossprod(ziLong*sqrt(ints1*wws))
          
          # The score vector of f(bi|others)
          score = score + assoc1*zTime[i,] - (1+alpha1)*s1derb/a1s1b
          
          # Information matrix
          inf = Linf + (1+alpha1)*(s1derbbt*a1s1b - alpha1*tcrossprod(s1derb))/a1s1b^2
          
        }
        
        if (deltaCause2[i]==1)
        {
          # \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \} with respect to b_{i}
          ints2 = exp( ui2Long + ( LxibetasLong + c(ziLong %*% mu1) )*assoc2 )
          s2b = 0.5*times[i]*sum( wws*ints2 )
          a2s2b = 1 + alpha2*s2b
          
          s2derb = 0.5*times[i]*assoc2*colSums( ziLong*(ints2*wws) )
          s2derbbt = (0.5*times[i]*assoc2^2)*crossprod(ziLong*sqrt(ints2*wws))
          
          # The score vector of f(bi|others)
          score = score + assoc2*zTime[i,] - (1+alpha2)*s2derb/a2s2b
          
          # Information matrix
          inf = Linf + (1+alpha2)*(s2derbbt*a2s2b - alpha2*tcrossprod(s2derb))/a2s2b^2
        }
        
        if (delta[i]==0)
        {
          # Gradient of \int_{0}^{Ti}\lambda_{01}(s)\exp\{ a1*m_{i}(t) \} with respect to b_{i}
          ints1 = exp( ui1Long + ( LxibetasLong + c(ziLong %*% mu1) )*assoc1 )
          s1derb = 0.5*times[i]*assoc1*colSums( ziLong*(ints1*wws) )
          s1derbbt = (0.5*times[i]*assoc1^2)*crossprod(ziLong*sqrt(ints1*wws))
          
          # Gradient of \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \} with respect to b_{i}
          ints2 = exp( ui2Long + ( LxibetasLong + c(ziLong %*% mu1) )*assoc2 )
          s2derb = 0.5*times[i]*assoc2*colSums( ziLong*(ints2*wws) )
          s2derbbt = (0.5*times[i]*assoc2^2)*crossprod(ziLong*sqrt(ints2*wws))
          
          # \int_{0}^{Ti}\lambda_{01}(s)\exp\{ a1*m_{i}(t) \}
          s1b = 0.5*times[i]*sum( wws*ints1 )
          a1s1b = 1 + alpha1*s1b
          
          # \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \} with respect to b_{i}
          s2b = 0.5*times[i]*sum( wws*ints2 )
          a2s2b = 1 + alpha2*s2b

          # The score vector of f(bi|others)
          st = a1s1b^(-1/alpha1) + a2s2b^(-1/alpha2) - 1
          scoreSurv = ( (a1s1b^(-(1+alpha1)/alpha1))*s1derb + (a2s2b^(-(1+alpha2)/alpha2))*s2derb )/st
          score = score - scoreSurv
          
          inf = Linf + ( -(1+alpha1)*a1s1b^(-(1+alpha1)/alpha1-1)*tcrossprod(s1derb) + a1s1b^(-(1+alpha1)/alpha1)*s1derbbt )/st
          inf = inf +  ( -(1+alpha2)*a2s2b^(-(1+alpha2)/alpha2-1)*tcrossprod(s2derb) + a2s2b^(-(1+alpha2)/alpha2)*s2derbbt )/st
          inf = inf + tcrossprod(scoreSurv)
          
          # If the CIFs are not bounded by 1 at mu1, set the information matrix equal to the information matrix 
          # from the longitudinal submodel
          if ( st < 0 )
          {
            inf = Linf
          }
          
        }
        
        # One NR iteration
        L = .Internal(La_chol(inf, pivot, tol))
        mu1 = .Internal( backsolve(L, .Internal(backsolve(L,score,q,TRUE,TRUE)) ,q,TRUE,FALSE )  ) + mu1
      }
      
      Linfmu1 = Linf %*% mu1
      # Propose one value from a multivariate normal with mean mu1 and covariance matrix
      # (D^{-1}+omega Z_{i}^T Z_{i})^{-1}
      epsilon = rnorm(q)
      bican = .Internal(backsolve(Llong,epsilon,q,TRUE,FALSE)) + mu1
      
      bdiff = bican - bi
      
      logp = sum(bdiff*bb)
      
      if (deltaCause1[i]==1)
      {
        # \int_{0}^{Ti}\lambda_{01}(s)\exp\{ a1*m_{i}(t) \}
        # at bican
        ints1can = exp( ui1Long + ( LxibetasLong + c(ziLong %*% bican) )*assoc1 )
        s1bcan = 0.5*times[i]*sum( wws*ints1can )
        a1s1bcan = 1 + alpha1*s1bcan
        
        # \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \}
        # at bican
        ints2can = exp( ui2Long + ( LxibetasLong + c(ziLong %*% bican) )*assoc2 )
        s2bcan = 0.5*times[i]*sum( wws*ints2can )
        a2s2bcan = 1 + alpha2*s2bcan
        
        # \int_{0}^{Ti}\lambda_{01}(s)\exp\{ a1*m_{i}(t) \}
        # at bi
        ints1 = exp( ui1Long + ( LxibetasLong + c(ziLong %*% bi) )*assoc1 )
        s1b = 0.5*times[i]*sum( wws*ints1 )
        a1s1b = 1 + alpha1*s1b
        
        logp = logp + assoc1*sum(bdiff*zTime[i,]) - (1+alpha1)*log(a1s1bcan/a1s1b)/alpha1 - sum(bdiff*Linfmu1)
        
        # Set the probability of acceptance equal to zero if the sum of the CIFs is not bounded by 1
        st = a1s1bcan^(-1/alpha1) + a2s2bcan^(-1/alpha2) - 1
        if ( st < 0)
        {
          logp = - Inf
        }
        
        
      }
      
      if (deltaCause2[i]==1)
      {
        # \int_{0}^{Ti}\lambda_{01}(s)\exp\{ a1*m_{i}(t) \}
        # at bican
        ints1can = exp( ui1Long + ( LxibetasLong + c(ziLong %*% bican) )*assoc1 )
        s1bcan = 0.5*times[i]*sum( wws*ints1can )
        a1s1bcan = 1 + alpha1*s1bcan
        
        # \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \}
        # at bican
        ints2can = exp( ui2Long + ( LxibetasLong + c(ziLong %*% bican) )*assoc2 )
        s2bcan = 0.5*times[i]*sum( wws*ints2can )
        a2s2bcan = 1 + alpha2*s2bcan
        
        # at bi
        ints2 = exp( ui2Long + ( LxibetasLong + c(ziLong %*% bi) )*assoc2 )
        s2b = 0.5*times[i]*sum( wws*ints2 )
        a2s2b = 1 + alpha2*s2b
        
        logp = logp + assoc2*sum(bdiff*zTime[i,]) - (1+alpha2)*log(a2s2bcan/a2s2b)/alpha2 - sum(bdiff*Linfmu1)
        
        # Set the probability of acceptance equal to zero if the sum of the CIFs is not bounded by 1
        st = a1s1bcan^(-1/alpha1) + a2s2bcan^(-1/alpha2) - 1
        if ( st < 0)
        {
          logp = - Inf
        }
        
      }
      
      if (delta[i]==0)
      {
        # \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \}
        ints1can = exp( ui1Long + ( LxibetasLong + c(ziLong %*% bican) )*assoc1 )
        ints1 = exp( ui1Long + ( LxibetasLong + c(ziLong %*% bi) )*assoc1 )
        
        # \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \}
        ints2can = exp( ui2Long + ( LxibetasLong + c(ziLong %*% bican) )*assoc2 )
        ints2 = exp( ui2Long + ( LxibetasLong + c(ziLong %*% bi) )*assoc2 )
        
        # \int_{0}^{Ti}\lambda_{01}(s)\exp\{ a1*m_{i}(t) \}
        s1bcan = 0.5*times[i]*sum( wws*ints1can )
        a1s1bcan = 1 + alpha1*s1bcan
        
        s1b = 0.5*times[i]*sum( wws*ints1 )
        a1s1b = 1 + alpha1*s1b
        
        # \int_{0}^{Ti}\lambda_{02}(s)\exp\{ a2*m_{i}(t) \}
        s2bcan = 0.5*times[i]*sum( wws*ints2can )
        a2s2bcan = 1 + alpha2*s2bcan
        
        s2b = 0.5*times[i]*sum( wws*ints2 )
        a2s2b = 1 + alpha2*s2b
        
        logp = logp + log(a1s1bcan^(-1/alpha1) + a2s2bcan^(-1/alpha2) - 1) - log(a1s1b^(-1/alpha1) + a2s2b^(-1/alpha2) - 1) - sum(bdiff*Linfmu1)
      }
      
      # Probability of acceptance
      unif = runif(1)
      
      if (!is.finite(logp))
      {
        logp = -Inf
      }
      
      if (log(unif)<logp)
      {
        b[i,] = bican
        accept_re[i] = accept_re[i] + 1
      }
    }
    
    # Adjusted dependent variable: y - z*b
    Lz_blong = .Internal(t.default(Lzmatrix*b[id,]))
    ytilde = y - .Internal(colSums(Lz_blong, q, n, FALSE))
    
    ############################################################
    ### Update the design matrices of the survival submodels ###
    ############################################################
    
    # Predictions of the marker value at the survival time of each subject
    zTimeb = rowSums(zTime*b)
    mtpred = timesLbeta + zTimeb
    
    # Predictions of the marker value at additional points (required to integrate the subdistribution hazard)
    zTimeLongb = rowSums(zTimeLong*b[ind,])
    mtpredLong = timesLbetaLong + zTimeLongb
    
    # Design matrices of the survival submodel
    if (lbetasCause1==1)
    {
      xCause1 = mtpred
      xCause1Long = mtpredLong
    } else{
      xCause1[,lbetasCause1] = mtpred
      xCause1Long[,lbetasCause1] = mtpredLong
    }
    
    if (lbetasCause2==1)
    {
      xCause2 = mtpred
      xCause2Long = mtpredLong
    } else{
      xCause2[,lbetasCause1] = mtpred
      xCause2Long[,lbetasCause1] = mtpredLong
    }
    
    # Combined design matrices of the survival submodel (including the B-splines basis)
    xCause1Comb = cbind(xCause1,splBasisCause1)
    xCause1CombLong = cbind(xCause1Long,splBasisCause1Long)
    
    xCause2Comb = cbind(xCause2,splBasisCause2)
    xCause2CombLong = cbind(xCause2Long,splBasisCause2Long)
    
    ############################################################
    ### Update the design matrices of the survival submodels ###
    ############################################################
    
    # Predictions of the marker value at the survival time of each subject
    zTimeb = rowSums(zTime*b)
    mtpred = timesLbeta + zTimeb
    
    # Predictions of the marker value at additional points (required to integrate the subdistribution hazard)
    zTimeLongb = rowSums(zTimeLong*b[ind,])
    mtpredLong = timesLbetaLong + zTimeLongb
    
    # Design matrices of the survival submodel
    if (lbetasCause1==1)
    {
      xCause1 = mtpred
      xCause1Long = mtpredLong
    } else{
      xCause1[,lbetasCause1] = mtpred
      xCause1Long[,lbetasCause1] = mtpredLong
    }
    
    if (lbetasCause2==1)
    {
      xCause2 = mtpred
      xCause2Long = mtpredLong
    } else{
      xCause2[,lbetasCause1] = mtpred
      xCause2Long[,lbetasCause1] = mtpredLong
    }
    
    # Combined design matrices of the survival submodel (including the B-splines basis)
    xCause1Comb = cbind(xCause1,splBasisCause1)
    xCause1CombLong = cbind(xCause1Long,splBasisCause1Long)
    
    xCause2Comb = cbind(xCause2,splBasisCause2)
    xCause2CombLong = cbind(xCause2Long,splBasisCause2Long)
    
    ##################################################
    ### Update the survival parameters for cause 1 ###
    ##################################################
    
    # Starting from the current value of the chain, perform BFGS updates
    if (iterMaxSurv_1==1)
    {
      mu = paramCause1
      gradient.t = fDer1(mu) + c(Sinvc0_1 %*% (mu - Smu0_1))
      mu.t = mu - solve(hessStart1,gradient.t)
    } else {
      mu = paramCause1
      gradient.t = fDer1(mu) + c(Sinvc0_1 %*% (mu - Smu0_1))
      approxInf = hessStart1
      
      mu.t = mu
      
      for (j in 1:iterMaxSurv_1)
      {
        mu = mu.t
        mu.t = try(mu - solve(approxInf,gradient.t),silent = T)
        if (class(mu.t) %in% "try-error")
        {
          mu.t = mu - solve(hessStart1,gradient.t)
        }
        #print(mu.t)
        
        if (j<iterMaxSurv_1)
        {
          # Evaluate the gradient
          gradient = gradient.t
          gradient.t = fDer1(mu.t) + c(Sinvc0_1 %*% (mu.t - Smu0_1))
          
          s = (mu.t - mu)
          psi = (gradient.t - gradient)
          
          # Update the Hessian matrix
          approxInf = approxInf + tcrossprod(psi)/sum(psi*s) - (approxInf%*%tcrossprod(s)%*%approxInf)/c(crossprod(s,approxInf) %*% s)
        }
        
      }
    }
    
    paramCause1.can = c(rmvt(n = 1,sigma = Sigmastart1,delta = mu.t,df = tdf))
    
    # Calculate posterior ratio
    logp = fsurvEval1(paramCause1.can) - fsurvEval1(paramCause1) 
    logp = logp + dmvnorm(paramCause1.can,mean = Smu0_1,sigma = Sc0_1,log = T) - dmvnorm(paramCause1,mean = Smu0_1,sigma = Sc0_1,log = T)
    
    # Proposal ratio
    logp = logp - dmvt(paramCause1.can,delta = mu.t,sigma = Sigmastart1,df = tdf,log = T)
    
    # Starting from the proposed value, perform BFGS updates
    if (iterMaxSurv_1==1)
    {
      mu = paramCause1.can
      gradient.t = fDer1(mu) + c(Sinvc0_1 %*% (mu - Smu0_1))
      mu.t = mu - solve(hessStart1,gradient.t)
    } else {
      mu = paramCause1.can
      gradient.t = fDer1(mu) + c(Sinvc0_1 %*% (mu - Smu0_1))
      approxInf = hessStart1
      
      mu.t = mu
      
      for (j in 1:iterMaxSurv_1)
      {
        mu = mu.t
        mu.t = try(mu - solve(approxInf,gradient.t),silent = T)
        if (class(mu.t) %in% "try-error")
        {
          mu.t = mu - solve(hessStart1,gradient.t)
        }
        #print(mu.t)
        
        if (j<iterMaxSurv_1)
        {
          # Evaluate the gradient
          gradient = gradient.t
          gradient.t = fDer1(mu.t) + c(Sinvc0_1 %*% (mu.t - Smu0_1))
          
          s = (mu.t - mu)
          psi = (gradient.t - gradient)
          
          # Update the Hessian matrix
          approxInf = approxInf + tcrossprod(psi)/sum(psi*s) - (approxInf%*%tcrossprod(s)%*%approxInf)/c(crossprod(s,approxInf) %*% s)
        }
        
      }
    }
    
    logp = logp + dmvt(paramCause1,delta = mu.t,sigma = Sigmastart1,df = tdf,log = T)
    
    # Draw u~U(0,1) to decide
    u = runif(1)
    
    if (!is.finite(logp))
    {
      logp = -Inf
    }
    
    if (log(u)<logp)
    {
      paramCause1 = paramCause1.can
      accept_surv1 = accept_surv1 + 1
    }
    
    ##################################################
    ### Update the survival parameters for cause 2 ###
    ##################################################
    
    # Starting from the current value of the chain, perform BFGS updates
    if (iterMaxSurv_2==1)
    {
      mu = paramCause2
      gradient.t = fDer2(mu) + c(Sinvc0_2 %*% (mu - Smu0_2))
      mu.t = mu - solve(hessStart2,gradient.t)
    } else {
      mu = paramCause2
      gradient.t = fDer2(mu) + c(Sinvc0_2 %*% (mu - Smu0_2))
      approxInf = hessStart2
      
      mu.t = mu
      
      for (j in 1:iterMaxSurv_2)
      {
        mu = mu.t
        mu.t = try(mu - solve(approxInf,gradient.t),silent = T)
        if (class(mu.t) %in% "try-error")
        {
          mu.t = mu - solve(hessStart2,gradient.t)
        }
        #print(mu.t)
        
        if (j<iterMaxSurv_2)
        {
          # Evaluate the gradient
          gradient = gradient.t
          gradient.t = fDer2(mu.t) + c(Sinvc0_2 %*% (mu.t - Smu0_2))
          
          s = (mu.t - mu)
          psi = (gradient.t - gradient)
          
          # Update the Hessian matrix
          approxInf = approxInf + tcrossprod(psi)/sum(psi*s) - (approxInf%*%tcrossprod(s)%*%approxInf)/c(crossprod(s,approxInf) %*% s)
        }
        
      }
    }
    
    paramCause2.can = c(rmvt(n = 1,sigma = Sigmastart2,delta = mu.t,df = tdf))
    
    # Calculate posterior ratio
    logp = fsurvEval2(paramCause2.can) - fsurvEval2(paramCause2) 
    logp = logp + dmvnorm(paramCause2.can,mean = Smu0_2,sigma = Sc0_2,log = T) - dmvnorm(paramCause2,mean = Smu0_2,sigma = Sc0_2,log = T)
    
    # Proposal ratio
    logp = logp - dmvt(paramCause2.can,delta = mu.t,sigma = Sigmastart2,df = tdf,log = T)
    
    # Starting from the proposed value, perform BFGS updates
    if (iterMaxSurv_2==1)
    {
      mu = paramCause2.can
      gradient.t = fDer2(mu) + c(Sinvc0_2 %*% (mu - Smu0_2))
      mu.t = mu - solve(hessStart2,gradient.t)
    } else {
      mu = paramCause2.can
      gradient.t = fDer2(mu) + c(Sinvc0_2 %*% (mu - Smu0_2))
      approxInf = hessStart2
      
      mu.t = mu
      
      for (j in 1:iterMaxSurv_2)
      {
        mu = mu.t
        mu.t = try(mu - solve(approxInf,gradient.t),silent = T)
        if (class(mu.t) %in% "try-error")
        {
          mu.t = mu - solve(hessStart2,gradient.t)
        }
        #print(mu.t)
        
        if (j<iterMaxSurv_2)
        {
          # Evaluate the gradient
          gradient = gradient.t
          gradient.t = fDer2(mu.t) + c(Sinvc0_2 %*% (mu.t - Smu0_2))
          
          s = (mu.t - mu)
          psi = (gradient.t - gradient)
          
          # Update the Hessian matrix
          approxInf = approxInf + tcrossprod(psi)/sum(psi*s) - (approxInf%*%tcrossprod(s)%*%approxInf)/c(crossprod(s,approxInf) %*% s)
        }
        
      }
    }
    
    logp = logp + dmvt(paramCause2,delta = mu.t,sigma = Sigmastart2,df = tdf,log = T)
    
    # Draw u~U(0,1) to decide
    u = runif(1)
    
    if (!is.finite(logp))
    {
      logp = -Inf
    }
    
    if (log(u)<logp)
    {
      paramCause2 = paramCause2.can
      accept_surv2 = accept_surv2 + 1
    }
    
    # paramSurv is the updated survival model parameter vector
    paramSurv = c(paramCause1,paramCause2)
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    assoc1 = betasCause1[lbetasCause1]
    assoc2 = betasCause2[lbetasCause2]
    #rm(paramCause1,paramCause2)
    
    if (lbetasCause1==1)
    {
      uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    } else{
      uCause1Long = c(splBasisCause1Long %*% gammasCause1 + as.matrix(xCause1Long[,-lbetasCause1]) %*% betasCause1[-lbetasCause1] )
    }
    
    if (lbetasCause2==1)
    {
      uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    } else{
      uCause2Long = c(splBasisCause2Long %*% gammasCause2 + as.matrix(xCause2Long[,-lbetasCause2]) %*% betasCause2[-lbetasCause2] )
    }
    
    # Update the missing failure indicators
    # For the non-doubly sampled individuals, simulate the ``true'' failure cause conditional on all the observed information and 
    # the current parameter values
    temp = fsurvEvalRdouble(paramSurv)
    logf12 = temp$DercincCause1 - temp$DercincCause2
    
    miscProb = log(pi11/(1 - pi22))*deltaCause1Obs + log((1 - pi11)/pi22)*deltaCause2Obs
    miscProb = miscProb[Rmis]
    
    prob1 = exp(logf12[Rmis] + miscProb)
    prob1 = prob1/(1+prob1)
    
    deltaCause1[Rmis] = rbinom(nRmis,size = 1,prob = prob1)
    deltaCause2[Rmis] = 1 - deltaCause1[Rmis]
    
    deltaCause1Long = rep( deltaCause1, each = length(pps) )
    deltaCause2Long = rep( deltaCause2, each = length(pps) )
    deltaLong = deltaCause1Long + deltaCause2Long
    
    # Failure indicator
    delta = deltaCause1 + deltaCause2
    indcens = which(delta == 0)
    indDelta1 = which(deltaCause1==1)
    indDelta2 = which(deltaCause2==1)
    
    # Update the misclassification parameters
    pi11 = rbeta(1,shape1 = a1 + sum(deltaCause1*deltaCause1Obs),shape2 = b1 + sum(deltaCause1*deltaCause2Obs) )
    pi22 = rbeta(1,shape1 = a2 + sum(deltaCause2*deltaCause2Obs),shape2 = b2 + sum(deltaCause2*deltaCause1Obs) )
    
    ###################
    ### Update beta ###
    ###################
    
    # C_{0}^{-1}mu_{0} + omega*sum x_{i}^{T}(y_{i}-z_{i}b_{i})
    right = .Internal(crossprod(Lx,ytilde))*omega + Linvc0mu0
    
    # Information matrix of beta based on the longitudinal model
    infBetaLong = Linvc0 + omega*Lxtx
    LinfBetaLong = .Internal(La_chol(infBetaLong, pivot, tol))
    
    mu = c(.Internal( backsolve(LinfBetaLong, .Internal(backsolve(LinfBetaLong,right,Lp,T,T)) ,Lp,T,F)  ))
    #gradient.t = fbetaDer(mu)
    approxInf = infBetaLong
    
    mu.t = mu
    
    # for (j in 1:1)
    # {
    #   mu = mu.t
    #   mu.t = mu - solve(approxInf,gradient.t)
    #   #print(mu.t)
    #   
    #   # Evaluate the gradient
    #   gradient = gradient.t
    #   gradient.t = fsurvDer(mu.t)
    #   
    #   s = (mu.t - mu)
    #   psi = (gradient.t - gradient)
    #   
    #   # Update the Hessian matrix
    #   approxInf = approxInf + tcrossprod(psi)/sum(psi*s) - (approxInf%*%tcrossprod(s)%*%approxInf)/c(crossprod(s,approxInf) %*% s)
    #   
    # }
    
    epsilon = rnorm(Lp)
    Lbetacan = c(.Internal(backsolve(LinfBetaLong,epsilon,Lp,TRUE,FALSE))) + mu.t
    
    logp = fbetaEval(Lbetacan) - fbetaEval(Lbeta) 
    #logp = logp + dmvnorm(Lbeta,mean = mu.t,sigma = solve(infBetaLong),log = T) - dmvnorm(Lbetacan,mean = mu.t,sigma = solve(infBetaLong),log = T)
    
    logp = logp + 0.5*sum( (LinfBetaLong %*% (Lbetacan - mu.t) )^2 ) - 0.5*sum( (LinfBetaLong %*% (Lbeta - mu.t) )^2 )
    
    # If the proposed value of Lbeta does not meet the constraint,
    # set the acceptance probability equal to zero
    if (!is.finite(logp))
    {
      logp = - Inf
    }
    
    u = runif(1)
    
    if (log(u)<logp)
    {
      Lbeta = Lbetacan
      accept_Lbeta = accept_Lbeta + 1
    }
    
    # We need to update all quantities involving Lbeta
    Lxbeta = c(Lx %*% Lbeta)
    
    # Populations residuals
    Lres = y - Lxbeta
    
    # Population marker values at the survival times
    timesLbeta = xTime %*% Lbeta
    timesLbetaLong = xTimeLong %*% Lbeta
    
    # Update omega
    omega = rgamma(1,shape = 0.5*n + lambda1, rate = lambba2 + 0.5*sum((ytilde - Lxbeta)^2) )
    
    # z_{i}^{T}(y_{i}-x_{i}\beta)*omega
    ztres = as.vector(ztbdiag %*% Lres)*omega
    ztres = .Internal( matrix(ztres, nrow = G, ncol = q, byrow = TRUE, dimnames = NULL,FALSE, FALSE ) )
    
    # Update Dinv
    Sigma = A + .Internal(crossprod(b,NULL))
    Sigma = .Internal(La_chol2inv( .Internal(La_chol(Sigma,pivot,tol)),q))
    
    Dinv = rWishart(1, G + df, Sigma)[,,1]
    
    if (it>0)
    {
      Ldraws[it,] = c(Lbeta,1/omega)
      Sdraws[it,] = paramSurv
      Mdraws[it,] = c(pi11,pi22)
      print(it)
      # Store the covariance matrix
      drawsD[,,it] = Dinv
      
      if (store.b == T)
      {
        draws.b[,,it] = b 
      }
    }
    
  } #MCMC loop ends
  
  for (i in 1:(ndraw*thin))
  {
    drawsD[,,i] = .Internal(La_chol2inv( .Internal(La_chol(drawsD[,,i],pivot,tol)),q))
  }
  
  # Prepare output
  Ldraws = Ldraws[(1:ndraw)*thin,]
  Sdraws = Sdraws[(1:ndraw)*thin,]
  Mdraws = Mdraws[(1:ndraw)*thin,]
  drawsD = drawsD[,,(1:ndraw)*thin]
  if (store.b == T)
  {
    draws.b = draws.b[,,(1:ndraw)*thin]
  }
  
  # Output list
  out = list()
  
  out$accept_Lbeta = round(100*accept_Lbeta/(ndraw*thin + nburn), 1)
  out$accept_re = round(100*accept_re/(ndraw*thin + nburn), 1)
  out$accept_surv1 = round(100*accept_surv1/(ndraw*thin + nburn), 1)
  out$accept_surv2 = round(100*accept_surv2/(ndraw*thin + nburn), 1)
  
  # Posterior samples
  out$Sdraws = Sdraws
  out$Ldraws = Ldraws
  out$Mdraws = Mdraws
  if (store.b==T)
  {
    out$draws.b = draws.b
  }
  
  sum = list()
  
  # Longitudinal parameters
  sumLong = matrix(nrow = Lp + 1,ncol = 6)
  rownames(sumLong) = c(names.Lbeta,"Variance")
  colnames(sumLong) = c("Post.mean","Post.sd","MC.se","2.5%","50%","97.5%")
  
  sumLong[,1] = colMeans(Ldraws)
  sumLong[,2] = apply(Ldraws,2,sd)
  sumLong[,3] = sumLong[,2]/sqrt(ndraw)
  sumLong[,4:6] = t(apply(Ldraws,2,quantile,c(0.025,0.5,0.975),na.rm = TRUE))
  
  # Longitudinal parameters
  sumMiscl = matrix(nrow = 2,ncol = 6)
  rownames(sumMiscl) = c("Pr(Kobs=1|Ktrue=1)","Pr(Kobs=2|Ktrue=2)")
  colnames(sumMiscl) = c("Post.mean","Post.sd","MC.se","2.5%","50%","97.5%")
  
  sumMiscl[,1] = colMeans(Mdraws)
  sumMiscl[,2] = apply(Mdraws,2,sd)
  sumMiscl[,3] = sumMiscl[,2]/sqrt(ndraw)
  sumMiscl[,4:6] = t(apply(Mdraws,2,quantile,c(0.025,0.5,0.975),na.rm = TRUE))
  
  sum$Longitudinal_Process = sumLong
  sum$sumMiscl = sumMiscl
  
  # Survival parameters
  sumSurvCause1 = matrix(nrow = lpCause1,ncol = 6)
  rownames(sumSurvCause1) = c(names(coef(fitCoxCause1)),"Assoc1",paste0("GammasCause1_",1:lgammasCause1))
  colnames(sumSurvCause1) = c("Post.mean","Post.sd","MC.se","2.5%","50%","97.5%")
  
  sumSurvCause1[,1] = colMeans(Sdraws[,1:lpCause1])
  sumSurvCause1[,2] = apply(Sdraws[,1:lpCause1],2,sd)
  sumSurvCause1[,3] = sumSurvCause1[,2]/sqrt(ndraw)
  sumSurvCause1[,4:6] = t(apply(Sdraws[,1:lpCause1],2,quantile,c(0.025,0.5,0.975),na.rm = TRUE))
  
  sumSurvCause2 = matrix(nrow = lpCause2,ncol = 6)
  rownames(sumSurvCause2) = c(names(coef(fitCoxCause2)),"Assoc2",paste0("GammasCause2_",1:lgammasCause2))
  colnames(sumSurvCause2) = c("Post.mean","Post.sd","MC.se","2.5%","50%","97.5%")
  
  sumSurvCause2[,1] = colMeans(Sdraws[,-c(1:lpCause1)])
  sumSurvCause2[,2] = apply(Sdraws[,-c(1:lpCause1)],2,sd)
  sumSurvCause2[,3] = sumSurvCause2[,2]/sqrt(ndraw)
  sumSurvCause2[,4:6] = t(apply(Sdraws[,-c(1:lpCause1)],2,quantile,c(0.025,0.5,0.975),na.rm = TRUE))
  
  sum$sumSurvCause1 = sumSurvCause1
  sum$sumSurvCause2 = sumSurvCause2
  
  # Covariance matrix of the random effects
  out$drawsD = drawsD
  
  sum$VarCovRE_mean = apply(drawsD,1:2,mean)
  sum$VarCovRE_median = apply(drawsD,1:2,median)
  sum$VarCovRE_sd = apply(drawsD,1:2,sd)
  
  drawsCorD = drawsD
  for (i in 1:ndraw)
  {
    drawsCorD[,,i] = cov2cor(drawsD[,,i])
  }
  
  sum$CorrRE_mean = apply(drawsCorD,1:2,mean)
  sum$CorrRE_median = apply(drawsCorD,1:2,median)
  sum$CorrRE_sd = apply(drawsCorD,1:2,sd)
  
  out$drawsCorD = drawsCorD
  out$summary = sum
  out$fitStart = fitStart
  out$cox1 = fitCoxCause1
  out$cox2 = fitCoxCause2
  out$fitlme = fitlme
  out$knotsCause1 = knotsCause1
  out$knotsCause2 = knotsCause2
  out$bknotsCause1 = bknotsCause1
  out$bknotsCause2 = bknotsCause2
  out$lpCause1 = lpCause1
  out$lpCause2 = lpCause2
  out$lbetasCause1 = lbetasCause1
  out$lbetasCause2 = lbetasCause2
  out$lgammasCause1 = lgammasCause1
  out$lgammasCause2 = lgammasCause2
  out$timeVar = timeVar
  out$alpha1 = alpha1
  out$alpha2 = alpha2
  out$prior = prior
  out$useGauleg = useGauleg
  out$Gauleg_points = Gauleg_points
  
  out$xTime = xTime
  out$xTimeLong = xTimeLong
  out$zTime = zTime
  out$zTimeLong = zTimeLong
  out$ylist = ylist
  out$xlist = xlist
  out$zlist = zlist
  out$zzlist = zzlist
  out$ni = ni
  out$G = G
  out$splBasisCause1 = splBasisCause1
  out$splBasisCause1Long = splBasisCause1Long
  out$splBasisCause2 = splBasisCause2
  out$splBasisCause2Long = splBasisCause2Long
  out$Lp = Lp 
  out$q = q
  out$lpCause1 = lpCause1 
  out$lbetasCause1 = lbetasCause1
  out$lbetasCause2 = lbetasCause2
  out$ids = ids
  out$deltaCause1 = deltaCause1
  out$deltaCause2 = deltaCause2
  out$xCause1 = xCause1
  out$xCause2 = xCause2
  out$times = times
  out$Rdouble = Rdouble
  out$deltaCause1Obs = deltaCause1Obs
  out$deltaCause2Obs = deltaCause2Obs
  
  return(out)
}

DICIncrCumInc_Misc = function(fitProp,nMC = 100,thinDIC = 10)
{
  # 7-15 Gauss-Kronrod quadradure points and weights
  gaussKronrod <-
    function (k = 15) {
      sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
              0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
              -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
              0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
      wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
      wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975, 
               0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
      if (k == 7) 
        list(sk = sk[1:7], wk = wk7)
      else
        list(sk = sk, wk = wk15)
    }
  
  gauleg <- function(n,a=-1,b=1) {# Gauss-Legendre: returns x,w so that
    #\int_a^b f(x) dx \doteq \sum w_i f(x_i)
    EPS <- 3.e-14
    m <- trunc((n+1)/2)
    xm <- 0.5*(b+a)
    xl <- 0.5*(b-a)
    x <- w <- rep(-1,n)
    for (i in 1:m) {
      z <- cos(pi*(i-.25)/(n+.5))
      z1 <- z+1
      while (abs(z-z1) > EPS) {
        p1 <- 1
        p2 <- 0
        for (j in 1:n) {# recursively evaluates pn(x)
          p3 <- p2
          p2 <- p1
          p1 <- ((2*j-1)*z*p2-(j-1)*p3)/j
        }
        pp <- n*(z*p1-p2)/(z*z-1)
        z1 <- z
        z <- z1-p1/pp #Newton iteration
      }
      x[i] <- xm-xl*z
      x[n+1-i] <- xm+xl*z
      w[i] <- 2*xl/((1-z*z)*pp*pp)
      w[n+1-i] <- w[i]
    }
    list(x=x,w=w)
  }
  
  useGauleg = fitProp$useGauleg
  Gauleg_points = fitProp$Gauleg_points
  
  if (useGauleg == F)
  {
    gk = gaussKronrod()
    pps = gk$sk
    wws = gk$wk
  } else {
    gk = gauleg(Gauleg_points)
    pps = gk$x
    wws = gk$w
  }
  
  # A function approximating the observed loglikelihood of the model 
  f = function(param)
  {
    # Fixed effects
    Lbeta = param[1:Lp]
    
    # Level-1 vairiance
    sigma2 = exp(param[Lp+1])
    
    # Covariance matrix of the random effect
    vechL = param[(Lp+1+1):(Lp+1+q*(q+1)/2)]
    Lt[lower.tri(Lt,diag = T)] = vechL
    L = t(Lt)
    D = crossprod(L)
    Dinv = solve(D)
    
    # Parameters of the survival model
    param.surv = param[(1+Lp+1+q*(q+1)/2):(lpCause1 + lpCause2 + Lp+1+q*(q+1)/2)]
    
    paramCause1 = param.surv[1:lpCause1]
    paramCause2 = param.surv[-c(1:lpCause1)]
    
    # Parameters for cause 1
    betasCause1 = paramCause1[1:lbetasCause1]
    gammasCause1 = paramCause1[-c(1:lbetasCause1)]
    
    # Parameters for cause 2
    betasCause2 = paramCause2[1:lbetasCause2]
    gammasCause2 = paramCause2[-c(1:lbetasCause2)]
    
    # Misclassification parameters
    pi11 = param[-c(1:(lpCause1 + lpCause2 + Lp+1+q*(q+1)/2))][1]
    pi22 = param[-c(1:(lpCause1 + lpCause2 + Lp+1+q*(q+1)/2))][2]
    
    # u_{k}(t)
    uCause1 = c(splBasisCause1 %*% gammasCause1)
    uCause1Long = c(splBasisCause1Long %*% gammasCause1)
    
    uCause2 = c(splBasisCause2 %*% gammasCause2)
    uCause2Long = c(splBasisCause2Long %*% gammasCause2)
    
    timesLbeta = c(xTime %*% Lbeta)
    timesLbetaLong = c(xTimeLong %*% Lbeta)
    
    # Function evaluating the survival likelihood for the ith subject
    # for each row of a matrix including values of the random effects
    fsurv = function(br)
    {
      # Transpose of the random effects matrix
      tbr = t(br)
      
      # The true marker value at the additional time points for each value of the random effects for the ith subject
      # Each row corresponds to the additional time points
      # and each column coorresponds to the values of the random effects
      mtimeLong = timesLbetaLong[ids[[i]]] + zTimeLong[ids[[i]],] %*% tbr
      mtime = timesLbeta[i] + br %*% ztTime[,i]
      
      # Cause 1
      if (lbetasCause1==1)
      {
        xbetasCause1 = c(mtime*betasCause1[lbetasCause1])
        xbetasCause1Long = mtimeLong*betasCause1[lbetasCause1]
      } else {
        # Not working!!!
        xbetasCause1 = c(mtime*betasCause1[lbetasCause1])
        xbetasCause1Long = mtimeLong*betasCause1[lbetasCause1]
      }
      
      linpredCause1 = xbetasCause1 + uCause1[i] + sum(xCause1[i,-lbetasCause1]*betasCause1[-lbetasCause1])
      linpredCause1Long = xbetasCause1Long + uCause1Long[ids[[i]]] + sum(xCause1[i,-lbetasCause1]*betasCause1[-lbetasCause1])
      explinpredCause1Long = exp(linpredCause1Long)
      int1 = 0.5*times[i]*colSums(wws*explinpredCause1Long)
      
      # Cause 2
      if (lbetasCause2==1)
      {
        xbetasCause2 = c(mtime*betasCause2[lbetasCause2])
        xbetasCause2Long = mtimeLong*betasCause2[lbetasCause2]
      } else {
        # Not working!!!
        xbetasCause2 = c(mtime*betasCause2[lbetasCause2])
        xbetasCause2Long = mtimeLong*betasCause2[lbetasCause2]
      }
      
      linpredCause2 = xbetasCause2 + uCause2[i] + sum(xCause2[i,-lbetasCause2]*betasCause2[-lbetasCause2])
      linpredCause2Long = xbetasCause2Long + uCause2Long[ids[[i]]] + sum(xCause2[i,-lbetasCause2]*betasCause2[-lbetasCause2])
      explinpredCause2Long = exp(linpredCause2Long)
      int2 = 0.5*times[i]*colSums(wws*explinpredCause2Long)
      
      # Cumulative incidence functions
      cincCause1 = 1 - ( 1 + alpha1*int1 )^(-1/alpha1)
      cincCause2 = 1 - ( 1 + alpha2*int2 )^(-1/alpha2)
      
      # Derivatives of the cumulative incidence functions (on the log scale)
      DercincCause1 = -((alpha1+1)/alpha1)*log(1 + alpha1*int1) + linpredCause1
      DercincCause2 = -((alpha2+1)/alpha2)*log(1 + alpha2*int2) + linpredCause2
      eDercincCause1 = exp(DercincCause1)
      eDercincCause2 = exp(DercincCause2)
      
      # Overall survival function
      st = 1 - cincCause1 - cincCause2
      
      # For those who failed but are not included in the double sampling
      if (length(which(Rdouble[i] == 0)) > 0)
      {
        out = deltaCause1Obs[i]*( eDercincCause1*pi11 + eDercincCause2*(1-pi22) )
        out = out + deltaCause2Obs[i]*( eDercincCause1*(1-pi11) + eDercincCause2*pi22 )
        out = log(out)
      }
      
      # For those included in the double sampling
      if (length(which(Rdouble[i] == 1)) > 0)
      {
        if (deltaCause1[i]==1)
        {
          out = DercincCause1 + log(pi11)*deltaCause1Obs[i] + log(1-pi11)*deltaCause2Obs[i]
        }
        
        if (deltaCause2[i]==1)
        {
          out = DercincCause2 + log(1-pi22)*deltaCause1Obs[i] + log(pi22)*deltaCause2Obs[i]
        }
      }
      
      if (deltaCause1[i] == 0 & deltaCause2[i] == 0)
      {
        # Cumulative incidence functions
        out = log(st)
      }
      
      jj = which(st > 0)
      res = rep(0,nrow(br))
      res[jj] = exp(out[jj])
      res
    }
    
    
    for (i in 1:G)
    {
      # Matrices/vectors of the ith subject
      yi = ylist[[i]]
      xi = xlist[[i]]
      zi = zlist[[i]]
      zzi = zzlist[[i]]
      nii = ni[i]
      resi = yi - c(xi %*% Lbeta)
      
      # Covariance matrix for the ith subject
      Vi = diag(nii)*sigma2 + tcrossprod(zi %*% Lt )
      Ui = .Internal(La_chol(Vi, pivot, tol))
      
      # t(yi-xi%*%Lbeta) %*% solve(Vi) %*% (yi-xi%*%Lbeta)
      Uim1tresi = .Internal(backsolve(Ui,resi,nii,TRUE,TRUE))
      restWres = sum(Uim1tresi^2)
      
      # Log-determinant of Vi
      logDetVi = 2*sum(log(diag(Ui)))
      
      # Posterior distribution of b_{i}|y_{i}
      Linf = Dinv + zzi/sigma2
      Bi = .Internal(La_chol(Linf, pivot, tol))
      temp = crossprod(zi,resi/sigma2)
      mub = c(.Internal( backsolve(Bi, .Internal(backsolve(Bi,temp,q,TRUE,TRUE)) ,q,TRUE,FALSE )  ))
      
      # Decide between mub and mubBar[i,]
      # val1 = (fsurv(rbind(mubBar[i,],mubBar[i,]))*dmvnorm(rbind(mubBar[i,],mubBar[i,]),mean = mub,sigma = solve(Linf),log = F))[1]
      # val2 = (fsurv(rbind(mub,mub))*dmvnorm(rbind(mub,mub),mean = mub,sigma = solve(Linf),log = F))[1]
      
      val1 = fsurv(rbind(mubBar[i,],mubBar[i,]))[1]*exp(-0.5*sum( (Bi %*% (mubBar[i,]-mub))^2 ) )
      val2 = fsurv(rbind(mub,mub))[1]
      
      
      if (val1>=val2)
      {
        #br = mnorm.prec(nMC,mu = mubBar[i,],prec = Linf)
        br = mnorm.CholPrec(nMC,mu = mubBar[i,],L = Bi)
        
        # Posterior values
        #val = fsurv(br)*dmvnorm(br,mean = mub,sigma = solve(Linf),log = F)/dmvnorm(br,mean = mubBar[i,],sigma = solve(Linf),log = F)
        val = fsurv(br)*exp( - 0.5*.colSums( (Bi %*% (t(br)-mub))^2 - (Bi %*% (t(br)-mubBar[i,]))^2,q,nMC,F) )
        
        if (mean(val>0)<0.05)
        {
          nMC = nMC*15
          #br = mnorm.prec(nMC,mu = mubBar[i,],prec = Linf)
          br = mnorm.CholPrec(nMC,mu = mubBar[i,],L = Bi)
          
          #val = fsurv(br)*dmvnorm(br,mean = mub,sigma = solve(Linf),log = F)/dmvnorm(br,mean = mubBar[i,],sigma = solve(Linf),log = F)
          val = fsurv(br)*exp( - 0.5*.colSums( (Bi %*% (t(br)-mub))^2 - (Bi %*% (t(br)-mubBar[i,]))^2,q,nMC,F) )
          
          nMC = nMCbup
        }
      } else{
        #br = mnorm.prec(nMC,mu = mub,prec = Linf)
        br = mnorm.CholPrec(nMC,mu = mub,L = Bi)
        
        # Posterior values
        val = fsurv(br)
        
        if (mean(val>0)<0.05)
        {
          nMC = nMC*15
          br = mnorm.CholPrec(nMC,mu = mub,L = Bi)
          val = fsurv(br)
          nMC = nMCbup
        }
      }
      
      # Survival contribution
      logSurv = log(mean(val))
      
      logl[i] = -0.5*(restWres + logDetVi + nii*log(2*pi)) + logSurv
    }
    return(-2*sum(logl))
  }
  
  # Function calculating quasi-Monte carlo draws from multivariate normal
  # Quasi Monte-Carlo method for multivariate normal
  mnorm.prec = function(n,mu,prec)
  {
    # n: number of draws
    # mu: mean vector
    # prec: precision matrix
    
    p = ncol(prec)
    
    pivot = FALSE
    tol = -1
    
    # Upper triangular matrix
    L = .Internal(La_chol(prec, pivot, tol))
    
    z = t.default(sobol(n,dim = p,normal = T))
    
    sample = .Internal(backsolve(L, z, k = p, upper.tri = TRUE, transpose = FALSE)) + mu
    
    return(t.default(sample))
  }
  
  mnorm.CholPrec = function(n,mu,L)
  {
    # n: number of draws
    # mu: mean vector
    # prec: precision matrix
    
    p = ncol(L)
    
    z = matrix(rnorm(n*p),nr = p,nc = n)
    
    sample = .Internal(backsolve(L, z, k = p, upper.tri = TRUE, transpose = FALSE)) + mu
    
    return(t.default(sample))
  }
  
  pivot = FALSE
  tol = -1
  
  # Matrices/vectors needed to approximate the DIC criterion
  xTime = fitProp$xTime 
  xTimeLong = fitProp$xTimeLong 
  zTime = fitProp$zTime 
  zTimeLong = fitProp$zTimeLong 
  ylist = fitProp$ylist 
  xlist = fitProp$xlist 
  zlist = fitProp$zlist 
  zzlist = fitProp$zzlist 
  ni = fitProp$ni 
  G = fitProp$G 
  splBasisCause1 = fitProp$splBasisCause1
  splBasisCause1Long = fitProp$splBasisCause1Long
  splBasisCause2 = fitProp$splBasisCause2
  splBasisCause2Long = fitProp$splBasisCause2Long
  Lp = fitProp$Lp  
  q = fitProp$q 
  lpCause1 = fitProp$lpCause1
  lpCause2 = fitProp$lpCause2
  lbetasCause1 = fitProp$lbetasCause1
  lbetasCause2 = fitProp$lbetasCause2
  ids = fitProp$ids
  deltaCause1 = fitProp$deltaCause1
  deltaCause2 = fitProp$deltaCause2
  xCause1 = fitProp$xCause1
  xCause2 = fitProp$xCause2
  times = fitProp$times
  alpha1 = fitProp$alpha1
  alpha2 = fitProp$alpha2
  deltaCause1Obs = fitProp$deltaCause1Obs
  deltaCause2Obs = fitProp$deltaCause2Obs
  Rdouble = fitProp$Rdouble
  mubBar = t(apply(fitProp$draws.b,1, rowMeans))
  
  Sdraws = fitProp$Sdraws
  Ldraws = fitProp$Ldraws
  drawsD = fitProp$drawsD
  Mdraws = fitProp$Mdraws
  ndraw = nrow(Ldraws)
  
  ztTime = t(zTime)
  logl = rep(NA,G)
  
  # Deviance at the posterior mean
  L = chol(fitProp$summary$VarCovRE_mean)
  Lt = t(L)
  Lbeta = fitProp$summary$Longitudinal_Process[1:Lp,1]
  sigma2 = fitProp$summary$Longitudinal_Process[-c(1:Lp),1]
  paramSurv = colMeans(Sdraws)
  pi11 = mean(Mdraws[,1])
  pi22 = mean(Mdraws[,2])
  
  param = c(Lbeta,log(sigma2),vech(Lt),paramSurv,pi11,pi22)
  names(param) = NULL
  nMCbup = nMC
  nMC = nMC*15
  DevatPostMean = f(param)
  PostDev = rep(NA,ndraw)
  nMC = nMCbup
  
  for(j in seq(1,ndraw,by = thinDIC))
  {
    L = chol(drawsD[,,j])
    Lt = t(L)
    sigma2 = Ldraws[j,-c(1:Lp)]
    Lbeta = Ldraws[j,1:Lp]
    paramSurv = Sdraws[j,]
    pi11 = Mdraws[j,1]
    pi22 = Mdraws[j,2]
    
    param = c(Lbeta,log(sigma2),vech(Lt),paramSurv,pi11,pi22)
    
    PostDev[j] = f(param)
    print(j)
  }
  eff = mean(PostDev,na.rm = T) - DevatPostMean
  dic = mean(PostDev,na.rm = T) + eff
  out = list()
  out$eff = eff
  out$dic = dic
  out$DevatPostMean = DevatPostMean
  out$PostDev = PostDev
  return(out)
}




predictLatentStates_IncrCumInc = function(fitProp,newdata1,newdata2,tt = seq(0,10,by = 2),nMC = 1000,seed = 12,thin = 1,nburn = 200)
{
  # if (class(fitProp)!=T)
  # {
  #   
  # }
  
  set.seed(seed)
  # 7-15 Gauss-Kronrod quadradure points and weights
  gaussKronrod <-
    function (k = 15) {
      sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
              0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
              -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
              0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
      wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
      wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975, 
               0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
      if (k == 7) 
        list(sk = sk[1:7], wk = wk7)
      else
        list(sk = sk, wk = wk15)
    }
  
  gauleg <- function(n,a=-1,b=1) {# Gauss-Legendre: returns x,w so that
    #\int_a^b f(x) dx \doteq \sum w_i f(x_i)
    EPS <- 3.e-14
    m <- trunc((n+1)/2)
    xm <- 0.5*(b+a)
    xl <- 0.5*(b-a)
    x <- w <- rep(-1,n)
    for (i in 1:m) {
      z <- cos(pi*(i-.25)/(n+.5))
      z1 <- z+1
      while (abs(z-z1) > EPS) {
        p1 <- 1
        p2 <- 0
        for (j in 1:n) {# recursively evaluates pn(x)
          p3 <- p2
          p2 <- p1
          p1 <- ((2*j-1)*z*p2-(j-1)*p3)/j
        }
        pp <- n*(z*p1-p2)/(z*z-1)
        z1 <- z
        z <- z1-p1/pp #Newton iteration
      }
      x[i] <- xm-xl*z
      x[n+1-i] <- xm+xl*z
      w[i] <- 2*xl/((1-z*z)*pp*pp)
      w[n+1-i] <- w[i]
    }
    list(x=x,w=w)
  }
  
  useGauleg = fitProp$useGauleg
  Gauleg_points = fitProp$Gauleg_points
  
  if (useGauleg == F)
  {
    gk = gaussKronrod()
    pps = gk$sk
    wws = gk$wk
  } else {
    gk = gauleg(Gauleg_points)
    pps = gk$x
    wws = gk$w
  }
  
  # Transformation of time needed for gaussKronrod integration
  t1 = 0.5*tt
  t1Long = rep(t1,each = length(pps))
  ppsLong = rep(pps,length(tt))
  wwsLong = rep(wws,length(tt))
  ind = rep(1:length(tt),each = length(pps))
  
  c1 = fitProp$alpha1
  c2 = fitProp$alpha2
  
  # Indices of id's
  ids = list()
  for (i in 1:length(tt)){ids[[i]] = which(ind==i)}
  
  # Additional time points used to integrate the subditribution hazard
  timesForInt = t1Long*ppsLong + t1Long
  
  fitcox1 = fitProp$cox1
  fitcox2 = fitProp$cox2
  
  w1 = as.matrix(model.matrix(formula(fitcox1),newdata1)[,-1]) 
  w2 = as.matrix(model.matrix(formula(fitcox2),newdata2)[,-1])
  
  knotsCause1 = fitProp$knotsCause1
  knotsCause2 = fitProp$knotsCause2
  bknotsCause1 = fitProp$bknotsCause1
  bknotsCause2 = fitProp$bknotsCause2
  
  # Construct the spline basis for cause 1
  splBasisCause1 = bSpline(timesForInt,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T)
  colnames(splBasisCause1) = NULL
  
  # Construct the spline basis for cause 2
  splBasisCause2 = bSpline(timesForInt,knots = knotsCause2,Boundary.knots = bknotsCause2,intercept = T)
  colnames(splBasisCause2) = NULL
  
  lpCause1 = fitProp$lpCause1
  lpCause2 = fitProp$lpCause2
  lbetasCause1 = fitProp$lbetasCause1
  lbetasCause2 = fitProp$lbetasCause2
  lgammasCause1 = fitProp$lgammasCause1
  lgammasCause2 = fitProp$lgammasCause2
  
  # Get the data from the LME object
  fitlme = fitProp$fitlme
  dataLme = fitlme$data
  
  # Time variable in the lme model
  timeVar = fitProp$timeVar
  
  # The id variable in the lme model
  id = as.numeric(as.character(fitlme$groups[,1]))
  
  ######################################################################
  ### Create the design matrices of the fixed and the random effects ###
  ### at the time points of interest                                 ###
  ######################################################################
  dataLmeid = dataLme[!duplicated(id),][1:length(tt),]
  dataLmeid[,timeVar] = tt
  
  xTime = model.matrix(formula(fitlme),data = dataLmeid)
  colnames(xTime) = NULL
  rownames(xTime) = NULL
  
  zTime = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeid)
  colnames(zTime) = NULL
  rownames(zTime) = NULL
  
  # Dimension of the random effects
  q = ncol(zTime)
  
  # Get the predicted values of the marker at additional time points (needed to integrate the subdistribution hazard)
  
  ##################################################################################
  ### Create the design matrices of the fixed and the random effects             ###
  ### at additional time points (needed to integrate the subdistribution hazard) ###
  ##################################################################################
  dataLmeidLong = dataLme[!duplicated(id),][1:length(timesForInt),]
  dataLmeidLong[,timeVar] = timesForInt
  
  xTimeLong = model.matrix(formula(fitlme),data = dataLmeidLong)
  colnames(xTimeLong) = NULL
  rownames(xTimeLong) = NULL
  
  zTimeLong = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeidLong)
  colnames(zTimeLong) = NULL
  rownames(zTimeLong) = NULL
  
  # Currently it is unknown why it's here
  dd = dataLmeidLong[1:length(pps),]
  ddid = dd[1:2,]
  
  #br = rmvnorm(10000,sigma = Dtrue)
  
  # Matrices to integrate the subdistribution hazard function
  int1 = matrix(NA,nr = length(tt),nc = nMC)
  int2 = matrix(NA,nr = length(tt),nc = nMC)
  
  # Creates matrix of CIF estimates for cause 1
  # with rows corresponding to the time points and the columns to the random effects
  cif1 = function(br)
  {
    mtpredLong = c(xTimeLong %*% Lbeta) + tcrossprod(zTimeLong,br)
    
    if (lbetasCause1==1)
    {
      xbetasCause1 = mtpredLong*betasCause1[lbetasCause1]
    } else {
      xbetasCause1 = c(w1 %*% betasCause1[-lbetasCause1]) + mtpredLong*betasCause1[lbetasCause1]
    }
    
    #xbetasCause1 = c(w1 %*% betasCause1[-lbetasCause1]) + mtpredLong*betasCause1[lbetasCause1]
    uCause1 = c(splBasisCause1 %*% gammasCause1)
    linpredCause1 = xbetasCause1 + uCause1
    expLinpredCause1 = exp(linpredCause1)
    
    qq1 = (wwsLong*t1Long)*expLinpredCause1
    
    for (i in 1:length(tt))
    {
      int1[i,] = .colSums(qq1[ids[[i]],],length(pps),nMC,F)
    }
    1 - (1 + c1*int1)^(-1/c1)
  }
  
  # Creates CIF estimates for cause 1
  # for a specific time point and random effects values
  cif1eval = function(time,br)
  {
    # Tranform times for integration
    t1 = 0.5*time
    t1Long = rep(t1,each = length(pps))
    timesForInt = t1Long*pps + t1Long
    #dd[,timeVar] = timesForInt
    
    temp = cbind(1,lspline(timesForInt,knots = c(1,5)))
    
    # The design matrices at the additional points 
    xTimeLong = temp
    
    zTimeLong = temp[,-4]
    
    # Construct the spline basis for cause 1
    splBasisCause1 = bSpline(timesForInt,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T)
    colnames(splBasisCause1) = NULL
    
    # The true marker value at all points needed for integration
    mtpredLong = c(xTimeLong %*% Lbeta) + zTimeLong %*% br
    
    if (lbetasCause1==1)
    {
      xbetasCause1 = mtpredLong*betasCause1[lbetasCause1]
    } else {
      xbetasCause1 = c(w1 %*% betasCause1[-lbetasCause1]) + mtpredLong*betasCause1[lbetasCause1]
    }
    
    #xbetasCause1 = c(w1 %*% betasCause1[-lbetasCause1]) + mtpredLong*betasCause1[lbetasCause1]
    uCause1 = c(splBasisCause1 %*% gammasCause1)
    linpredCause1 = xbetasCause1 + uCause1
    expLinpredCause1 = exp(linpredCause1)
    
    qq1 = (wws*t1Long)*expLinpredCause1
    
    1 - (1 + c1*sum(qq1))^(-1/c1)
  }
  
  
  # Creates vector of CIF estimates for cause 1
  # for a specific time point for all random effects values
  cif1evalVec = function(k,br)
  {
    # The true marker value at all points needed for integration
    mtpredLong = c(xTimeList[[k]] %*% Lbeta) + tcrossprod(zTimeList[[k]],br)
    
    if (lbetasCause1==1)
    {
      xbetasCause1 = mtpredLong*betasCause1[lbetasCause1]
    } else {
      xbetasCause1 = c(w1 %*% betasCause1[-lbetasCause1]) + mtpredLong*betasCause1[lbetasCause1]
    }
    
    #xbetasCause1 = c(w1 %*% betasCause1[-lbetasCause1]) + mtpredLong*betasCause1[lbetasCause1]
    uCause1 = c(splBasisCause1List[[k]] %*% gammasCause1)
    linpredCause1 = xbetasCause1 + uCause1
    expLinpredCause1 = exp(linpredCause1)
    
    qq1 = (wws*tt[k]/2)*expLinpredCause1
    
    1 - (1 + c1*colSums(qq1))^(-1/c1)
  }
  
  # Creates matrix of CIF estimates for cause 2
  # with rows corresponding to the time points and the columns to the random effects
  cif2 = function(br)
  {
    mtpredLong = c(xTimeLong %*% Lbeta) + tcrossprod(zTimeLong,br)
    
    if (lbetasCause2==1)
    {
      xbetasCause2 = mtpredLong*betasCause2[lbetasCause2]
    } else {
      xbetasCause2 = c(w2 %*% betasCause2[-lbetasCause2]) + mtpredLong*betasCause2[lbetasCause2]
    }
    
    #xbetasCause2 = c(w2 %*% betasCause2[-lbetasCause2]) + mtpredLong*betasCause2[lbetasCause2]
    uCause2 = c(splBasisCause2 %*% gammasCause2)
    linpredCause2 = xbetasCause2 + uCause2
    expLinpredCause2 = exp(linpredCause2)
    
    qq2 = (wwsLong*t1Long)*expLinpredCause2
    
    for (i in 1:length(tt))
    {
      int2[i,] = .colSums(qq2[ids[[i]],],length(pps),nMC,F)
    }
    1 - (1 + c2*int2)^(-1/c2)
  }
  
  # Creates CIF estimates for cause 2
  # for a specific time point and random effects values
  cif2eval = function(time,br)
  {
    # Tranform times for integration
    t1 = 0.5*time
    t1Long = rep(t1,each = length(pps))
    timesForInt = t1Long*pps + t1Long
    dd[,timeVar] = timesForInt
    
    temp = cbind(1,lspline(timesForInt,knots = c(1,5)))
    
    # The design matrices at the additional points 
    xTimeLong = temp
    
    zTimeLong = temp[,-4]
    
    
    # Construct the spline basis for cause 2
    splBasisCause2 = bSpline(timesForInt,knots = knotsCause2,Boundary.knots = bknotsCause2,intercept = T)
    colnames(splBasisCause2) = NULL
    
    # The true marker value at all points needed for integration
    mtpredLong = c(xTimeLong %*% Lbeta) + zTimeLong %*% br
    
    if (lbetasCause2==1)
    {
      xbetasCause2 = mtpredLong*betasCause2[lbetasCause2]
    } else {
      xbetasCause2 = c(w2 %*% betasCause2[-lbetasCause2]) + mtpredLong*betasCause2[lbetasCause2]
    }
    
    #xbetasCause2 = c(w2 %*% betasCause2[-lbetasCause2]) + mtpredLong*betasCause2[lbetasCause2]
    uCause2 = c(splBasisCause2 %*% gammasCause2)
    linpredCause2 = xbetasCause2 + uCause2
    expLinpredCause2 = exp(linpredCause2)
    
    qq2 = (wws*t1Long)*expLinpredCause2
    
    1 - (1 + c2*sum(qq2))^(-1/c2)
  }
  
  # Creates vector of CIF estimates for cause 2
  # for a specific time point for all random effects values
  cif2evalVec = function(k,br)
  {
    # The true marker value at all points needed for integration
    mtpredLong = c(xTimeList[[k]] %*% Lbeta) + tcrossprod(zTimeList[[k]],br)
    
    if (lbetasCause2==1)
    {
      xbetasCause2 = mtpredLong*betasCause2[lbetasCause2]
    } else {
      xbetasCause2 = c(w2 %*% betasCause2[-lbetasCause2]) + mtpredLong*betasCause2[lbetasCause2]
    }
    
    #xbetasCause2 = c(w2 %*% betasCause2[-lbetasCause2]) + mtpredLong*betasCause2[lbetasCause2]
    uCause2 = c(splBasisCause2List[[k]] %*% gammasCause2)
    linpredCause2 = xbetasCause2 + uCause2
    expLinpredCause2 = exp(linpredCause2)
    
    qq2 = (wws*tt[k]/2)*expLinpredCause2
    
    1 - (1 + c2*colSums(qq2))^(-1/c2)
  }
  
  ff = function(time)
  {
    cif1eval(time,br[k,]) + cif2eval(time,br[k,]) - 1
  }
  
  # ffder = function(time)
  # {
  #   cif1evalDer(time,br[k,]) + cif2evalDer(time,br[k,]) 
  # }
  
  ndraw = nrow(fitProp$Sdraws)
  
  Sdraws = fitProp$Sdraws[1:(ndraw/thin)*thin,]
  Ldraws = fitProp$Ldraws[1:(ndraw/thin)*thin,]
  drawsD = fitProp$drawsD[,,1:(ndraw/thin)*thin]
  Sdraws1 = Sdraws[,1:lpCause1]
  Sdraws2 = Sdraws[,-c(1:lpCause1)]
  
  ndraw = nrow(Sdraws)
  
  ################################################################################################################
  ### Posterior draws for the population CIFs conditional on the baseline state                                ###
  ### The rows correspond to posterior draws, the columns to the times at which the marginal CIFs are computed ###
  ### and the third dimension of the array to the baseline state                                               ###
  ################################################################################################################
  ss1 = array(NA,dim = c(ndraw,length(tt),7) )
  ss2 = array(NA,dim = c(ndraw,length(tt),7) )
  
  ss1Overall = matrix(NA,nr = ndraw,nc = length(tt))
  ss2Overall = matrix(NA,nr = ndraw,nc = length(tt))
  
  states = c(-Inf,sqrt(c(50,100,200,250,350,500)))
  states2 = c(sqrt(c(50,100,200,250,350,500)),Inf)
  
  # Latent state probabilities for each time point
  state1Overall = matrix(nr = ndraw,nc = length(tt))
  state2Overall = matrix(nr = ndraw,nc = length(tt))
  state3Overall = matrix(nr = ndraw,nc = length(tt))
  state4Overall = matrix(nr = ndraw,nc = length(tt))
  state5Overall = matrix(nr = ndraw,nc = length(tt))
  state6Overall = matrix(nr = ndraw,nc = length(tt))
  state7Overall = matrix(nr = ndraw,nc = length(tt))
  
  ############################################################################################
  ### Posterior draws for the transition probabilities to states from given baseline state ###
  ############################################################################################
  # The rows correspond to posterior draws
  # the columns to time points
  # the third dimension of the array to the state to which transition is made
  # the component of the list corresponds to the baseline state
  
  stateList = list()
  for(i in 1:7)
  {
    stateList[[i]] = array(NA,dim = c(ndraw,length(tt),7) )
  }
  
  # Marginal probabilities for each state
  margStateProbs = matrix(nr = 7,nc = length(tt))
  
  # Marginal transition probabilities from given baseline states
  transStateProbs = list()
  
  for (i in 1:7)
  {
    transStateProbs[[i]] = matrix(NA,nr = 7,ncol = length(tt)-1)
  }
  
  # Design matrices for each time point
  xTimeList = list()
  zTimeList = list()
  splBasisCause1List = list()
  splBasisCause2List = list()
  
  for (k in 1:length(tt))
  {
    xTimeList[[k]] = xTimeLong[ids[[k]],]
    zTimeList[[k]] = zTimeLong[ids[[k]],]
    splBasisCause1List[[k]] = splBasisCause1[ids[[k]],]
    splBasisCause2List[[k]] = splBasisCause2[ids[[k]],]
  }
  
  # Starting values for the random effects meeting the constraints (based on the posterior means of the parameters)
  # Needed in Gibbs sampling/ HMC sampling
  startingValues <- expand.grid(state_start = paste(round(states^2),round(states2^2)),state_stop = paste(round(states^2),round(states2^2)),
                                k = 2:length(tt))
  startingValues$l1 = NA
  startingValues$l1[startingValues$state_start=="Inf 50"] = 1
  startingValues$l1[startingValues$state_start=="50 100"] = 2
  startingValues$l1[startingValues$state_start=="100 200"] = 3
  startingValues$l1[startingValues$state_start=="200 250"] = 4
  startingValues$l1[startingValues$state_start=="250 350"] = 5
  startingValues$l1[startingValues$state_start=="350 500"] = 6
  startingValues$l1[startingValues$state_start=="500 Inf"] = 7
  
  startingValues$l2 = NA
  startingValues$l2[startingValues$state_stop=="Inf 50"] = 1
  startingValues$l2[startingValues$state_stop=="50 100"] = 2
  startingValues$l2[startingValues$state_stop=="100 200"] = 3
  startingValues$l2[startingValues$state_stop=="200 250"] = 4
  startingValues$l2[startingValues$state_stop=="250 350"] = 5
  startingValues$l2[startingValues$state_stop=="350 500"] = 6
  startingValues$l2[startingValues$state_stop=="500 Inf"] = 7
  
  br = rmvnorm(10^4,sigma = apply(drawsD,1:2,mean)  )
  mtpred = c(xTime %*% colMeans(Ldraws)[1:ncol(xTime)] ) + tcrossprod(zTime,br)
  
  startingValues = startingValues[order(startingValues$k,startingValues$l1,startingValues$l2),]
  startingValues = cbind(startingValues,matrix(NA,nr = nrow(startingValues),nc = q))
  names(startingValues)[-c(1:5)] = paste0("startb",1:q)
  startingValues$it = 1:nrow(startingValues)
  
  for (i in startingValues$it)
  {
    l1 = startingValues$l1[i]
    l2 = startingValues$l2[i]
    k = startingValues$k[i]
    
    jj = which(mtpred[1,] > states[l1] & mtpred[1,] < states2[l1] & mtpred[k,] > states[l2] & mtpred[k,] < states2[l2])[1]
    startingValues[i,paste0("startb",1:q)] = br[jj,]
  }
  
  # Used to obtain feasible starting values for the HMC algorithm
  omegastar = mean(1/Ldraws[,ncol(Ldraws)])*10^6
  ss = sqrt(c(25,75,150,225,300,425,600))
  
  for (j in 1:nrow(Sdraws))
  {
    # The jth draw from from the posterior distribution of the parameters
    betasCause1 = Sdraws1[j,1:lbetasCause1]
    gammasCause1 = Sdraws1[j,(lbetasCause1+1):lpCause1]
    Lbeta = Ldraws[j,-ncol(Ldraws)]
    D = drawsD[,,j]
    Dinv = solve(D)
    
    betasCause2 = Sdraws2[j,1:lbetasCause2]
    gammasCause2 = Sdraws2[j,(lbetasCause2+1):lpCause2]
    
    # Marginal mean and variance of the `true' marker values at time points tt
    mu = c(xTime %*% Lbeta)
    vars = rep(NA,length(tt))
    
    for (k in 1:length(tt))
    {
      vars[k] = t(zTime[k,]) %*% D %*% zTime[k,]
    }
    
    # Marginal probabilities of being in a specific state at time points tt.
    for (k in 1:length(tt))
    {
      margStateProbs[,k] = pnorm(states2,mean = mu[k],sd = sqrt(vars[k])) - pnorm(states,mean = mu[k],sd = sqrt(vars[k]))
    }
    
    # Marginal transition probabilities from given baseline states
    for (i in 1:7)
    {
      for (jj in 2:length(tt))
      {
        Sigma <- rbind( c(vars[jj],zTime[1,]%*% D %*% zTime[jj,] ) ,c(zTime[jj,]%*% D %*% zTime[1,],vars[1]) )
        mus <- c(mu[jj],mu[1])
        
        for (l in 1:7)
        {
          transStateProbs[[i]][l,jj-1] <- pmvnorm(lower = c(states[l],states[i]),upper = c(states2[l],states2[i]),mean = mus,sigma = Sigma)
        }
      }
    }
    
    # Starting values for the Hamiltonian Monte Carlo technique
    # ok = F
    # kk = 0 
    # while(ok == F)
    # {
    #   kk = kk + 1
    #   
    #   br = rmvnorm(10^4,sigma = 10*D)
    #   mtpred = c(xTime %*% Lbeta ) + tcrossprod(zTime,br)
    #   
    #   for (i in startingValues$it)
    #   {
    #     l1 = startingValues$l1[i]
    #     l2 = startingValues$l2[i]
    #     k = startingValues$k[i]
    #     
    #     jj = which(mtpred[1,] > states[l1] & mtpred[1,] < states2[l1] & mtpred[k,] > states[l2] & mtpred[k,] < states2[l2])[1]
    #     
    #     startingValues$startb1[i] = br[jj,1]
    #     startingValues$startb2[i] = br[jj,2]
    #     startingValues$startb3[i] = br[jj,3]
    #   }
    #   
    #   ok <- !any(is.na(startingValues$startb1))
    # }
    
    for (i in startingValues$it)
    {
      l1 = startingValues$l1[i]
      l2 = startingValues$l2[i]
      k = startingValues$k[i]
      
      zi = zTime[c(1,k),]
      
      bi = omegastar*solve(Dinv + crossprod(zi)*omegastar)%*%crossprod(zi,ss[c(l1,l2)]-mu[c(1,k)]) 
      
      startingValues[i,paste0("startb",1:q)] = bi
      
    }
    
    #####################################################################
    ### Draws from N(0,D), contrainted such that m_{i}(0) < \sqrt{50},### 
    ### using Hamiltonian Monte Carlo                                 ###
    #####################################################################
    br <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = c(sqrt(40)-mu[1],rep(0,q-1)),
               f = - t(as.matrix(zTime[1,])), g = sqrt(50) - mu[1], burn.in = nburn)
    a1 = cif1(br)
    a2 = cif2(br)
    
    # Overall survival function for group 1
    st = 1 - a1 - a2
    
    # Samples from the prior that resulted in negative overall survival function
    indNegSt = which(st[nrow(st),] < 0)
    for (k in indNegSt)
    {
      indTime = which(st[,k]<0)[1]
      
      # Calculate \tau_{i}, the upper bound for the overall survival time
      fit = uniroot(f = ff,lower = tt[indTime-1],upper = tt[indTime])
      taui = fit$root
      
      # taui = 0.5*(tt[indTime-1] + tt[indTime])
      # for (i in 1:10)
      # {
      #   taui = taui - ff(taui)/ffder(taui)
      #   print(taui)
      # }
      
      a1[tt>taui,k] = cif1eval(taui,br[k,])
      a2[tt>taui,k] = cif2eval(taui,br[k,])
      st[tt>taui,k] = 0
    }
    
    ss1[j,,1] = rowMeans(a1)
    ss2[j,,1] = rowMeans(a2)
    
    #############################################################################################
    ### Draws from N(0,D), contrainted such that m_{i}(0) < \sqrt{100} & m_{i}(0) > \sqrt{50},### 
    ### using Hamiltonian Monte Carlo                                                         ###
    #############################################################################################
    br <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = c(sqrt(70)-mu[1],rep(0,q-1)),
               f = rbind(zTime[1,],-zTime[1,]), g = c(mu[1]-sqrt(50),sqrt(100)-mu[1]), burn.in = nburn)
    
    a1 = cif1(br)
    a2 = cif2(br)
    
    # Overall survival function for group 1
    st = 1 - a1 - a2
    
    # Samples from the prior that resulted in negative overall survival function
    indNegSt = which(st[nrow(st),] < 0)
    for (k in indNegSt)
    {
      indTime = which(st[,k]<0)[1]
      
      # Calculate \tau_{i}, the upper bound for the overall survival time
      fit = uniroot(f = ff,lower = tt[indTime-1],upper = tt[indTime])
      taui = fit$root
      
      # taui = 0.5*(tt[indTime-1] + tt[indTime])
      # for (i in 1:10)
      # {
      #   taui = taui - ff(taui)/ffder(taui)
      #   print(taui)
      # }
      
      a1[tt>taui,k] = cif1eval(taui,br[k,])
      a2[tt>taui,k] = cif2eval(taui,br[k,])
      st[tt>taui,k] = 0
    }
    
    ss1[j,,2] = rowMeans(a1)
    ss2[j,,2] = rowMeans(a2)
    
    #############################################################################################
    ### Draws from N(0,D), contrainted such that m_{i}(0) < \sqrt{200} & m_{i}(0) > \sqrt{100}### 
    ### using Hamiltonian Monte Carlo                                                         ###
    #############################################################################################
    br <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = c(sqrt(150)-mu[1],rep(0,q-1)),
               f = rbind(zTime[1,],-zTime[1,]), g = c(mu[1]-sqrt(100),sqrt(200)-mu[1]), burn.in = nburn)
    
    a1 = cif1(br)
    a2 = cif2(br)
    
    # Overall survival function for group 1
    st = 1 - a1 - a2
    
    # Samples from the prior that resulted in negative overall survival function
    indNegSt = which(st[nrow(st),] < 0)
    for (k in indNegSt)
    {
      indTime = which(st[,k]<0)[1]
      
      # Calculate \tau_{i}, the upper bound for the overall survival time
      fit = uniroot(f = ff,lower = tt[indTime-1],upper = tt[indTime])
      taui = fit$root
      
      # taui = 0.5*(tt[indTime-1] + tt[indTime])
      # for (i in 1:10)
      # {
      #   taui = taui - ff(taui)/ffder(taui)
      #   print(taui)
      # }
      
      a1[tt>taui,k] = cif1eval(taui,br[k,])
      a2[tt>taui,k] = cif2eval(taui,br[k,])
      st[tt>taui,k] = 0
    }
    
    ss1[j,,3] = rowMeans(a1)
    ss2[j,,3] = rowMeans(a2)
    
    #############################################################################################
    ### Draws from N(0,D), contrainted such that m_{i}(0) < \sqrt{250} & m_{i}(0) > \sqrt{200}### 
    ### using Hamiltonian Monte Carlo                                                         ###
    #############################################################################################
    br <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = c(sqrt(220)-mu[1],rep(0,q-1)),
               f = rbind(zTime[1,],-zTime[1,]), g = c(mu[1]-sqrt(200),sqrt(250)-mu[1]), burn.in = nburn)
    
    a1 = cif1(br)
    a2 = cif2(br)
    
    # Overall survival function for group 1
    st = 1 - a1 - a2
    
    # Samples from the prior that resulted in negative overall survival function
    indNegSt = which(st[nrow(st),] < 0)
    for (k in indNegSt)
    {
      indTime = which(st[,k]<0)[1]
      
      # Calculate \tau_{i}, the upper bound for the overall survival time
      fit = uniroot(f = ff,lower = tt[indTime-1],upper = tt[indTime])
      taui = fit$root
      
      # taui = 0.5*(tt[indTime-1] + tt[indTime])
      # for (i in 1:10)
      # {
      #   taui = taui - ff(taui)/ffder(taui)
      #   print(taui)
      # }
      
      a1[tt>taui,k] = cif1eval(taui,br[k,])
      a2[tt>taui,k] = cif2eval(taui,br[k,])
      st[tt>taui,k] = 0
    }
    
    ss1[j,,4] = rowMeans(a1)
    ss2[j,,4] = rowMeans(a2)
    
    #############################################################################################
    ### Draws from N(0,D), contrainted such that m_{i}(0) < \sqrt{350} & m_{i}(0) > \sqrt{250}### 
    ### using Hamiltonian Monte Carlo                                                         ###
    #############################################################################################
    br <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = c(sqrt(300)-mu[1],rep(0,q-1)),
               f = rbind(zTime[1,],-zTime[1,]), g = c(mu[1]-sqrt(250),sqrt(350)-mu[1]), burn.in = nburn)
    
    a1 = cif1(br)
    a2 = cif2(br)
    
    # Overall survival function for group 1
    st = 1 - a1 - a2
    
    # Samples from the prior that resulted in negative overall survival function
    indNegSt = which(st[nrow(st),] < 0)
    for (k in indNegSt)
    {
      indTime = which(st[,k]<0)[1]
      
      # Calculate \tau_{i}, the upper bound for the overall survival time
      fit = uniroot(f = ff,lower = tt[indTime-1],upper = tt[indTime])
      taui = fit$root
      
      # taui = 0.5*(tt[indTime-1] + tt[indTime])
      # for (i in 1:10)
      # {
      #   taui = taui - ff(taui)/ffder(taui)
      #   print(taui)
      # }
      
      a1[tt>taui,k] = cif1eval(taui,br[k,])
      a2[tt>taui,k] = cif2eval(taui,br[k,])
      st[tt>taui,k] = 0
    }
    
    ss1[j,,5] = rowMeans(a1)
    ss2[j,,5] = rowMeans(a2)
    
    #############################################################################################
    ### Draws from N(0,D), contrainted such that m_{i}(0) < \sqrt{500} & m_{i}(0) > \sqrt{350}### 
    ### using Hamiltonian Monte Carlo                                                         ###
    #############################################################################################
    br <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = c(sqrt(400)-mu[1],rep(0,q-1)),
               f = rbind(zTime[1,],-zTime[1,]), g = c(mu[1]-sqrt(350),sqrt(500)-mu[1]), burn.in = nburn)
    
    a1 = cif1(br)
    a2 = cif2(br)
    
    # Overall survival function for group 1
    st = 1 - a1 - a2
    
    # Samples from the prior that resulted in negative overall survival function
    indNegSt = which(st[nrow(st),] < 0)
    for (k in indNegSt)
    {
      indTime = which(st[,k]<0)[1]
      
      # Calculate \tau_{i}, the upper bound for the overall survival time
      fit = uniroot(f = ff,lower = tt[indTime-1],upper = tt[indTime])
      taui = fit$root
      
      # taui = 0.5*(tt[indTime-1] + tt[indTime])
      # for (i in 1:10)
      # {
      #   taui = taui - ff(taui)/ffder(taui)
      #   print(taui)
      # }
      
      a1[tt>taui,k] = cif1eval(taui,br[k,])
      a2[tt>taui,k] = cif2eval(taui,br[k,])
      st[tt>taui,k] = 0
    }
    
    ss1[j,,6] = rowMeans(a1)
    ss2[j,,6] = rowMeans(a2)
    
    ######################################################################
    ### Draws from N(0,D), contrainted such that m_{i}(0) > \sqrt{500} ### 
    ### using Hamiltonian Monte Carlo                                  ###
    ######################################################################
    br <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = c(sqrt(550)-mu[1],rep(0,q-1)),
               f = rbind(zTime[1,]), g = c(mu[1]-sqrt(500)), burn.in = nburn)
    
    a1 = cif1(br)
    a2 = cif2(br)
    
    # Overall survival function for group 1
    st = 1 - a1 - a2
    
    # Samples from the prior that resulted in negative overall survival function
    indNegSt = which(st[nrow(st),] < 0)
    for (k in indNegSt)
    {
      indTime = which(st[,k]<0)[1]
      
      # Calculate \tau_{i}, the upper bound for the overall survival time
      fit = uniroot(f = ff,lower = tt[indTime-1],upper = tt[indTime])
      taui = fit$root
      
      # taui = 0.5*(tt[indTime-1] + tt[indTime])
      # for (i in 1:10)
      # {
      #   taui = taui - ff(taui)/ffder(taui)
      #   print(taui)
      # }
      
      a1[tt>taui,k] = cif1eval(taui,br[k,])
      a2[tt>taui,k] = cif2eval(taui,br[k,])
      st[tt>taui,k] = 0
    }
    
    ss1[j,,7] = rowMeans(a1)
    ss2[j,,7] = rowMeans(a2)
    
    # mtpred = c(xTime %*% Lbeta) + tcrossprod(zTime,br)
    # drwss1 = mtpred < sqrt(50)
    # drwss2 = mtpred > sqrt(50) &  mtpred < sqrt(100)
    # drwss3 = mtpred > sqrt(100) &  mtpred < sqrt(200)
    # drwss4 = mtpred > sqrt(200) &  mtpred < sqrt(250)
    # drwss5 = mtpred > sqrt(250) &  mtpred < sqrt(350)
    # drwss6 = mtpred > sqrt(350) &  mtpred < sqrt(500)
    # drwss7 = mtpred > sqrt(500)
    # 
    # state1[j,] = rowMeans( drwss1*st )
    # state2[j,] = rowMeans( drwss2*st )
    # state3[j,] = rowMeans( drwss3*st )
    # state4[j,] = rowMeans( drwss4*st )
    # state5[j,] = rowMeans( drwss5*st )
    # state6[j,] = rowMeans( drwss6*st )
    # state7[j,] = rowMeans( drwss7*st )
    
    for (k in 2:length(tt))
    {
      
      for (l in 1:7)
      {
        
        ####################################################################
        ### m_{i}(0) < sqrt(50) & m_{i}(t) \in S_{l} (in the l-th state) ###
        ####################################################################
        startPoint = as.numeric(startingValues[startingValues$l1==1 & startingValues$l2 == l & startingValues$k == k,
                                               c("startb1","startb2","startb3")])
        
        # If the marginal probability of the state is too low, use Gibbs sampling
        # otherwise we use HMC sampling
        
        if (transStateProbs[[1]][l,k-1] < 1e-6)
        {
          # Gibbs sampling: we carry out a simple approach checking convergence of the Gibbs sampling
          ok = F
          jj = 0 
          while(ok == F)
          {
            jj = jj + 1
            if (l == 1)
            {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],zTime[k,]),
                                             offset = c(states2[1] - mu[1],states2[l] - mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else if (l==7) {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[k,]),
                                             offset = c(states2[1] - mu[1], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],zTime[k,],-zTime[k,]),
                                             offset = c(states2[1] - mu[1],states2[l] - mu[k], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            }
            
            #diag(cov(samp))
            ok = !all(diag(cov(samp)) < 1e-4)
          }
          
        } else{
          
          # Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          if (l == 1)
          {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(-zTime[1,],-zTime[k,]), 
                         g = c( states2[1]-mu[1], states2[l]-mu[k]), burn.in = nburn)
          } else if (l==7) {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(-zTime[1,],zTime[k,]), 
                         g = c( states2[1]-mu[1], mu[k]-states[l] ), burn.in = nburn)
          } else {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(-zTime[1,],-zTime[k,],zTime[k,]), 
                         g = c( states2[1]-mu[1], states2[l]-mu[k], mu[k]-states[l] ), burn.in = nburn)
          }
          
          
        }
        
        a1 = cif1evalVec(k,samp)
        a2 = cif2evalVec(k,samp)
        
        # Overall survival function for group 1
        st = 1 - a1 - a2
        st[st<0] = 0
        
        stateList[[1]][j,k,l] <- mean(st)*transStateProbs[[1]][l,k-1]/margStateProbs[1,1]
        
        #############################################################################################
        ### m_{i}(0) < \sqrt{100} & m_{i}(0) > \sqrt{50} & m_{i}(t) \in S_{l} (in the l-th state) ###
        #############################################################################################
        startPoint = as.numeric(startingValues[startingValues$l1==2 & startingValues$l2 == l & startingValues$k == k,
                                               c("startb1","startb2","startb3")])
        
        # If the marginal probability of the state is too low, use Gibbs sampling
        # otherwise we use HMC sampling
        
        if (transStateProbs[[2]][l,k-1] < 1e-6)
        {
          # Gibbs sampling: we carry out a simple approach checking convergence of the Gibbs sampling
          ok = F
          jj = 0 
          while(ok == F)
          {
            jj = jj + 1
            if (l == 1)
            {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,]),
                                             offset = c( states2[2] - mu[1],- states[2] + mu[1],states2[l] - mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else if (l == 7){
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],-zTime[k,]),
                                             offset = c( states2[2] - mu[1],- states[2] + mu[1], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,],-zTime[k,]),
                                             offset = c( states2[2] - mu[1],- states[2] + mu[1],states2[l] - mu[k], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            }
            
            #diag(cov(samp))
            ok = !all(diag(cov(samp)) < 1e-4)
          }
          
        } else{
          
          # Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          if (l == 1)
          {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,]), 
                         g = c(mu[1]-states[2],states2[2]-mu[1], states2[l]-mu[k]), burn.in = nburn)
          } else if (l==7) {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],zTime[k,]), 
                         g = c(mu[1]-states[2],states2[2]-mu[1], mu[k]-states[l] ), burn.in = nburn)
          } else {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,],zTime[k,]), 
                         g = c(mu[1]-states[2],states2[2]-mu[1], states2[l]-mu[k], mu[k]-states[l] ), burn.in = nburn)
          }
          
        }
        
        a1 = cif1evalVec(k,samp)
        a2 = cif2evalVec(k,samp)
        
        # Overall survival function for group 1
        st = 1 - a1 - a2
        st[st<0] = 0
        stateList[[2]][j,k,l] = mean(st)*transStateProbs[[2]][l,k-1]/margStateProbs[2,1]
        
        ##############################################################################################
        ### m_{i}(0) < \sqrt{200} & m_{i}(0) > \sqrt{100} & m_{i}(t) \in S_{l} (in the l-th state) ###
        ##############################################################################################
        startPoint = as.numeric(startingValues[startingValues$l1==3 & startingValues$l2 == l & startingValues$k == k,
                                               c("startb1","startb2","startb3")])
        
        
        # If the marginal probability of the state is too low, use Gibbs sampling
        # otherwise we use HMC sampling
        
        if (transStateProbs[[3]][l,k-1] < 1e-6)
        {
          # Gibbs sampling: we carry out a simple approach checking convergence of the Gibbs sampling
          ok = F
          jj = 0 
          while(ok == F)
          {
            jj = jj + 1
            if (l == 1)
            {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,]),
                                             offset = c( states2[3] - mu[1],- states[3] + mu[1],states2[l] - mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else if (l == 7){
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],-zTime[k,]),
                                             offset = c( states2[3] - mu[1],- states[3] + mu[1], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,],-zTime[k,]),
                                             offset = c( states2[3] - mu[1],- states[3] + mu[1],states2[l] - mu[k], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            }
            
            #diag(cov(samp))
            ok = !all(diag(cov(samp)) < 1e-4)
          }
          
        } else{
          
          # Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          if (l == 1)
          {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,]), 
                         g = c(mu[1]-states[3],states2[3]-mu[1], states2[l]-mu[k]), burn.in = nburn)
          } else if (l==7) {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],zTime[k,]), 
                         g = c(mu[1]-states[3],states2[3]-mu[1], mu[k]-states[l] ), burn.in = nburn)
          } else {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,],zTime[k,]), 
                         g = c(mu[1]-states[3],states2[3]-mu[1], states2[l]-mu[k], mu[k]-states[l] ), burn.in = nburn)
          }
          
        }
        
        a1 = cif1evalVec(k,samp)
        a2 = cif2evalVec(k,samp)
        
        # Overall survival function for group 1
        st = 1 - a1 - a2
        st[st<0] = 0
        stateList[[3]][j,k,l] = mean(st)*transStateProbs[[3]][l,k-1]/margStateProbs[3,1]
        
        ##############################################################################################
        ### m_{i}(0) < \sqrt{250} & m_{i}(0) > \sqrt{200} & m_{i}(t) \in S_{l} (in the l-th state) ###
        ##############################################################################################
        startPoint = as.numeric(startingValues[startingValues$l1==4 & startingValues$l2 == l & startingValues$k == k,
                                               c("startb1","startb2","startb3")])
        
        # If the marginal probability of the state is too low, use Gibbs sampling
        # otherwise we use HMC sampling
        
        if (transStateProbs[[4]][l,k-1] < 1e-6)
        {
          # Gibbs sampling: we carry out a simple approach checking convergence of the Gibbs sampling
          ok = F
          jj = 0 
          while(ok == F)
          {
            jj = jj + 1
            if (l == 1)
            {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,]),
                                             offset = c( states2[4] - mu[1],- states[4] + mu[1],states2[l] - mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else if (l == 7){
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],-zTime[k,]),
                                             offset = c( states2[4] - mu[1],- states[4] + mu[1], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,],-zTime[k,]),
                                             offset = c( states2[4] - mu[1],- states[4] + mu[1],states2[l] - mu[k], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            }
            
            #diag(cov(samp))
            ok = !all(diag(cov(samp)) < 1e-4)
          }
          
        } else{
          
          # Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          if (l == 1)
          {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,]), 
                         g = c(mu[1]-states[4],states2[4]-mu[1], states2[l]-mu[k]), burn.in = nburn)
          } else if (l==7) {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],zTime[k,]), 
                         g = c(mu[1]-states[4],states2[4]-mu[1], mu[k]-states[l] ), burn.in = nburn)
          } else {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,],zTime[k,]), 
                         g = c(mu[1]-states[4],states2[4]-mu[1], states2[l]-mu[k], mu[k]-states[l] ), burn.in = nburn)
          }
          
        }
        
        a1 = cif1evalVec(k,samp)
        a2 = cif2evalVec(k,samp)
        
        # Overall survival function for group 1
        st = 1 - a1 - a2
        st[st<0] = 0
        stateList[[4]][j,k,l] = mean(st)*transStateProbs[[4]][l,k-1]/margStateProbs[4,1]
        
        ##############################################################################################
        ### m_{i}(0) < \sqrt{350} & m_{i}(0) > \sqrt{250} & m_{i}(t) \in S_{l} (in the l-th state) ###
        ##############################################################################################
        startPoint = as.numeric(startingValues[startingValues$l1==5 & startingValues$l2 == l & startingValues$k == k,
                                               c("startb1","startb2","startb3")])
        
        # If the marginal probability of the state is too low, use Gibbs sampling
        # otherwise we use HMC sampling
        
        if (transStateProbs[[5]][l,k-1] < 1e-6)
        {
          # Gibbs sampling: we carry out a simple approach checking convergence of the Gibbs sampling
          ok = F
          jj = 0 
          while(ok == F)
          {
            jj = jj + 1
            if (l == 1)
            {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,]),
                                             offset = c( states2[5] - mu[1],- states[5] + mu[1],states2[l] - mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else if (l == 7){
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],-zTime[k,]),
                                             offset = c( states2[5] - mu[1],- states[5] + mu[1], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,],-zTime[k,]),
                                             offset = c( states2[5] - mu[1],- states[5] + mu[1],states2[l] - mu[k], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            }
            
            #diag(cov(samp))
            ok = !all(diag(cov(samp)) < 1e-4)
          }
          
        } else{
          
          # Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          if (l == 1)
          {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,]), 
                         g = c(mu[1]-states[5],states2[5]-mu[1], states2[l]-mu[k]), burn.in = nburn)
          } else if (l==7) {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],zTime[k,]), 
                         g = c(mu[1]-states[5],states2[5]-mu[1], mu[k]-states[l] ), burn.in = nburn)
          } else {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,],zTime[k,]), 
                         g = c(mu[1]-states[5],states2[5]-mu[1], states2[l]-mu[k], mu[k]-states[l] ), burn.in = nburn)
          }
          
        }
        
        a1 = cif1evalVec(k,samp)
        a2 = cif2evalVec(k,samp)
        
        # Overall survival function for group 1
        st = 1 - a1 - a2
        st[st<0] = 0
        stateList[[5]][j,k,l] = mean(st)*transStateProbs[[5]][l,k-1]/margStateProbs[5,1]
        
        ##############################################################################################
        ### m_{i}(0) < \sqrt{500} & m_{i}(0) > \sqrt{350} & m_{i}(t) \in S_{l} (in the l-th state) ###
        ##############################################################################################
        startPoint = as.numeric(startingValues[startingValues$l1==6 & startingValues$l2 == l & startingValues$k == k,
                                               c("startb1","startb2","startb3")])
        
        # If the marginal probability of the state is too low, use Gibbs sampling
        # otherwise we use HMC sampling
        
        if (transStateProbs[[6]][l,k-1] < 1e-6)
        {
          # Gibbs sampling: we carry out a simple approach checking convergence of the Gibbs sampling
          ok = F
          jj = 0 
          while(ok == F)
          {
            jj = jj + 1
            if (l == 1)
            {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,]),
                                             offset = c( states2[6] - mu[1],- states[6] + mu[1],states2[l] - mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else if (l == 7){
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],-zTime[k,]),
                                             offset = c( states2[6] - mu[1],- states[6] + mu[1], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else {
              samp = sample_from_constraints(linear_part = rbind(zTime[1,],-zTime[1,],zTime[k,],-zTime[k,]),
                                             offset = c( states2[6] - mu[1],- states[6] + mu[1],states2[l] - mu[k], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            }
            
            #diag(cov(samp))
            ok = !all(diag(cov(samp)) < 1e-4)
          }
          
        } else{
          
          # Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          if (l == 1)
          {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,]), 
                         g = c(mu[1]-states[6],states2[6]-mu[1], states2[l]-mu[k]), burn.in = nburn)
          } else if (l==7) {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],zTime[k,]), 
                         g = c(mu[1]-states[6],states2[6]-mu[1], mu[k]-states[l] ), burn.in = nburn)
          } else {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[1,],-zTime[k,],zTime[k,]), 
                         g = c(mu[1]-states[6],states2[6]-mu[1], states2[l]-mu[k], mu[k]-states[l] ), burn.in = nburn)
          }
          
        }
        
        a1 = cif1evalVec(k,samp)
        a2 = cif2evalVec(k,samp)
        
        # Overall survival function for group 1
        st = 1 - a1 - a2
        st[st<0] = 0
        stateList[[6]][j,k,l] = mean(st)*transStateProbs[[6]][l,k-1]/margStateProbs[6,1]
        
        ######################################################################
        ### m_{i}(0) > \sqrt{500} & m_{i}(t) \in S_{l} (in the l-th state) ###
        ######################################################################
        startPoint = as.numeric(startingValues[startingValues$l1==7 & startingValues$l2 == l & startingValues$k == k,
                                               c("startb1","startb2","startb3")])
        
        # If the marginal probability of the state is too low, use Gibbs sampling
        # otherwise we use HMC sampling
        
        if (transStateProbs[[7]][l,k-1] < 1e-6)
        {
          # Gibbs sampling: we carry out a simple approach checking convergence of the Gibbs sampling
          ok = F
          jj = 0 
          while(ok == F)
          {
            jj = jj + 1
            if (l == 1)
            {
              samp = sample_from_constraints(linear_part = rbind(-zTime[1,],zTime[k,]),
                                             offset = c( - states[7] + mu[1],states2[l] - mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else if (l == 7){
              samp = sample_from_constraints(linear_part = rbind(-zTime[1,],-zTime[k,]),
                                             offset = c( - states[7] + mu[1], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            } else {
              samp = sample_from_constraints(linear_part = rbind(-zTime[1,],zTime[k,],-zTime[k,]),
                                             offset = c( - states[7] + mu[1],states2[l] - mu[k], - states[l] + mu[k]),
                                             mean_param = rep(0,q),
                                             covariance = D,
                                             initial_point = startPoint,
                                             ndraw = nMC,
                                             burnin = nburn)
            }
            
            #diag(cov(samp))
            ok = !all(diag(cov(samp)) < 1e-4)
          }
          
        } else{
          
          # Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          # using Hamiltonian Monte Carlo
          if (l == 1)
          {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[k,]), 
                         g = c(mu[1]-states[7], states2[l]-mu[k]), burn.in = nburn)
          } else if (l==7) {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],zTime[k,]), 
                         g = c(mu[1]-states[7], mu[k]-states[l] ), burn.in = nburn)
          } else {
            samp <- rtmg(n = nMC, M = Dinv, r = rep(0,q), initial = startPoint,
                         f = rbind(zTime[1,],-zTime[k,],zTime[k,]), 
                         g = c(mu[1]-states[7], states2[l]-mu[k], mu[k]-states[l] ), burn.in = nburn)
          }
          
        }
        
        a1 = cif1evalVec(k,samp)
        a2 = cif2evalVec(k,samp)
        
        # Overall survival function for group 1
        st = 1 - a1 - a2
        st[st<0] = 0
        stateList[[7]][j,k,l] = mean(st)*transStateProbs[[7]][l,k-1]/margStateProbs[7,1]
        
      }
      
    }
    
    # We keep estimates for the population CIF for cause 1 conditional on baseline states
    # based on the remaining Importance sampling Monte Carlo estimates!!!
    for (i in 1:7)
    {
      ss1[j,,i] <- 1 - ss2[j,,i] - rowSums(stateList[[i]][j,,])
    }
    
    ###############################
    ### Population average CIFs ###
    ###############################
    ss1Overall[j,] = colSums(t(ss1[j,,])*margStateProbs[,1])
    ss2Overall[j,] = colSums(t(ss2[j,,])*margStateProbs[,1])
    
    #####################################
    ### Marginal latent marker states ###
    #####################################
    
    # State 1
    state1Overall[j,] = stateList[[1]][j,,1]*margStateProbs[1,1]
    for(i in 2:7)
    {
      state1Overall[j,] = state1Overall[j,] + stateList[[i]][j,,1]*margStateProbs[i,1]
    }
    state1Overall[j,1] = margStateProbs[1,1]
    
    # State 2
    state2Overall[j,] = stateList[[1]][j,,2]*margStateProbs[1,1]
    for(i in 2:7)
    {
      state2Overall[j,] = state2Overall[j,] + stateList[[i]][j,,2]*margStateProbs[i,1]
    }
    state2Overall[j,1] = margStateProbs[2,1]
    
    # State 3
    state3Overall[j,] = stateList[[1]][j,,3]*margStateProbs[1,1]
    for(i in 2:7)
    {
      state3Overall[j,] = state3Overall[j,] + stateList[[i]][j,,3]*margStateProbs[i,1]
    }
    state3Overall[j,1] = margStateProbs[3,1]
    
    # State 4
    state4Overall[j,] = stateList[[1]][j,,4]*margStateProbs[1,1]
    for(i in 2:7)
    {
      state4Overall[j,] = state4Overall[j,] + stateList[[i]][j,,4]*margStateProbs[i,1]
    }
    state4Overall[j,1] = margStateProbs[4,1]
    
    # State 5
    state5Overall[j,] = stateList[[1]][j,,5]*margStateProbs[1,1]
    for(i in 2:7)
    {
      state5Overall[j,] = state5Overall[j,] + stateList[[i]][j,,5]*margStateProbs[i,1]
    }
    state5Overall[j,1] = margStateProbs[5,1]
    
    # State 6
    state6Overall[j,] = stateList[[1]][j,,6]*margStateProbs[1,1]
    for(i in 2:7)
    {
      state6Overall[j,] = state6Overall[j,] + stateList[[i]][j,,6]*margStateProbs[i,1]
    }
    state6Overall[j,1] = margStateProbs[6,1]
    
    # State 7
    state7Overall[j,] = stateList[[1]][j,,7]*margStateProbs[1,1]
    for(i in 2:7)
    {
      state7Overall[j,] = state7Overall[j,] + stateList[[i]][j,,7]*margStateProbs[i,1]
    }
    state7Overall[j,1] = margStateProbs[7,1]
    
    
    # 1 - ss1Overall[j,] - ss2Overall[j,]
    # state1Overall[j,] + state2Overall[j,] + state3Overall[j,] + state4Overall[j,] + state5Overall[j,] + state6Overall[j,] + state7Overall[j,]
    
    print(j)
  }
  
  ### Checking
  round( 100*colMeans( 1 - ss1[,,1] - ss2[,,1] ),3)
  round( 100*colMeans( apply(stateList[[1]],1:2,sum) ),3)
  
  round( 100*apply( 1 - ss1[,,7] - ss2[,,7] ,2,sd),3)
  round( 100*apply( apply(stateList[[7]],1:2,sum) ,2,sd),3)
  
  round( 100*colMeans( ss1[,,1] ),3)
  round( 100*colMeans( 1 - ss2[,,1] - apply(stateList[[1]],1:2,sum) ),3)
  
  round( 100*apply( ss1[,,1],2,sd ),3)
  round( 100*apply( 1 - ss2[,,1] - apply(stateList[[1]],1:2,sum),2,sd ),3)
  ###
  
  ss1Overall[,1] = 0
  
  # Results for population CIFs
  sumCif1 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  sumCif2 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  
  # Results for latent marker probabilities
  sumS1 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  sumS2 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  sumS3 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  sumS4 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  sumS5 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  sumS6 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  sumS7 = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
  
  sumCif1$Mean = colMeans(ss1Overall)
  sumCif1$Median = apply(ss1Overall,2,median)
  sumCif1$SD = apply(ss1Overall,2,sd)
  sumCif1[,c("LB","UB")] = t(apply(ss1Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumCif2$Mean = colMeans(ss2Overall)
  sumCif2$Median = apply(ss2Overall,2,median)
  sumCif2$SD = apply(ss2Overall,2,sd)
  sumCif2[,c("LB","UB")] = t(apply(ss2Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumS1$Mean = colMeans(state1Overall)
  sumS1$Median = apply(state1Overall,2,median)
  sumS1$SD = apply(state1Overall,2,sd)
  sumS1[,c("LB","UB")] = t(apply(state1Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumS2$Mean = colMeans(state2Overall)
  sumS2$Median = apply(state2Overall,2,median)
  sumS2$SD = apply(state2Overall,2,sd)
  sumS2[,c("LB","UB")] = t(apply(state2Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumS3$Mean = colMeans(state3Overall)
  sumS3$Median = apply(state3Overall,2,median)
  sumS3$SD = apply(state3Overall,2,sd)
  sumS3[,c("LB","UB")] = t(apply(state3Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumS4$Mean = colMeans(state4Overall)
  sumS4$Median = apply(state4Overall,2,median)
  sumS4$SD = apply(state4Overall,2,sd)
  sumS4[,c("LB","UB")] = t(apply(state4Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumS5$Mean = colMeans(state5Overall)
  sumS5$Median = apply(state5Overall,2,median)
  sumS5$SD = apply(state5Overall,2,sd)
  sumS5[,c("LB","UB")] = t(apply(state5Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumS6$Mean = colMeans(state6Overall)
  sumS6$Median = apply(state6Overall,2,median)
  sumS6$SD = apply(state6Overall,2,sd)
  sumS6[,c("LB","UB")] = t(apply(state6Overall,2,quantile,probs = c(0.025,0.975)))
  
  sumS7$Mean = colMeans(state7Overall)
  sumS7$Median = apply(state7Overall,2,median)
  sumS7$SD = apply(state7Overall,2,sd)
  sumS7[,c("LB","UB")] = t(apply(state7Overall,2,quantile,probs = c(0.025,0.975)))
  
  ##########################################################################
  ### Results for population averaged CIFs conditional on baseline state ###
  ##########################################################################
  
  # Cause 1
  sumCif1List = list()
  for (i in 1:7)
  {
    sumCif1List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumCif1List[[i]]$Mean = colMeans(ss1[,,i])
    sumCif1List[[i]]$Median = apply(ss1[,,i],2,median)
    sumCif1List[[i]]$SD = apply(ss1[,,i],2,sd)
    sumCif1List[[i]][,c("LB","UB")] = t(apply(ss1[,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  # Cause 2
  sumCif2List = list()
  for (i in 1:7)
  {
    sumCif2List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumCif2List[[i]]$Mean = colMeans(ss2[,,i])
    sumCif2List[[i]]$Median = apply(ss2[,,i],2,median)
    sumCif2List[[i]]$SD = apply(ss2[,,i],2,sd)
    sumCif2List[[i]][,c("LB","UB")] = t(apply(ss2[,,i],2,quantile,probs = c(0.025,0.975)))
  }
  
  # Transition probabilities to latent marker states by baseline state 1
  sumS1List = list()
  for (i in 1:7)
  {
    sumS1List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumS1List[[i]]$Mean = colMeans(stateList[[1]][,,i])
    sumS1List[[i]]$Median = apply(stateList[[1]][,,i],2,median)
    sumS1List[[i]]$SD = apply(stateList[[1]][,,i],2,sd)
    sumS1List[[i]][,c("LB","UB")] = t(apply(stateList[[1]][,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  # Transition probabilities to latent marker states by baseline state 2
  sumS2List = list()
  for (i in 1:7)
  {
    sumS2List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumS2List[[i]]$Mean = colMeans(stateList[[2]][,,i])
    sumS2List[[i]]$Median = apply(stateList[[2]][,,i],2,median)
    sumS2List[[i]]$SD = apply(stateList[[2]][,,i],2,sd)
    sumS2List[[i]][,c("LB","UB")] = t(apply(stateList[[2]][,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  # Transition probabilities to latent marker states by baseline state 3
  sumS3List = list()
  for (i in 1:7)
  {
    sumS3List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumS3List[[i]]$Mean = colMeans(stateList[[3]][,,i])
    sumS3List[[i]]$Median = apply(stateList[[3]][,,i],2,median)
    sumS3List[[i]]$SD = apply(stateList[[3]][,,i],2,sd)
    sumS3List[[i]][,c("LB","UB")] = t(apply(stateList[[3]][,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  # Transition probabilities to latent marker states by baseline state 4
  sumS4List = list()
  for (i in 1:7)
  {
    sumS4List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumS4List[[i]]$Mean = colMeans(stateList[[4]][,,i])
    sumS4List[[i]]$Median = apply(stateList[[4]][,,i],2,median)
    sumS4List[[i]]$SD = apply(stateList[[4]][,,i],2,sd)
    sumS4List[[i]][,c("LB","UB")] = t(apply(stateList[[4]][,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  # Transition probabilities to latent marker states by baseline state 5
  sumS5List = list()
  for (i in 1:7)
  {
    sumS5List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumS5List[[i]]$Mean = colMeans(stateList[[5]][,,i])
    sumS5List[[i]]$Median = apply(stateList[[5]][,,i],2,median)
    sumS5List[[i]]$SD = apply(stateList[[5]][,,i],2,sd)
    sumS5List[[i]][,c("LB","UB")] = t(apply(stateList[[5]][,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  # Transition probabilities to latent marker states by baseline state 6
  sumS6List = list()
  for (i in 1:7)
  {
    sumS6List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumS6List[[i]]$Mean = colMeans(stateList[[6]][,,i])
    sumS6List[[i]]$Median = apply(stateList[[6]][,,i],2,median)
    sumS6List[[i]]$SD = apply(stateList[[6]][,,i],2,sd)
    sumS6List[[i]][,c("LB","UB")] = t(apply(stateList[[6]][,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  # Transition probabilities to latent marker states by baseline state 7
  sumS7List = list()
  for (i in 1:7)
  {
    sumS7List[[i]] = data.frame(times = tt,Mean  = NA,Median = NA,SD = NA, LB = NA,UB = NA)
    sumS7List[[i]]$Mean = colMeans(stateList[[7]][,,i])
    sumS7List[[i]]$Median = apply(stateList[[7]][,,i],2,median)
    sumS7List[[i]]$SD = apply(stateList[[7]][,,i],2,sd)
    sumS7List[[i]][,c("LB","UB")] = t(apply(stateList[[7]][,,i],2,quantile,probs = c(0.025,0.975),na.rm = T))
  }
  
  sumS1List[[1]] + sumS1List[[2]] + sumS1List[[3]] + sumS1List[[4]] + sumS1List[[5]] + sumS1List[[6]] + sumS1List[[7]]
  1 - sumCif1List[[1]] - sumCif2List[[1]]
  
  
  sumS1$Mean + sumS2$Mean + sumS3$Mean + sumS4$Mean + sumS5$Mean + sumS6$Mean + sumS7$Mean 
  1 - sumCif1$Mean - sumCif2$Mean
  
  
  out = list()
  out$ss1 = ss1
  out$ss2 = ss2
  out$stateList = stateList
  
  out$sumCif1 = sumCif1
  out$sumCif2 = sumCif2
  out$sumS1 = sumS1
  out$sumS2 = sumS2
  out$sumS3 = sumS3
  out$sumS4 = sumS4
  out$sumS5 = sumS5
  out$sumS6 = sumS6
  out$sumS7 = sumS7
  
  out$sumCif1List = sumCif1List
  out$sumCif2List = sumCif2List
  out$sumS1List = sumS1List
  out$sumS2List = sumS2List
  out$sumS3List = sumS3List
  out$sumS4List = sumS4List
  out$sumS5List = sumS5List
  out$sumS6List = sumS6List
  out$sumS7List = sumS7List
  
  out$ss1Overall = ss1Overall
  out$ss2Overall = ss2Overall
  out$state1Overall = state1Overall
  out$state2Overall = state2Overall
  out$state3Overall = state3Overall
  out$state4Overall = state4Overall
  out$state5Overall = state5Overall
  out$state6Overall = state6Overall
  out$state7Overall = state7Overall
  
  out$nMC = nMC
  out$seed = seed
  out$nburn = nburn
  out$thin = thin
  out$dimsamp = dim(samp)
  out$Sdraws = Sdraws
  return(out)
  
}


