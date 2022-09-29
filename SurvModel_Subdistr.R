########################################################
### Restricted cubic splines for the baseline hazard ###
########################################################


survModel_Subdistr_addConstr = function(fitCoxCause1,fitCoxCause2,fitlme,Lbeta = fixed.effects(fitlme),reffects = as.matrix(ranef(fitlme)),
                                        nknotsCause1 = 4,nknotsCause2 = 4,timeVar = "times",
                                        maxit = 50,epsilon = 1e-6,alpha = 0.05,startIter = 5,full.knotsCause1 = NULL,full.knotsCause2 = NULL,
                                        useGauleg = F,Gauleg_points = 30)
{
  # Check if this is a coxph object
  if ( !( class(fitCoxCause1)=="coxph" & class(fitCoxCause2)=="coxph") )
  {
    break
  }
  
  if ( (!is.null(nknotsCause1) &  !is.null(full.knotsCause1)) |  (!is.null(nknotsCause2) &  !is.null(full.knotsCause2)) )
  {
    break
  }
  
  if ( ncol(fitCoxCause1$y) != 2 | ncol(fitCoxCause2$y) != 2  )
  {
    break
  }
  
  if ( sum(fitCoxCause1$y[,1] - fitCoxCause2$y[,1]) != 0 )
  {
    break
  }
  
  # Number of regression coefficients
  pCause1 = ncol(model.matrix(fitCoxCause1)) 
  nCause1 = nrow(model.matrix(fitCoxCause1))
  pCause2 = ncol(model.matrix(fitCoxCause2)) 
  nCause2 = nrow(model.matrix(fitCoxCause2))
  
  if (ncol(fitCoxCause1$y) == 2)
  {
    times = fitCoxCause1$y[,1]
    deltaCause1 = fitCoxCause1$y[,2]
    
    names(times)=NULL;names(deltaCause1)=NULL
  }
  
  if (ncol(fitCoxCause2$y) == 2)
  {
    #times = fitCoxCause1$y[,1]
    deltaCause2 = fitCoxCause2$y[,2]
    
    names(times)=NULL;names(deltaCause2)=NULL
  }
  
  # Failure indicator
  delta = deltaCause1 + deltaCause2
  indcens = which(delta == 0)
  indDelta1 = which(deltaCause1==1)
  indDelta2 = which(deltaCause2==1)
  
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
  
  # The design matrix for cause 1
  xCause1 = model.matrix(fitCoxCause1)
  colnames(xCause1) = NULL
  rownames(xCause1) = NULL
  
  # The design matrix for cause 2
  xCause2 = model.matrix(fitCoxCause2)
  colnames(xCause2) = NULL
  rownames(xCause2) = NULL
  
  #################################################################
  ### Obtaining starting values for the spline basis parameters ###
  #################################################################
  
  delta1_temp = deltaCause1
  delta1_temp[deltaCause2==1] = 2
  
  delta2_temp = deltaCause2
  delta2_temp[deltaCause1==1] = 2
  
  fit1 = crr(ftime = times,fstatus = delta1_temp,cov1 = xCause1)
  fit2 = crr(ftime = times,fstatus = delta2_temp,cov1 = xCause2)
  
  # Baseline cumulative subdistribution hazards
  baseCsubhaz1 = -log(1 - predict.crr(fit1,cov1 = 0)[,2])
  baseCsubhaz2 = -log(1 - predict.crr(fit2,cov1 = 0)[,2])
  
  basesubhaz1 = diff(baseCsubhaz1)
  basesubhaz2 = diff(baseCsubhaz2)
  
  times1 = predict.crr(fit1,cov1 = 0)[-1,1]
  times2 = predict.crr(fit2,cov1 = 0)[-1,1]

  y1 = log( basesubhaz1 )
  #y1 = log(l1*exp(-kappa1*times1))
  
  plot(times1,y1,type = "l")
  
  spltemp1 = bSpline(times1,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T)
  res1 = lm(y1 ~ -1 + spltemp1)
  lines(times1,fitted(res1),col = "red")
  
  y2 = log( basesubhaz2 )
  #y2 = log(l2*exp(-kappa2*times2))
  
  plot(times2,y2,type = "l")
  
  spltemp2 = bSpline(times2,knots = knotsCause2,Boundary.knots = bknotsCause2,intercept = T)
  res2 = lm(y2 ~ -1 + spltemp2)
  lines(times2,fitted(res2),col = "red")
  
  rm(list = c("baseCsubhaz1","baseCsubhaz2","basesubhaz1","basesubhaz2","delta1_temp","delta2_temp",
              "times1","times2","spltemp1","spltemp2","y1","y2") )
  
  # Parameter vector for cause 1 
  betasCause1 = fit1$coef
  gammasCause1 = coef(res1)
  lbetasCause1 = length(betasCause1)
  lgammasCause1 = length(gammasCause1)
  paramCause1 = c(betasCause1,gammasCause1)
  names(paramCause1) = NULL
  lpCause1 = length(paramCause1)
  
  # Parameter vector for cause 2 
  betasCause2 = fit2$coef
  gammasCause2 = coef(res2)
  lbetasCause2 = length(betasCause2)
  lgammasCause2 = length(gammasCause2)
  paramCause2 = c(betasCause2,gammasCause2)
  names(paramCause2) = NULL
  lpCause2 = length(paramCause2)
  
  paramCause1 = c(betasCause1,gammasCause1)
  paramCause2 = c(betasCause2,gammasCause2)
  
  param = c(paramCause1,paramCause2)
  
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
  
  # Comparing results of Gauss-Kronrod with the integrate function of R
  ff = function(x){
    temp = bSpline(x,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T) %*% gammasCause1
    
    ss = cbind(1,lspline(x,knots = c(0.5,2))) %*% Lbeta
    
    ss = ss + cbind(1,lspline(x,knots = c(0.5,2)))[,-4] %*% reffects[i,]
    
    c(exp(temp + sum(xCause1[i,-lbetasCause1]*betasCause1[-lbetasCause1]) + betasCause1[lbetasCause1]*ss))
  }
  # i = which(times==10)[1]
  # integrate(ff,0,times[i])
  # 
  # t1 = 0.5*times[i]
  # t1*sum(wws*ff(t1*pps+t1))
  
  t1 = 0.5*times
  t1Long = rep(t1,each = length(pps))
  ppsLong = rep(pps,nCause1)
  wwsLong = rep(wws,nCause1)
  ind = rep(1:nCause1,each = length(pps))
  
  # Additional time points used to integrate the subditribution hazard
  timesForInt = t1Long*ppsLong + t1Long
  
  # Construct the spline basis for cause 1
  splBasisCause1Long = bSpline(timesForInt,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T)
  splBasisCause1 = bSpline(times,knots = knotsCause1,Boundary.knots = bknotsCause1,intercept = T)

  # Construct the spline basis for cause 2
  splBasisCause2Long = bSpline(timesForInt,knots = knotsCause2,Boundary.knots = bknotsCause2,intercept = T)
  splBasisCause2 = bSpline(times,knots = knotsCause2,Boundary.knots = bknotsCause2,intercept = T)

  # Make predictions for the marker values at the survival times
  id = as.numeric(as.character(fitlme$groups[,1]))
  dataLmeid = fitlme$data[!duplicated(id),]
  dataLmeid[,timeVar] = times
  xTime = model.matrix(formula(fitlme),data = dataLmeid)
  zTime = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeid)
  mtpred = xTime %*% Lbeta + rowSums(zTime*reffects)
  xCause1[,length(betasCause1)] = mtpred
  xCause2[,length(betasCause2)] = mtpred
  
  # Make predictions for the marker values at additional time points (used to integrate the subdistribution hazard) 
  dataLmeid = dataLmeid[ind,]
  dataLmeid[,timeVar] = timesForInt
  
  xTime = model.matrix(formula(fitlme),data = dataLmeid)
  zTime = model.matrix(formula(fitlme$modelStruct$reStr)[[1]],dataLmeid)
  mtpred = xTime %*% Lbeta + rowSums(zTime*reffects[ind,])
  
  xCause1Long = as.matrix(xCause1[ind,])
  xCause1Long[,length(betasCause1)] = mtpred
  
  xCause2Long = as.matrix(xCause2[ind,])
  xCause2Long[,length(betasCause2)] = mtpred
  
  # Indices of id's
  ids = list()
  for (i in 1:nCause1){ids[[i]] = which(ind==i)}
  
  # Combined design matrices of the survival submodel (including the B-splines basis)
  xCause1Comb = cbind(xCause1,splBasisCause1)
  xCause1CombLong = cbind(xCause1Long,splBasisCause1Long)
  
  xCause2Comb = cbind(xCause2,splBasisCause2)
  xCause2CombLong = cbind(xCause2Long,splBasisCause2Long)
  
  # Failure indicators for the `extended` design matrices of the survival submodels
  nLong = nCause1*length(pps)
  deltaCause1Long = rep( deltaCause1, each = length(pps) )
  deltaCause2Long = rep( deltaCause2, each = length(pps) )
  deltaLong = deltaCause1Long + deltaCause2Long
  
  # The objective function 
  f = function(param,computeDeriv = F)
  {
    paramCause1 = param[1:lpCause1]
    paramCause2 = param[-c(1:lpCause1)]
    
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
    qq = t1Long*wwsLong*explinpredCause1Long
    int1 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq[x]))))
    
    # int_{0}^{T_{i}} exp{ u_{k}(s) + a_{k}m_{i}(s)ds }
    qq = t1Long*wwsLong*explinpredCause2Long
    int2 = length(pps)*unlist(lapply(ids,function(x) .Internal(mean(qq[x]))))
    
    # Cumulative incidence functions
    cincCause1 = 1 - exp(-int1)
    cincCause2 = 1 - exp(-int2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = - int1 + linpredCause1
    DercincCause2 = - int2 + linpredCause2
    
    logl1 = sum(DercincCause1[indDelta1]) + sum(DercincCause2[indDelta2]) 

    logl2 = log(1-cincCause1-cincCause2)[indcens]
    
    out = - logl1 - sum(logl2)

    if (computeDeriv)
    {
      
    }
    return(out)
    
  }

  # The objective function: returns NaNs when the sum of the CIFs is greater than 1
  fsurvEval = function(paramSurv,computeDeriv = F)
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
    cincCause1 = 1 - exp(-int1)
    cincCause2 = 1 - exp(-int2)
    
    # Derivatives of the cumulative incidence functions
    DercincCause1 = - int1 + linpredCause1
    DercincCause2 = - int2 + linpredCause2
    
    logl = sum(DercincCause1[indDelta1]) + sum(DercincCause2[indDelta2]) + sum(log(1-cincCause1-cincCause2)*(1-delta))
    return(-logl)
    
  }
  
  # Function specifying the boundness constraints: For each individual the sum of the CIFs
  # should be less than 1 at the observed survival time
  hin = function(paramSurv,computeDeriv = F)
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
    cincCause1 = 1 - exp(-int1)
    cincCause2 = 1 - exp(-int2)
    
    1 - cincCause1 - cincCause2
  }
  
  # The gradient of the objective function
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
    cincCause1 = 1 - exp(-int1)
    cincCause2 = 1 - exp(-int2)
    
    # Compute the score vector
    qq1Mat = xCause1CombLong*qq1
    qq2Mat = xCause2CombLong*qq2
    
    tt1 = (1-cincCause1)/(1-cincCause1-cincCause2)
    tt2 = (1-cincCause2)/(1-cincCause1-cincCause2)
    aa1 = rep( tt1, each = length(pps) )
    aa2 = rep( tt2, each = length(pps) )
    
    score1 = .colSums(deltaCause1*xCause1Comb,nCause1,lpCause1,F) - .colSums( deltaCause1Long*qq1Mat,nLong,lpCause1,F) 
    score1 = score1 - .colSums( (1 - deltaLong)*qq1Mat*aa1,nLong,lpCause1,F)
    
    score2 = .colSums(deltaCause2*xCause2Comb,nCause1,lpCause2,F) - .colSums( deltaCause2Long*qq2Mat,nLong,lpCause2,F ) 
    score2 = score2 - .colSums( (1 - deltaLong)*qq2Mat*aa2,nLong,lpCause2,F)
    
    # Score vector
    gradient = c(score1,score2)
    
    return(-gradient)
  }
  
  # The jacobian of the nonlinear constraints
  hin.jac = function(paramSurv,iter = 5)
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
    cincCause1 = 1 - exp(-int1)
    cincCause2 = 1 - exp(-int2)
    
    # Compute the score vector
    qq1Mat = xCause1CombLong*qq1
    qq2Mat = xCause2CombLong*qq2
    
    aa1 = rep( 1-cincCause1, each = length(pps) )
    aa2 = rep( 1-cincCause2, each = length(pps) )
    
    mat = cbind(qq1Mat*aa1,qq2Mat*aa2)
    
    for (i in 1:nCause1)
    {
      jacMat[i,] = -.colSums(mat[ids[[i]],],length(pps),lpCause1 + lpCause2,F) 
    }
    jacMat
  }
  
  
  jacMat = matrix(nr = nCause1, nc = lpCause1 + lpCause2)
  
  # To get plausible starting values
  nn = is.finite(f(param))
  tt = 0
  
  while(nn == F)
  {
    tt = tt + 1
    
    # Descease parameters for the baseline cumulative incidence functions
    param[(lbetasCause1+1):lpCause1] = param[(lbetasCause1+1):lpCause1] - 0.5
    param[lpCause1+lbetasCause2+1] = param[lpCause1+lbetasCause2+1] - 0.5
    nn = is.finite(f(param))
  }
  
  cd <-
    function (x, f, ..., eps = 1e-05) {
      n <- length(x)
      res <- numeric(n)
      ex <- pmax(abs(x), 1)
      for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[i] <- diff.f / diff.x
      }
      res
    }
  
  cd.vec <-
    function (x, f, ..., eps = 1e-04) {
      n <- length(x)
      res <- matrix(0, n, n)
      ex <- pmax(abs(x), 1)
      for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[, i] <- diff.f / diff.x
      }
      0.5 * (res + t(res))
    }
  
  grad.param = function(x){cd(x,f)}
  hess.param = function(x){cd.vec(x,grad.param)}
  
  # Get starting values without using the boundness constraints for individuals who have failed
  fit = optim(param,fn = f,gr = fsurvDer,control = list(trace = 1,maxit = 200,reltol = 1e-12),method = "BFGS")
  param = fit$par
  
  # To get plausible starting values
  nn = is.finite(fsurvEval(param))
  tt = 0
  
  while(nn == F)
  {
    tt = tt + 1
    
    # Descease parameters for the baseline cumulative incidence functions
    param[(lbetasCause1+1):lpCause1] = param[(lbetasCause1+1):lpCause1] - 0.05
    param[lpCause1+lbetasCause2+1] = param[lpCause1+lbetasCause2+1] - 0.05
    nn = is.finite(fsurvEval(param))
  }
  
  # Incorporate the boundess constraint on the CIFs for all subjects, using the alabama package
  ans <- alabama::constrOptim.nl(par = param, fn = fsurvEval, gr = fsurvDer, hin = hin,hin.jac = hin.jac,control.outer = list(trace = T)) 
  param = ans$par
  
  # # 
  # # fit = nlm(f,param,iterlim = 200,print.level = 2)
  # # param = fit$estimate
  # 
  # #fit = nlminb(param,f,gradient = grad.param,control = list(trace = 10,eval.max = 500,iter.max = 350))
  # #param = fit$par
  # fit$hess = hess.param(param)
  # 
  # # Store the history of the NR algorithm
  # # History of NR - algorithm
  # diff.coef = Inf
  # 
  # it = 0
  # 
  # NR.hist = matrix(nrow = maxit + 1,ncol = 3)
  # colnames(NR.hist) = c("Iteration","Diff_coef","logl")
  # NR.hist[1,] = c(it,diff.coef,-f(param))
  # 
  # # History of coefficients
  # coef.histCause1 = matrix(nrow = maxit + 1,ncol = 1 + length(param[1:lpCause1]))
  # names.coefCause1 = c( c(names(coef(fitCoxCause1)),paste0("gammas",1:lgammasCause1) )) 
  # colnames(coef.histCause1) = c("Iteration",names.coefCause1)
  # coef.histCause1[1,] = c(it,param[1:lpCause1])
  # 
  # # History of coefficients
  # coef.histCause2 = matrix(nrow = maxit + 1,ncol = 1 + length(param[-c(1:lpCause1)]))
  # names.coefCause2 = c( c(names(coef(fitCoxCause2)),paste0("gammas",1:lgammasCause2) )) 
  # colnames(coef.histCause2) = c("Iteration",names.coefCause2)
  # coef.histCause2[1,] = c(it,param[-c(1:lpCause1)])
  # 
  # param.t = param
  # # Start of NR algorithm
  # 
  # while (diff.coef > epsilon & it < maxit)
  # {
  #   # Update iteration index
  #   it = it+1
  #   
  #   param = param.t
  # 
  #   #########################################
  #   ### Update parameters using 1 NR step ###
  #   #########################################
  #   
  #   #res = f(param,computeDeriv = T)
  #   
  #   param.t = try(param - solve(hess.param(param),grad.param(param)),silent = T)
  #   
  #   if (class(param.t) %in% "try-error")
  #   {
  #     fit$param = fit$par
  #     fit$lpcause1 = lpCause1
  #     fit$lpcause2 = lpCause2
  #     return(fit)
  #   }
  #   
  #   # Difference in coefficients
  #   diff.coef = sqrt( sum( (param-param.t)^2 ) )
  #   
  #   NR.hist[it+1,] = c(it,diff.coef,f(param.t))
  #   coef.histCause1[it+1,] = c(it,param.t[1:lpCause1])
  #   coef.histCause2[it+1,] = c(it,param.t[-c(1:lpCause1)])
  #   
  # }
  # param = param.t
  # #res = f(param,computeDeriv = T)
  hess = jacobian(fsurvDer,param)
  covm = solve(hess)
  
  sumCause1 = matrix(nr = lpCause1,nc = 4)
  rownames(sumCause1) = c( c(names(coef(fitCoxCause1)),paste0("gammas",1:lgammasCause1) ))
  colnames(sumCause1) = c("Estimate","SE","LB","UB")
  
  sumCause1[,1] = param[1:lpCause1]
  sumCause1[,2] = sqrt(diag(covm)[1:lpCause1])
  sumCause1[,3] = sumCause1[,1] - qnorm(1-alpha/2)*sumCause1[,2]
  sumCause1[,4] = sumCause1[,1] + qnorm(1-alpha/2)*sumCause1[,2]
  
  sumCause2 = matrix(nr = lpCause2,nc = 4)
  rownames(sumCause2) = c( c(names(coef(fitCoxCause2)),paste0("gammas",1:lgammasCause2) )) 
  colnames(sumCause2) = c("Estimate","SE","LB","UB")
  
  sumCause2[,1] = param[-c(1:lpCause1)]
  sumCause2[,2] = sqrt(diag(covm)[-c(1:lpCause1)])
  sumCause2[,3] = sumCause2[,1] - qnorm(1-alpha/2)*sumCause2[,2]
  sumCause2[,4] = sumCause2[,1] + qnorm(1-alpha/2)*sumCause2[,2]
  
  # End NR iterations
  # Exclude na's
  # NR.hist = NR.hist[1:(it+1),]
  # coef.histCause1 = coef.histCause1[1:(it+1),]
  # coef.histCause2 = coef.histCause2[1:(it+1),]
  
  # Prepare output
  out = list()
  
  if (F)
  {
    out$note = paste("Convergence failed at",epsilon,
                     "Euclidean distance for coefficients",
                     ", may need more iterations than",maxit)
  } else {
    out$note = ans
    out$not2 = fit
  
  }
  
  #out$it = it
  # out$NR.hist = NR.hist
  # out$coef.histCause1 = coef.histCause1
  # out$coef.histCause2 = coef.histCause2
  out$summaryCause1 = sumCause1
  out$summaryCause2 = sumCause2
  out$var = covm
  out$hess = hess
  out$knotsCause1 = knotsCause1
  out$knotsCause2 = knotsCause2
  out$bknotsCause1 = bknotsCause1
  out$bknotsCause2 = bknotsCause2
  out$param = param
  out$lpcause1 = lpCause1
  out$lpcause2 = lpCause2
  
  return(out)
}

