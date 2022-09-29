############################################################
### Code for simulating one dataset from the SPM-1 model ###
### and fitting both the SPM-1 and SPM-2 model           ###
############################################################

set.seed(15)
# The working directory
setwd("C:/works/Methodological/jointCumINCpaper/ExampleCode/programs")

if (!require("mvtnorm")) install.packages("mvtnorm")
if (!require("lspline")) install.packages("lspline")
# library(mvtnorm)
# library(lspline)

##############################
### Simulation of the data ###
##############################

# Number of individuals in the dataset
N = 2000

# Percentage of doubly sampled patients (among non-censored observations)
Rper = 0.15

# Probability of correct clasification of cause 1
pi11 = 0.75

# Probability of correct clasification of cause 2
pi22 = 0.90

########################################################################################################
### Assumed population-averaged evolution: linear splines with knots at 1 and 5 years since baseline.### 
### the same for the random effects but no random effect for the slope after 5 years                 ###  
########################################################################################################

# Fixed effects
Lbeta.true = c(12.85,6.03,0.77,0)

# Covariance matrix of the random effects
Dtrue = matrix(c(25.09,-5.08,-2.65,-5.08,10.18,1.09,-2.65,1.09,0.85),nr = 3,nc = 3,byrow = F)

# Measurement times (identical for all individuals)
ni = length(seq(0,10,by = 0.5))
times = rep(seq(0,10,by = 0.5),N)


# ID variable: individual identifier
id = rep(1:N,each = ni)
ord = rep(1:ni,N)

# Data frame with the full data (i.e. in the absence of failure/censoring)
dataFull = data.frame(id = id,times = times,ord = ord)
maxtime = max(dataFull$times)

# Design matrices of the fixed and random effects
x = cbind(1,lspline(times,knots = c(1,5)))
z = x[,-4] # No random effect for the slope after 5 years
n = nrow(x)

# Simulation of the random effects
b = rmvnorm(n = N,sigma = Dtrue)

# The ``true'' marker values for all patients
zb = rowSums(b[id,]*z)
mt = x %*% Lbeta.true + zb

# Within-individual errors
epsilon = rnorm(n,sd = sqrt(8.14))

# The observed marker values
sqcd4 = mt + epsilon
dataFull$sqcd4 = sqcd4
dataFull$mt = mt

# Parameters of the subdistribution hazard models
p1 = 0.62
p2 = 0.70
phi1 = 0.25
phi2 = 0.13

# Association parameters
alpha1 = - 0.16
alpha2 = - 0.02

# Binary covariate affecting the subdistribution hazards
group = rbinom(N,size = 1,prob = 0.5)
bg1 = 0.15
bg2 = -0.15

# Function evaluating the subdistribution hazard for cause 1 at any time point
subhaz1 = function(time)
{
  aa = cbind(1,lspline(time,knots = c(1,5)))
  mtime = aa %*% Lbeta.true
  mtime = mtime + aa[,-4] %*% b[i,]
  
  base = p1*phi1*exp(-phi1*time)/(1-p1*(1-exp(-phi1*time)))
  
  base*exp(alpha1*(mtime-7) + bg1*group[i])
}

# Function evaluating the CIF for cause 1 at any time point
F1 = function(time)
{
  # Subdistribution cumulative hazard
  csubhaz = integrate(subhaz1,0,time)$val
  
  # Cumulative incidence
  1 - exp(-csubhaz)
}
F1 = Vectorize(F1)

# Function evaluating the subdistribution hazard for cause 2 at any time point
subhaz2 = function(time)
{
  aa = cbind(1,lspline(time,knots = c(1,5)))
  mtime = aa %*% Lbeta.true
  mtime = mtime + aa[,-4] %*% b[i,]
  
  base = p2*phi2*exp(-phi2*time)/(1-p2*(1-exp(-phi2*time)))
  
  base*exp(alpha2*(mtime-7) + bg2*group[i])
}

# Function evaluating the CIF for cause 1 at any time point
F2 = function(time)
{
  # Subdistribution cumulative hazard
  csubhaz = integrate(subhaz2,0,time)$val
  
  # Cumulative incidence
  1 - exp(-csubhaz)
}
F2 = Vectorize(F2)

# Data frame with information on failure times/events
dropoutData = data.frame(time = rep(NA,N), status = NA,group = group)

# Function used to simulate the overall failure time
ff = function(time)
{
  log( F1(time) + F2(time) ) - log(u[i])
}

# Simulate the overall failure time and then the ``true'' event type!
dropoutData$time = NA
u = runif(N)

# Simulate first the overall failure time and then the event type condition on the failure time
for (i in 1:N)
{ 
  # Simulation of the failure time
  fit = try( uniroot(ff,lower = .Machine$double.eps,upper = maxtime,tol = .Machine$double.eps),silent = T)
  
  # If an error occurs, it is implied that there is no solution of S(t)=u within [0,maxtime]
  # thus, data should be right-censored at maxtime
  if ( class(fit) %in% "try-error")
  {
    dropoutData$time[i] = maxtime
    dropoutData$status[i] = 0
  } else{
    
    if (fit$root==0)
    {
      fit$root = .Machine$double.eps
    }
    dropoutData$time[i] = fit$root
    
    # The probability of the``true'' event type for cause 1 is equal to the ratio of the cause-specific hazard for cause 1
    # to the sum of both cause-specific hazards
    # It can be shown that this ratio remains the same if we replace the cause-specific hazard with the derivatives of the corresponding CIFs
    derF1 = (1-F1(fit$root))*subhaz1(fit$root)
    derF2 = (1-F2(fit$root))*subhaz2(fit$root)
    pcause1 = derF1/(derF1+derF2)
    
    # Simulate the cause of failure
    cause1 = rbinom(n = 1,size = 1,prob = pcause1)
    dropoutData$status[i] = 1*cause1 + 2*(1-cause1)
    
  }
  print(i)
}

# We then simulate the ``reported'' event type condition on the ``true'' event type
# pi11 = Pr(Kobs=1|Ktrue=1)
# pi22 = Pr(Kobs=2|Ktrue=2)
dropoutData$statusObs = NA
dropoutData$statusObs[dropoutData$status==1] = ifelse(rbinom(sum(dropoutData$status==1),size = 1,prob = pi11)==1,yes = 1,no = 2)
dropoutData$statusObs[dropoutData$status==2] = ifelse(rbinom(sum(dropoutData$status==2),size = 1,prob = pi22)==1,yes = 2,no = 1)
dropoutData$statusObs[is.na(dropoutData$statusObs)] = 0

# Noninformative right-censoring time simulated from the exponential distribution
# We assume that right censoring is always correctly classified, i.e. Kobs = 0 <=> Ktrue = 0
dropoutData$cens = rexp(N,rate = 0.025)
dropoutData$delta = 1*(dropoutData$time<dropoutData$cens)
dropoutData$time = pmin(dropoutData$time,dropoutData$cens)
dropoutData$status = dropoutData$status*dropoutData$delta
dropoutData$statusObs = dropoutData$statusObs*dropoutData$delta

# Indicator of a doubly sampled patient
# Double sampling is performed in a random subset of the individuals observed to fail from any cause, with probability Rper
dropoutData$Robs = NA
dropoutData$Robs[dropoutData$status!=0] = rbinom(sum(dropoutData$status!=0),size = 1,prob = Rper)

# The status variable that would be available in real data applications,
# i.e. for doubly sampled individuals, the ``true'' status variable is also available
dropoutData$statusAnalysis = dropoutData$statusObs
dropoutData$statusAnalysis[which(dropoutData$Robs==1)] = dropoutData$status[which(dropoutData$Robs==1)]

# Merge with the longitudinal dataset
dataFull$etimes = rep(dropoutData$time,each = ni)
dataFull$mis = dataFull$times > dataFull$etimes
dataFull$status = dropoutData$status[dataFull$id]
dataFull$statusObs = dropoutData$statusObs[dataFull$id]
dataFull$statusAnalysis = dropoutData$statusAnalysis[dataFull$id]
dataFull$Robs = dropoutData$Robs[dataFull$id]
dataFull$group = rep(dropoutData$group,each = ni)

dataobs = dataFull[dataFull$mis==F,]
dataobs$status[which(dataobs$Robs==0)] <- NA

# Remove all objects
aa <- ls()
aa <- aa[aa!="dataobs"]
rm(list = aa)
rm(aa)

##########################
### Fitting the models ###
##########################

# Import functions for M1
source("SurvModel_Subdistr.R")
source("JointModel_Subdistr.R")

# Import functions for M2
source("SurvModel_IncrCumInc.R")
source("JointModel_IncrCumInc.R")

# List with the assumed prior distributions 
prior = list()
prior$Lc0 = diag(4)*100
prior$Lmu0 = rep(0,nrow(prior$Lc0))
prior$Smu0 = rep(0,17)
prior$Sc0 = diag(17)*100
prior$df = 3
prior$A = 3*diag(c(25,5,5))
prior$lambda1 = 0.01
prior$lambba2 = 0.01
prior$a1 = 1
prior$a2 = 1
prior$b1 = 1
prior$b2 = 1

# List with the initial values of the MCMC algorithm
startingValues = list()
startingValues$Lbeta = c(12.85,6.03,0.79,0.01)
startingValues$omega = 0.10
startingValues$D = matrix(c(25.26,-5.62,-2.44,-5.62,9.93,1.13,-2.44,1.13,0.85),nr = 3,nc = 3)

if (!require("nlme")) install.packages("nlme")
if (!require("survival")) install.packages("survival")
if (!require("Matrix")) install.packages("Matrix")
if (!require("cmprsk")) install.packages("cmprsk")
if (!require("splines2")) install.packages("splines2")
if (!require("numDeriv")) install.packages("numDeriv")
if (!require("matrixcalc")) install.packages("matrixcalc")
if (!require("tmg")) install.packages("https://cran.r-project.org/src/contrib/Archive/tmg/tmg_0.3.tar.gz", repo=NULL, type="source")
if (!require("restrictedMVN")) install.packages("restrictedMVN")
if (!require("RcppEigen")) install.packages("RcppEigen")
if (!require("alabama")) install.packages("alabama")

# library(nlme)
# library(survival)
# library(Matrix)
# library(cmprsk)
# library(splines2)
# library(numDeriv)
# library(matrixcalc)
# library(tmg)
# library(restrictedMVN)

# Fitting the linear mixed model
fitlme = try(lme(sqcd4 ~ lspline(times,knots = c(1,5)),random = ~ lspline(times,knots = c(1,5))[,-3]|id,data = dataobs,method = "ML",
                 control = list(apVar = T,opt = "optim" , returnObject = F, maxIter = 100, msMaxIter = 100, niterEM = 150,msVerbose = T)),
             silent = T)
summary(fitlme)

# Derivatives of the design matrices of the fixed and the random effects at the observed survival times
dataid = dataobs[dataobs$ord==1,]

# Fit the cox models
fitCoxCause1 = coxph(Surv(etimes,statusObs == 1) ~ group,data = dataid,x = T)
summary(fitCoxCause1)

fitCoxCause2 = coxph(Surv(etimes,statusObs == 2) ~ group,data = dataid,x = T)
summary(fitCoxCause2)

set.seed(15)
# Fit the SPM-1 model (subdistribution hazard)
fitProp_Subd = try(jointModel_Subdistr_Misc(fitlme,fitCoxCause1,fitCoxCause2,prior,startingValues = startingValues,
                                            ndraw = 3500,thin = 3,nburn = 200,iterMaxSurv_1 = 2,iterMaxSurv_2 = 1,store.b = T,
                                            nknotsCause1 = 2,nknotsCause2 = 3,scaleVarSurv = 1,useGauleg = T,Gauleg_points = 30),
                   silent = T)


round(fitProp_Subd$summary$Longitudinal_Process,3)
round(fitProp_Subd$summary$sumSurvCause1,3)
round(fitProp_Subd$summary$sumSurvCause2,3)
round(fitProp_Subd$summary$sumMiscl,3)

# Marginal DIC criterion
fitDICsubd = DICSubdistr_Misc(fitProp_Subd,nMC = 200,thinDIC = 7)
fitDICsubd$dic
fitDICsubd$eff

####################################################################
### Predict population CIFs and latent marker states for group 1 ###
####################################################################
newdata1 = data.frame(group = c(1),etimes = 0,statusObs = 0)
newdata2 = data.frame(group = c(1),etimes = 0,statusObs = 0)

tt = seq(0,10,by = 2)
fit1 = predictLatentStates_Subdistr(fitProp_Subd,newdata1,newdata2,tt = tt,nMC = 1000,thin = 10)

####################################################################
### Predict population CIFs and latent marker states for group 0 ###
####################################################################
newdata1 = data.frame(group = c(0),etimes = 0,statusObs = 0)
newdata2 = data.frame(group = c(0),etimes = 0,statusObs = 0)

fit0 = predictLatentStates_Subdistr(fitProp_Subd,newdata1,newdata2,tt = tt,nMC = 1000,thin = 10)

round(100*fit1$sumCif1,1)
round(100*fit1$sumCif2,1)
round(100*fit1$sumS1,1)
round(100*fit1$sumS2,1)
round(100*fit1$sumS3,1)
round(100*fit1$sumS4,1)
round(100*fit1$sumS5,1)
round(100*fit1$sumS6,1)
round(100*fit1$sumS7,1)

set.seed(15)
# Fit the SPM-2 model
fitProp_IncrCumInc = try(jointModel_IncrCumInc_Misc(fitlme,fitCoxCause1,fitCoxCause2,prior,startingValues = startingValues,
                                                    ndraw = 3500,thin = 3,nburn = 200,
                                                    alpha1 = 1,alpha2 = 1,iterMaxSurv_1 = 2,iterMaxSurv_2 = 1,store.b = T,
                                                    nknotsCause1 = 2,nknotsCause2 = 3,scaleVarSurv = 1,useGauleg = T,Gauleg_points = 30),
                         silent = T)

round(fitProp_IncrCumInc$summary$Longitudinal_Process,3)
round(fitProp_IncrCumInc$summary$sumSurvCause1,3)
round(fitProp_IncrCumInc$summary$sumSurvCause2,3)
round(fitProp_IncrCumInc$summary$sumMiscl,3)

fitDICIncrCinc = DICIncrCumInc_Misc(fitProp_IncrCumInc,nMC = 200,thinDIC = 7)
fitDICIncrCinc$eff
fitDICIncrCinc$dic

####################################################################
### Predict population CIFs and latent marker states for group 1 ###
####################################################################
newdata1 = data.frame(group = c(1),etimes = 0,statusObs = 0)
newdata2 = data.frame(group = c(1),etimes = 0,statusObs = 0)

fit1 = predictLatentStates_IncrCumInc(fitProp_IncrCumInc,newdata1,newdata2,tt = tt,nMC = 1000,thin = 10)

####################################################################
### Predict population CIFs and latent marker states for group 0 ###
####################################################################
newdata1 = data.frame(group = c(0),etimes = 0,statusObs = 0)
newdata2 = data.frame(group = c(0),etimes = 0,statusObs = 0)

fit0 = predictLatentStates_IncrCumInc(fitProp_IncrCumInc,newdata1,newdata2,tt = tt,nMC = 1000,thin = 10)

round(100*fit1$sumCif1,1)
round(100*fit1$sumCif2,1)
round(100*fit1$sumS1,1)
round(100*fit1$sumS2,1)
round(100*fit1$sumS3,1)
round(100*fit1$sumS4,1)
round(100*fit1$sumS5,1)
round(100*fit1$sumS6,1)
round(100*fit1$sumS7,1)
