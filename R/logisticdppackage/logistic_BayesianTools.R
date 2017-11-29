library('BayesianTools')


sampleSize = 30
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
y <-  1 * x + 1*x^2 + rnorm(n=sampleSize,mean=0,sd=10)
plot(x,y, main="Test Data")

likelihood1 <- function(param){
  pred = param[1] + param[2]*x + param[3] * x^2
  singlelikelihoods = dnorm(y, mean = pred, sd = 1/(param[4]^2), log = T)
  return(sum(singlelikelihoods))  
}


setUp1 <- createBayesianSetup(likelihood1, lower = c(-5,-5,-5,0.01), upper = c(5,5,5,30))

settings = list(iterations = 15000,  message = FALSE)
out <- runMCMC(bayesianSetup = setUp1, sampler = "Metropolis", settings = settings)

gelmanDiagnostics(out)
summary(out)

plot(out, start = 1000)
correlationPlot(out, start = 1000)
marginalPlot(out, start = 1000)

MAP(out)
DIC(out)
WAIC(out) # requires special definition of the likelihood, see help
marginalLikelihood(out)





library(mixtools)
library(reshape2)
library(DPpackage)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggthemes)


source('GenData.R')
source('logisticFunctions.R')
source('initialize_parms_log.R')
dyn.load('splineFunctionsC.so')

### Generate Data
xsig=.707; ysig=2.121
temp=genData(type=6,nobs=500)
xobs=temp$xobs; yobs=temp$yobs; xcens=temp$xcens; ycens=temp$ycens; 
trueFit=temp$true; xtrue=temp$xtrue; densTrue=temp$densTrue; xgrid=temp$xgrid
nobs=length(xobs)

parms=initialize_parms_log(xobs,yobs,xcens,ycens)
coefs=parms$coefs
xtrue=parms$xtrue
ytrue=getYtrue(coefs,xtrue)

nobs=length(xobs)
onevec <- c(rep(1,nobs))
xcensl=rep(0,nobs)
xcensl[xcens==-1] = 1
xcensu=rep(0,nobs)
xcensu[xcens==1] = 1
ycensl=rep(0,nobs)
ycensl[ycens==-1] = 1
ycensu=rep(0,nobs)
ycensu[ycens==1] = 1

parms=updateDens(xtrue,state=0,first=TRUE,niter=1000,xgrid)
# mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state;
# dens=parms$densEst; k=length(p); groups=parms$groups
k=parms$k; densEst=parms$densEst; state=parms$state


likelihood1=function(param){
  
  
  mb = (2*param[3]*param[4])/(param[3]+param[4])
  fx = 1/(1+exp(-mb*(param[2]-xtrue)))
  ytrue=param[1]*(fx*exp(param[3]*(param[2]-xtrue))+(1-fx)*exp(param[4]*
           (param[2]-xtrue)))/(1+fx*exp(param[3]*(param[2]-xtrue))+(1-fx)*
                                 exp(param[4]*(param[2]-xtrue)))
  
  
  pmic <- pnorm(xobs,xtrue,xsig)*(1-xcensu) + onevec*xcensu - 
    pnorm(xobs-1,xtrue,xsig)*(1-xcensl)
  pmic[pmic<.001]=.0001
  pmic[is.nan(pmic)]=.00000001
  
  
  pdia <- pnorm(yobs+0.5,ytrue,ysig)*(1-ycensu) + onevec*ycensu -
    pnorm(yobs-0.5,ytrue,ysig)*(1-ycensl)
  pdia[pdia<.001]=.0001
  pdia[is.nan(pdia)]=.00000001
  
  ### MIC-DIA
  llike <- sum(log(pmic) + log(pdia))
  
  # ### population
  # pxden=log(approx(xgrid,dens,xout=xtrue,rule=2)$y)
  
  ### prior
  priorcoefs <- sum(dlnorm(param,0,10,log=TRUE))
  
  return(llike+priorcoefs)
  
}

#### Y

likelihood1=function(param){
  
  
  mb = (2*param[3]*param[4])/(param[3]+param[4])
  fx = 1/(1+exp(-mb*(param[2]-xtrue)))
  ytrue=param[1]*(fx*exp(param[3]*(param[2]-xtrue))+(1-fx)*exp(param[4]*
         (param[2]-xtrue)))/(1+fx*exp(param[3]*(param[2]-xtrue))+(1-fx)*
         exp(param[4]*(param[2]-xtrue)))
  
  pdia <- pnorm(yobs+0.5,ytrue,ysig)*(1-ycensu) + onevec*ycensu -
    pnorm(yobs-0.5,ytrue,ysig)*(1-ycensl)
  pdia[pdia<.001]=.0001
  pdia[is.nan(pdia)]=.00000001
  
  
  
  
  # llik = sum(dnorm(ytrue, mean = yobs, sd = ysig, log = T))
  
  return(sum(log(pdia)))
  
}

setUp1 <- createBayesianSetup(likelihood1, lower = c(20,0,0,0), upper = c(60,5,5,5))
settings = list(iterations = 25000,  message = FALSE)
out <- runMCMC(bayesianSetup = setUp1, sampler = "Metropolis", settings = settings)



summary(out)
plot(out, start = 1000)


coefs=MAP(out)$parametersMAP

mb = (2*coefs[3]*coefs[4])/(coefs[3]+coefs[4])
fx = 1/(1+exp(-mb*(coefs[2]-xgrid)))
ytrue=coefs[1]*(fx*exp(coefs[3]*(coefs[2]-xgrid))+(1-fx)*exp(coefs[4]*
     (coefs[2]-xgrid)))/(1+fx*exp(coefs[3]*(coefs[2]-xgrid))+(1-fx)*
                           exp(coefs[4]*(coefs[2]-xgrid)))
plot(xobs,yobs,pch=16)
lines(xgrid,trueFit)

lines(xgrid,ytrue,col=2)



#### X

likelihood1=function(xtrue){
  
  

  
  pmic <- pnorm(xobs,xtrue,xsig)*(1-xcensu) + onevec*xcensu - 
    pnorm(xobs-1,xtrue,xsig)*(1-xcensl)
  pmic[pmic<.001]=.0001
  pmic[is.nan(pmic)]=.00000001
  
  
  pxden=log(approx(xgrid,densEst,xout=xtrue,rule=2)$y)
  # pxden[is.nan(pxden)]=.00000001
  
  
  
  # llik = sum(dnorm(ytrue, mean = yobs, sd = ysig, log = T))
  
  return(sum(log(pmic))+sum(pxden))
  
}

setUp1 <- createBayesianSetup(likelihood1,lower=xobs-1,upper=xobs+1)
settings = list(iterations = 2000,  message = FALSE)
out <- runMCMC(bayesianSetup = setUp1, sampler = "DREAM", settings = settings)



summary(out)
plot(out, start = 1000)


coefs=MAP(out)$parametersMAP

mb = (2*coefs[3]*coefs[4])/(coefs[3]+coefs[4])
fx = 1/(1+exp(-mb*(coefs[2]-xgrid)))
ytrue=coefs[1]*(fx*exp(coefs[3]*(coefs[2]-xgrid))+(1-fx)*exp(coefs[4]*
                                                               (coefs[2]-xgrid)))/(1+fx*exp(coefs[3]*(coefs[2]-xgrid))+(1-fx)*
                                                                                     exp(coefs[4]*(coefs[2]-xgrid)))
plot(xobs,yobs,pch=16)
lines(xgrid,trueFit)

lines(xgrid,ytrue,col=2)





