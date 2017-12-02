
getLLK=function(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig)
{
  nobs=length(xobs)
  onevec <- c(rep(1,nobs))
  xcensl=rep(0,nobs)
  xcensl[xcens==-1] = 1
  xcensu=rep(0,nobs)
  xcensu[xcens==1] = 1

  pmic <- pnorm(xobs,xtrue,xsig)*(1-xcensu) + onevec*xcensu -
    pnorm(xobs-1,xtrue,xsig)*(1-xcensl)
  pmic[pmic<.001]=.0001

  ycensl=rep(0,nobs)
  ycensl[ycens==-1] = 1
  ycensu=rep(0,nobs)
  ycensu[ycens==1] = 1

  pdia <- pnorm(yobs+0.5,ytrue,ysig)*(1-ycensu) + onevec*ycensu -
    pnorm(yobs-0.5,ytrue,ysig)*(1-ycensl)
  pdia[pdia<.001]=.0001

  llike <- log(pmic) + log(pdia)

  return(llike)
}


getPOPLLK=function(x,mu,sigma,p,k,groups)
{

  pxden <- log(p[groups])+dnorm(x,mu[groups],sigma[groups],log=T)
  return(pxden)

}


updateDens=function(xtrue,state,first,niter,xgrid){

  nburn <- niter/2
  nsave <- niter
  nskip <- 10
  ndisplay <- 500
  mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
  priorDP <- list(alpha=.1,m1=rep(0,1),psiinv2=solve(diag(0.5,1)),
                  nu1=4,nu2=4,tau1=1,tau2=100)

  if(first==TRUE){
    fit <- DPdensity(y=xtrue,prior=priorDP,mcmc=mcmc,status=TRUE,grid=xgrid)
  }else{
    fit <- DPdensity(y=xtrue,prior=priorDP,mcmc=mcmc,status=FALSE,state=state,grid=xgrid)
  }
  state=fit$state
  k=fit$state$ncluster
  groups=fit$state$ss

  mu=rep(NA,k); sigma=rep(NA,k); p=rep(NA,k);
  for(i in 1:k){
    mu[i]=mean(xtrue[groups==i])
    sigma[i]=sd(xtrue[groups==i])
    p[i]=sum(groups==i)/length(groups)
  }
  densEst=fit$dens/sum(fit$dens)


  return(list(mu=mu,sigma=sigma,p=p,state=state,densEst=densEst,groups=groups))
}
