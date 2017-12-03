bayesian_mon_errors_logistic_diag=function(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid,xsig=.707,ysig=2.121,
                                           numIter=26000,burnin=10000,thin=4){

  fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  MICDens=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  coefMat=matrix(nrow=numIter,ncol=4)
  acceptCoef=rep(NA,numIter)
  xtrue_sav=matrix(nrow=numIter-burnin,ncol=length(xobs))
  nobs=length(xobs)


  for(iter in 1:numIter){

    if(iter%%1000==0) cat('Iteration: ',iter,'\n')

    ## Update Density
    if(iter %% 500==0 | iter==1){
      if(iter==1){
        parms=updateDens(xtrue,state=0,first=TRUE,niter=1000,xgrid)
        mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state;
        dens=parms$densEst; k=length(p); groups=parms$groups
      }else{
        parms=updateDens(xtrue,state,first=FALSE,niter=125,xgrid)
        mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state;
        dens=parms$densEst; k=length(p); groups=parms$groups
      }
    }

    ### Update xtrue (and corresponding ytrue)
    if(iter %% 15==0 | iter==1){
      parms=updateXtrueLog(xtrue,ytrue,coefs,xcens,ycens,xsig,ysig,xobs,yobs,
           nobs,mu,sigma,k,p,groups)
      xtrue=parms$xtrue; ytrue=as.numeric(parms$ytrue)
    }

    ### xupdate Coefs
    parms=updateLogCoef(iter,coefs,xobs,yobs,xtrue,xcens,ycens,ytrue,xsig,ysig,coefMat)
    acceptCoef[iter]=parms$accept; coefs=parms$coefs; ytrue=parms$ytrue


    ###save
    coefMat[iter,]=coefs
    if(iter>burnin){
      fit=getylogtrue(coefs,xgrid)
      fitMat[iter-burnin,]=fit
      MICDens[iter-burnin,]=dens
      xtrue_sav[iter-burnin,]=xtrue
    }

  }

  ### Thin
  MICDens=MICDens[seq(1,nrow(MICDens),by=thin),]
  fitMat=fitMat[seq(1,nrow(fitMat),by=thin),]

  acceptCoef=acceptCoef[burnin:numIter]
  coefMat=coefMat[burnin:numIter,]

  return(list(MICDens=MICDens,fitMat=fitMat,acceptCoef=acceptCoef,coefMat=coefMat,xtrue_sav=xtrue_sav))

}





