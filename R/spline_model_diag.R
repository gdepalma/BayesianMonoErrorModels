bayesian_mon_errors_spline_diag=function(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid,lowept,upperept,designMatrixGrid,
                                    xsig=.707,ysig=2.121,numIter=26000,burnin=10000,thin=4){

  fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  MICDens=matrix(nrow=numIter-burnin,ncol=length(xgrid))
  coefMat=matrix(nrow=numIter,ncol=length(coefs))
  acceptCoef=rep(NA,numIter)
  xtrue_sav=matrix(nrow=numIter-burnin,ncol=length(xobs))
  smoothAccept=rep(NA,numIter)
  smoothParam_sav=rep(NA,numIter)
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
        parms=updateDens(xtrue,state,first=FALSE,niter=250,xgrid)
        mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state;
        dens=parms$densEst; k=length(p); groups=parms$groups
      }
    }

    ### Update xtrue (and corresponding ytrue)
    parms=updateSplineXtrue(xtrue,ytrue,mu,sigma,p,k,groups,knotseq,bases,lowept,upperept,coefs,xcens,
                            ycens,xsig,ysig,xobs,yobs,nobs)
    xtrue=parms$xtrue; ytrue=as.numeric(parms$ytrue); acceptCoef[iter]=parms$acceptCoef;



    ### Update Coefs
    parms=updateSplineCoefs(smoothParam,coefs,iter,bases,xtrue,
                            ytrue,knotseq,xobs,yobs,xcens,ycens,xsig,ysig,coefMat)
    coefs=parms$coefs; ytrue=as.numeric(parms$ytrue)

    ### update smoothness parameter
    parms=updateSmoothParm(smoothParam,coefs)
    smoothParam=parms$smooth; smoothAccept[iter]=parms$accept

    smoothParam_sav[iter,]=smoothParam
    coefMat[iter,]=log(coefs)
    if(iter>burnin){
      fitMat[iter-burnin,]=coefs%*%t(designMatrixGrid)
      MICDens[iter-burnin,]=densEst
      xtrue_sav[iter-burnin,]=xtrue
    }

  }

  ### Thin
  MICDens=MICDens[seq(1,nrow(MICDens),by=thin),]
  fitMat=fitMat[seq(1,nrow(fitMat),by=thin),]

  acceptCoef=acceptCoef[burnin:numIter]
  smoothAccept=smoothAccept[burnin:numIter]
  smoothParam=smoothParam[burnin:numIter]


  return(list(MICDens=MICDens,fitMat=fitMat,acceptCoef=acceptCoef,smoothAccept=smoothAccept,
              smoothParam=smoothParam,coefMat=coefMat))

}





