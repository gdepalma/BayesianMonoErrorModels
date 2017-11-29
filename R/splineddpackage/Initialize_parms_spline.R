initialize_parms_spline=function(xobs,yobs,xcens,ycens,xgrid){
  
  ## Initial xtrue
  xtrue=rep(NA,length(xobs))
  for(i in 1:length(xobs)){
    if(xcens[i]==0 & xcens[i]==0){
      xtrue[i]=xobs[i]-runif(1,0,1)
    }else if(xcens[i]==-1 & xcens[i]==0){
      xtrue[i]=xobs[i]-runif(1,.5,1.5)
    }else if(xcens[i]==0 & xcens[i]==1){
      xtrue[i]=xobs[i]+runif(1,0,1)
    }
  }
  
  ### Initalize Isplines
  lowept=min(xobs)-.5
  upperept=max(xobs)+.5
  if(sum(xcens==1)>0 | sum(ycens==1)>0)
    upperept=max(xobs)+4
  if(sum(xcens==-1)>0 | sum(ycens==-1)>0)
    lowept=min(xobs)-4
  dist=1
  parms=Ispline(seq(min(xobs),max(xobs),by=dist),lowept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq
  
  #### Initial spline coefficients
  designMatrix=getIsplineC(xtrue,knotseq,bases)
  coefs=findInitialSpline(xtrue,bases,knotseq,yobs,designMatrix)
  ytrue=coefs%*%t(designMatrix)
  
  designMatrixGrid=getIsplineC(xgrid,knotseq,bases)
  
  
  
  return(list(coefs=coefs,xtrue=xtrue,bases=bases,knotseq=knotseq,ytrue=ytrue,
              lowept=lowept,upperept=upperept,designMatrixGrid=designMatrixGrid))
}