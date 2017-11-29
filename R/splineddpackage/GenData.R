
Ispline=function(intknots,lowept,upperept)
{
  k=3
  #determine knot sequence
  knotseq=c(rep(lowept,k+1),intknots,rep(upperept,k+1))
  numbases=length(knotseq)-k-2
  
  #create matrix of bases
  bases=matrix(NA,nrow=numbases,ncol=2)
  for(i in 1:numbases) bases[i,]=c(knotseq[i+1],knotseq[i+k+1])
  
  return(list(bases=bases,knotseq=knotseq))
}


getIsplineC=function(xtrue,knotseq,bases){
  
  numBases=nrow(bases)
  lxtrue=length(xtrue)
  
  mat=rep(0,numBases*lxtrue)
  
  storage.mode(mat) <- "double"
  storage.mode(knotseq) <- "double"
  storage.mode(xtrue) <- "double"
  storage.mode(numBases) <- "integer"
  storage.mode(lxtrue) <- "integer"
  temp=.C("getIspline",xtrue,lxtrue,knotseq,mat,numBases)
  designMatrix=matrix(temp[[4]],ncol=numBases,nrow=lxtrue)
  return(designMatrix)
  
}



genData=function(type,nobs=600){
  
  
  xcens=rep(0,nobs)
  ycens=rep(0,nobs)
  
  xsig=.707; ysig=2.121
  xgrid=seq(-12,12,length=1000)
  
  
  #1
  if(type==1){
    popmn=c(-6,3); popstd=c(2,.7); popprob=c(.8,.2)
    xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
    xobs=ceiling(xtrue+rnorm(nobs,0,xsig))
    
    ytrue=30-2*xtrue
    true=30-2*xgrid 
  }
  
  #2
  if(type==2){
    popmn=c(-4,0,4); popstd=c(.7,.7,.7); popprob=c(1/3,1/3,1/3)
    coef=c(29,1.17,.48)
    xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
    xobs=ceiling(xtrue+rnorm(nobs,0,xsig))
    ytrue=coef[1]*(exp(coef[2]-coef[3]*xtrue)/(1+exp(coef[2]-coef[3]*xtrue)))
    true=coef[1]*(exp(coef[2]-coef[3]*xgrid)/(1+exp(coef[2]-coef[3]*xgrid)))
  }
  
  #3
  if(type==3){
    
    popmn=c(-4,0,4); popstd=c(.7,.7,.7); popprob=c(1/3,1/3,1/3)
    xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
    xobs=ceiling(xtrue+rnorm(nobs,0,xsig))
    
    icoefsT=c(1,1,20,1,20,1)
    parms=Ispline(c(-3,0,1),-6,6)
    basesT=parms$bases;
    knotseqT=parms$knotseq
    designMatrix=getIsplineC(xtrue,knotseqT,basesT)
    ytrue = icoefsT%*%t(designMatrix)+5
    designMatrix=getIsplineC(xgrid,knotseqT,basesT)
    true = icoefsT%*%t(designMatrix)+5  
  }
  
  #4
  if(type==4){
    
    popmn=c(-3,0,3); popstd=c(1,1,1); popprob=c(.5,.3,.2)
    xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
    xobs=ceiling(xtrue+rnorm(nobs,0,xsig))
    
    
    icoefsT=c(1,10,1,25,1,1,10,1)
    parms=Ispline(c(-4,-2,0,2,4),-6,6)
    basesT=parms$bases;
    knotseqT=parms$knotseq
    designMatrix=getIsplineC(xtrue,knotseqT,basesT)
    ytrue = icoefsT%*%t(designMatrix)
    designMatrix=getIsplineC(xgrid,knotseqT,basesT)
    true = icoefsT%*%t(designMatrix)
  }
  if(type==5){
    popmn=c(-4.6,-2,1); popstd=c(1.1,1.5,1.5); popprob=c(.6,.2,.2)
    xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
    xobs=ceiling(xtrue+rnorm(nobs,0,xsig))
    
    coef=c(35,1.17,.1,1.2)
    mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
    fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
    ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
        (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
        exp(coef[4]*(coef[2]-xtrue)))
    fx = 1/(1+exp(-mb*(coef[2]-xgrid)))
    true=coef[1]*(fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*exp(coef[4]*
          (coef[2]-xgrid)))/(1+fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*
           exp(coef[4]*(coef[2]-xgrid)))
  }
  
  ### censored values
  if (type==6){
    popmn=c(-4.6,-2,1); popstd=c(1.1,1.5,1.5); popprob=c(.6,.2,.2)
    xtrue=rnormmix(n=nobs,lambda=popprob,mu=popmn,sigma=popstd)
    xobs=ceiling(xtrue+rnorm(nobs,0,xsig))
    
    coef=c(35,1.17,.1,1.2)
    mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
    fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
    ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
        (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
       exp(coef[4]*(coef[2]-xtrue)))
    fx = 1/(1+exp(-mb*(coef[2]-xgrid)))
    true=coef[1]*(fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*exp(coef[4]*
       (coef[2]-xgrid)))/(1+fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*
        exp(coef[4]*(coef[2]-xgrid)))
    
  }
  
  yobs=round(ytrue+rnorm(nobs,0,ysig))
  
  densTrue=rep(0,length(xgrid))
  for(i in 1:length(xgrid))
    densTrue[i]=sum(popprob*dnorm(xgrid[i],popmn,popstd))
  densTrue=densTrue/sum(densTrue)
  
  return(list(xobs=xobs,yobs=as.numeric(yobs),xcens=xcens,ycens=ycens,
              trueFit=true,xtrue=xtrue,densTrue=densTrue,xgrid=xgrid))
}
