initialize_parms_log=function(xobs,yobs,xcens,ycens){
  
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
  
  #### Initial Logistic coefficients
  
  ### Need to add try catch
  nls_logistic=function(xtrue,coef){
    mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
    fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
    ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
        (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
       exp(coef[4]*(coef[2]-xtrue)))
    return(ytrue)
  }
  fit=nls(yobs~nls_logistic(xtrue,coef),start=list(coef=c(30,.2,.5,1)),
          lower=c(0.01,0.01,0.01,0.01),algorithm='port',trace=FALSE,control=list(maxiter=1000,tol = 1e-04))
  coefs <- coef(fit)
  
  
  return(list(coefs=coefs,xtrue1=xtrue))
}