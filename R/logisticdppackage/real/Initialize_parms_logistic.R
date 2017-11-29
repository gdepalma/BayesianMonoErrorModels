initialize_parms_logistic=function(xobs,yobs,xcens){
  
  xtrue = rnorm(length(xobs),xobs-.5,xsig)
  xtrue[xcens==1]=xtrue[xcens==1]+1
  xtrue[xcens==-1]=xtrue[xcens==-1]-1
  
  #### Initial Logistic coefficients
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
  
  
  return(list(xtrue=xtrue,coefs=coefs))
}