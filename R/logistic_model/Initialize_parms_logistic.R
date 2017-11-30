initialize_parms_logistic=function(xobs,yobs,xcens){

  ## Initialize xtrue
  xtrue=rep(NA,length(xobs))

  n=sum(xcens==0)
  n1=sum(xcens==-1)
  n2=sum(xcens==1)

  xtrue[xcens==0]=xobs[xcens==0]-runif(n,0,1)
  xtrue[xcens==-1]=xobs[xcens==-1]-runif(n1,.5,2)
  xtrue[xcens==1]=xobs[xcens==1]+runif(n2,.5,1.5)

  #### Initial Logistic coefficients
  nls_logistic=function(coef1,coef2,coef3,coef4){
    mb = (2*coef3*coef4)/(coef3+coef4)
    fx = 1/(1+exp(-mb*(coef2-xtrue)))
    ytrue=coef1*(fx*exp(coef3*(coef2-xtrue))+(1-fx)*exp(coef4*
          (coef2-xtrue)))/(1+fx*exp(coef3*(coef2-xtrue))+(1-fx)*
                             exp(coef4*(coef2-xtrue)))
    return(ytrue)
  }
  fit=tryCatch({
    nlsLM(yobs~nls_logistic(coef1,coef2,coef3,coef4),start=list(coef1=max(yobs)-2,coef2=.2,coef3=.1,coef4=.1))
  },
  error=function(cond) {
    message(cond)
    stop("Logistic initial parameters not converging.")
  }
  )
  coefs=coef(fit)


  return(list(xtrue=xtrue,coefs=coefs))
}
