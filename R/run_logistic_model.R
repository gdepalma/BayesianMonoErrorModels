
run_logistic_model=function(xobs,yobs,xcens,ycens){

  ### Initialize
  nobs=length(xobs)
  xgrid=seq(min(xobs)-1,max(xobs)+1,length=1000)
  parms=initialize_parms_logistic(xobs,yobs,xcens)
  xtrue=parms$xtrue
  coefs=parms$coefs
  ytrue=getylogtrue(coefs,xtrue)

  ### Run Model
  parms=bayesian_mon_errors_logistic(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid)
  MICDens=parms$MICDens; fitMat=parms$fitMat

  return(list(MICDens=MICDens,fitMat=fitMat))
}

