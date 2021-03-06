
a1=read_csv(file='test/data1.csv')

parms=parse_file(a1,logMIC=FALSE)
xobs=parms$xobs; yobs=parms$yobs; xcens=parms$xcens; ycens=parms$ycens

### Initialize
nobs=length(xobs)
xgrid=seq(min(xobs)-1,max(xobs)+1,length=1000)
parms=initialize_parms_logistic(xobs,yobs,xcens)
xtrue=parms$xtrue
coefs=parms$coefs
ytrue=getylogtrue(coefs,xtrue)


### Run
parms=bayesian_mon_errors_logistic_diag(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid)
MICDens=parms$MICDens; fitMat=parms$fitMat; acceptCoef=parms$acceptCoef
coefMat=parms$coefMat; xtrue_sav=parms$xtrue_sav
logistic_diagnostic(xgrid,xobs,yobs,xcens,ycens,MICDens,fitMat,acceptCoef,xtrue_sav,coefMat)

parms=bayesian_mon_errors_logistic(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid)
MICDens=parms$MICDens; fitMat=parms$fitMat

