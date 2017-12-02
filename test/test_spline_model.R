
a1=read_csv(file='test/data1.csv')

parms=parse_file(a1,logMIC=FALSE)
xobs=parms$xobs; yobs=parms$yobs; xcens=parms$xcens; ycens=parms$ycens

### Initialize
nobs=length(xobs)
xgrid=seq(min(xobs)-1,max(xobs)+1,length=1000)
parms=initialize_parms_spline(xobs,yobs,xcens,ycens,xgrid)
xtrue=parms$xtrue
coefs=parms$coefs
ytrue=as.numeric(parms$ytrue)
lowept=parms$lowept
upperept=parms$upperept
designMatrixGrid=parms$designMatrixGrid

### Run
parms=bayesian_mon_errors_spline_diag(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid,lowept,upperept,designMatrixGrid)
MICDens=parms$MICDens; fitMat=parms$fitMat; acceptCoef=parms$acceptCoef
coefMat=parms$coefMat; xtrue_sav=parms$xtrue_sav


spline_diagnostic(xgrid,xobs,yobs,xcens,ycens,MICDens,fitMat,acceptCoef,xtrue_sav,coefMat)

# parms=bayesian_mon_errors_spline(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid,lowept,upperept)


