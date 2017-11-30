
a1=read_csv(file='../data1.csv')

parms=parse_file(a1,FALSE)
xobs=parms$xobs; yobs=parms$yobs; xcens=parms$xcens; ycens=parms$ycens

### Initialize
nobs=length(xobs)
xgrid=seq(min(xobs)-1,max(xobs)+1,length=1000)
parms=initialize_parms_logistic(xobs,yobs,xcens)
xtrue=parms$xtrue
coefs=parms$coefs
ytrue=getylogtrue1(coefs,xtrue)


### Run
parms=bayesian_mon_errors_logistic(xobs,yobs,xcens,ycens,coefs,xtrue,ytrue,xgrid)

### plots
  