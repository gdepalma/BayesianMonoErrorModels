

#plots
par(mfrow=c(1,2))
plot(xgrid,densTrue,type='l')
lines(xgrid,MICDens,col=2)
plot(xgrid,trueFit,type='l')
lines(xgrid,fit,col=2)
# plot(acceptCoef,main='Accept Coef')
# plot(smoothAccept,main='smoothAccept Coef')

par(mfrow=c(2,2))
plot(post,main='Posterior',type='l')
plot(smoothSplineSave,main='Smoothness Parameter',type='l')
plot(sigmaM,main='xsig',type='l')
plot(sigmaD,main='ysig',type='l')

# par(mfrow=c(1,1))
# plot(scaleSigSave,main='scaleSig',type='l')

a1=data.frame(coefMat)
a1$grid=1:numIter
a2=melt(a1,id='grid')
plt=ggplot(a2,aes(x=grid,y=value))+geom_line()+facet_wrap(~ variable,scales = "free_y")
print(plt)

### Save fit and MIC Dens within Data Range
tempXgrid=xgrid[xgrid>min(xobs) & xgrid<max(xobs)]
MICDensTemp=MICDens[xgrid>min(xobs) & xgrid<max(xobs)]
fitTemp=fit[xgrid>min(xobs) & xgrid<max(xobs)]
densTrueTemp=densTrue[xgrid>min(xobs) & xgrid<max(xobs)]
trueFitTemp=trueFit[xgrid>min(xobs) & xgrid<max(xobs)]
par(mfrow=c(1,2))
plot(tempXgrid,densTrueTemp,type='l')
lines(tempXgrid,MICDensTemp,col=2)
plot(tempXgrid,trueFitTemp,type='l')
lines(tempXgrid,fitTemp,col=2)

### Compute DIA Breakpoints
M1Test=-1; M2Test=1
weights=MICDensTemp/sum(MICDensTemp)
parms=findDIAC(yobs,tempXgrid,weights,fitTemp,M1Test,M2Test,xsig,ysig,3,20)
print(parms)


