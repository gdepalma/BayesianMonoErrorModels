
  ### Likelihoods

a1=data_frame(iter=1:length(llk),llk,popllk,prior,post) %>%
    gather(type,val,-iter)
ggplot(a1,aes(x=iter,y=val))+geom_line()+facet_wrap(~type,scales='free_y')+
  labs(y='',x='Iteration')

### Accept rate for coefficients update
# plot(acceptCoef)
# mean(acceptCoef)

plot(kSave)

### Coefficients

coefMatSave1=as.data.frame(coefMatSave) %>% mutate(iter=1:nrow(coefMatSave)) %>%
  gather(coef,val,-iter)
ggplot(coefMatSave1,aes(x=iter,y=val))+geom_line()+facet_wrap(~coef,scales='free_y')+
  labs(y='',x='Iteration')

### xtrue
a1=data.frame(xobs,xtrue) %>% gather(type,val)
ggplot(data=a1, aes(val,..density..)) + facet_wrap(~type,ncol=1)+
  geom_histogram(binwidth = 0.5)+labs(x='')


#### Fits
xobs1=xobs
yobs1=yobs
xobs[xcens==1 & xobs==max(xobs)]=max(xobs)+1
xobs[xcens==-1 & xobs==min(xobs)]=min(xobs)-1
yobs[ycens==1 & yobs==max(yobs)]=max(yobs)+1
yobs[ycens==-1 & yobs==min(yobs)]=min(yobs)-1
a1=data.frame(table(xobs,yobs))
a1$xobs=as.numeric(as.character(a1$xobs))
a1$yobs=as.numeric(as.character(a1$yobs))
a1=a1[a1$Freq>0,]

### Get posterior Logistic Fit

# thin
coefMatSave=coefMatSave[seq(1,nrow(coefMatSave),by=4),]
MICDens=MICDens[seq(1,nrow(MICDens),by=4),]
fitMat=fitMat[seq(1,nrow(fitMat),by=4),]


densTrueDat=data.frame(xgrid,densTrue=densTrue*1000)
fitTrueDat=data.frame(xgrid,trueFit=as.numeric(trueFit))
# fit025DAT=data.frame(xgrid,gx_lower=getYtrue(apply(coefMatSave,2,function(x) quantile(x,probs=0.025)),xgrid))
# fitMAPDAT=data_frame(xgrid,gx_median=getYtrue(apply(coefMatSave,2,median),xgrid))
# fit975DAT=data.frame(xgrid,gx_upper=getYtrue(apply(coefMatSave,2,function(x) quantile(x,probs=0.975)),xgrid))
MICDensDat=data.frame(xgrid,densMedian=apply(MICDens,2,median)*1000)
fit025DAT=data.frame(xgrid,gx_lower=apply(fitMat,2,function(x) quantile(x,probs=0.025)))
fitMAPDAT=data_frame(xgrid,gx_median=apply(fitMat,2,median))
fit975DAT=data.frame(xgrid,gx_upper=apply(fitMat,2,function(x) quantile(x,probs=0.975)))

ggplot(a1,aes(x=xobs,y=yobs,label=Freq))+geom_text(size=3.2,color='black')+
  geom_line(data=fitMAPDAT,aes(x=xgrid,y=gx_median,label=NULL),color='deepskyblue4')+
  geom_line(data=fit025DAT,aes(x=xgrid,y=gx_lower,label=NULL),linetype=2,color='deepskyblue4')+
  geom_line(data=fit975DAT,aes(x=xgrid,y=gx_upper,label=NULL),linetype=2,color='deepskyblue4')+
  geom_line(data=MICDensDat,aes(x=xgrid,y=densMedian,label=NULL),linetype=2)+
  geom_line(data=densTrueDat,aes(x=xgrid,y=densTrue,label=NULL),linetype=1)+
  geom_line(data=fitTrueDat,aes(x=xgrid,y=trueFit,label=NULL),linetype=1)+
  scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                     labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                     limits = c(min(xobs1)-1,max(xobs1)+1))+
  scale_y_continuous(breaks = seq(min(yobs1)-1,max(yobs1)+1,by=1),
                     labels = c(paste("<",min(yobs1),sep=''),seq(min(yobs1),max(yobs1),by=1), paste(">",max(yobs1),sep='')),
                     limits = c(min(yobs1)-1,max(yobs1)+1))+
  theme_fivethirtyeight()+
  theme(axis.line = element_line(colour = "black"),
        legend.position='none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=8),
        axis.title=element_text(size=11),
        plot.title=element_text(size=15))+
  labs(title='Logistic Model',y='DIA (mm)',x="")

