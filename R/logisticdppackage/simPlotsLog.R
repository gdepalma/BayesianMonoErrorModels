library(gridExtra)
library(grid)
library(ggthemes)
library(dplyr)
library(ggplot2)
library(tidyr)
# dyn.load('splineFunctionsC.so')
source(paste(getwd(),'/GenData.R',sep=''))


setwd(paste0(getwd(),'/Runs'))
load('LogisticGen2_four.Rdata')
MICDensEst=tmp$saveMICDensEst
FitEst=tmp$saveFitEst

xgrid=seq(-12,12,length=1000)


getTrue=function(type=5){
  #### Type 1 True
  if(type==1){
    popmn=c(-6,3); popstd=c(2,.7); popprob=c(.8,.2)
    fitTrue=30-2*xgrid 
    
    # Type 2
  }else if(type==2){
    popmn=c(-4,0,4); popstd=c(.7,.7,.7); popprob=c(1/3,1/3,1/3)
    coef=c(29,1.17,.48)
    fitTrue=coef[1]*(exp(coef[2]-coef[3]*xgrid)/(1+exp(coef[2]-coef[3]*xgrid)))
  
  
  ### Type 3
  }else if(type==3){
    popmn=c(-4,0,4); popstd=c(.7,.7,.7); popprob=c(1/3,1/3,1/3)
    icoefsT=c(1,1,20,1,20,1)
    parms=Ispline(c(-3,0,1),-6,6)
    basesT=parms$bases;
    knotseqT=parms$knotseq
    designMatrix=getIsplineC(xgrid,knotseqT,basesT)
    fitTrue = icoefsT%*%t(designMatrix)+5  
    
  }else if(type==4){
    popmn=c(-4,0,4); popstd=c(3,3,3); popprob=c(1/3,1/3,1/3)
    icoefsT=c(1,10,1,25,1,1,10,1)
    parms=Ispline(c(-4,-2,0,2,4),-6,6)
    basesT=parms$bases;
    knotseqT=parms$knotseq
    designMatrix=getIsplineC(xgrid,knotseqT,basesT)
    fitTrue = icoefsT%*%t(designMatrix)
  
  ### Type 5
  }else if(type==5){
    popmn=c(-4.6,-2,1); popstd=c(1.1,1.5,1.5); popprob=c(.6,.2,.2)
    coef=c(35,1.17,.1,1.2)
    mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
    fx = 1/(1+exp(-mb*(coef[2]-xgrid)))
    fitTrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*exp(coef[4]*
                 (coef[2]-xgrid)))/(1+fx*exp(coef[3]*(coef[2]-xgrid))+(1-fx)*
                  exp(coef[4]*(coef[2]-xgrid)))
  }
  
  
  densTrue=rep(0,length(xgrid))
  for(i in 1:length(xgrid))
    densTrue[i]=sum(popprob*dnorm(xgrid[i],popmn,popstd))
  densTrue=densTrue/sum(densTrue)
  
  
  return(list(fitTrue=fitTrue,densTrue=densTrue))
  
}

trueFits=getTrue(2)


### Don't include all data
point1=-6
point2=6
points=which(xgrid>=point1 & xgrid<=point2)
xgrid=xgrid[points]

tmpMIC=matrix(0,nrow=200,ncol=length(points))
tmpFit=matrix(0,nrow=200,ncol=length(points))
for(i in 1:200){
  tmpMIC[i,]=MICDensEst[i,points]
  tmpFit[i,]=FitEst[i,points]
}
MICDensEst=tmpMIC
FitEst=tmpFit

### Normalize MIC Dens Est
MICDensEst=apply(MICDensEst,1,function(x) x/sum(x))

### MIC plot
a1=data.frame(cbind(xgrid,MICDensEst))
a1=gather(data = a1, key = variable, value = value, -xgrid)
a2=data.frame(xgrid,densTrue=trueFits$densTrue[points])
plt1=ggplot(a1,aes(x=xgrid,y=value,group=variable))+
  geom_line(color='red',alpha=.1)+
  scale_x_continuous(breaks=seq(-12,12,by=2))+
  geom_line(data=a2,aes(x=xgrid,y=densTrue,group=1),color='black',size=.6)+
  theme_fivethirtyeight()+
  xlab(expression(paste('MIC (',log[2],' ug/mL)')))+
  labs(y='Density')+
  theme(axis.title = element_text(),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=11),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=11))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

  
                             

### Fit plot
a1=data.frame(cbind(xgrid,t(FitEst)))
a1=gather(data = a1, key = variable, value = value,-xgrid)
a2=data.frame(xgrid,fitTrue=trueFits$fitTrue[points])
plt2=ggplot(a1,aes(x=xgrid,y=value,group=variable))+
  geom_line(color='red',alpha=.1)+
  scale_x_continuous(breaks=seq(-12,12,by=2))+
  geom_line(data=a2,aes(x=xgrid,y=fitTrue,group=1),color='black',size=.6)+
  theme_fivethirtyeight()+
  labs(x='',y='DIA (mm)',title='Four-parameter Logistic')+
  theme(axis.title = element_text(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=11))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

### Combine Fits
plt1 <- ggplot_gtable(ggplot_build(plt1))
plt2 <- ggplot_gtable(ggplot_build(plt2))
maxWidth = unit.pmax(plt1$widths[2:3], plt2$widths[2:3])
plt1$widths[2:3] <- maxWidth
plt2$widths[2:3] <- maxWidth

png(file = "fit2_4.png", width = 8, height = 6, units = 'in',res=300)
grid.arrange(plt2, plt1, ncol=1, heights=3:2)
dev.off()

#### Fit Statistics
fitTrue=trueFits$fitTrue[points]
densTrue=trueFits$densTrue[points]
densTrue=densTrue/sum(densTrue)*100000

getStatistics=function(DensEst,FitEst){

  sseDens=0
  sseFit=0
  for(i in length(xgrid)){
    sseDens=sseDens+(densTrue[i]-DensEst[i])^2
    sseFit=sseFit+(fitTrue[i]-FitEst[i])^2
  }
  return(list(sseDens=sseDens,sseFit=sseFit))
}

sseDens=rep(0,200)
sseFit=rep(0,200)
for(i in 1:200){
  parms=getStatistics(MICDensEst[,i]*100000,FitEst[i,])
  sseDens[i]=parms$sseDens
  sseFit[i]=parms$sseFit
}

round(mean(sseDens),2); round(median(sseDens),2); round(sd(sseDens),2)
round(mean(sseFit),2); round(median(sseFit),2); round(sd(sseFit),2)
# 
# 
# 
