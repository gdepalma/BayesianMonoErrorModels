diagnostics=function(xgrid,xobs,yobs,xcens,ycens,acceptCoef,MICDens,fitMat,coefMat){
  
  # ### Posteriors
  MICDens_1=apply(MICDens,2,median)
  MICDens_1=MICDens_1/sum(MICDens_1)
  fitMat_1=apply(fitMat,2,median)
  
  plot(xgrid,MICDens_1)
  plot(xgrid,fitMat_1)
  
  
  
  print(mean(acceptCoef))
  a1=melt(data.frame(coefMat))
  a1$iter=rep(1:nrow(coefMat),ncol(coefMat))
  ggplot(a1,aes(x=iter,y=value))+geom_line()+facet_wrap(~variable,scales='free_y')
  
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


  ### MIC Density
  densDat=data_frame(xgrid,y=apply(MICDens,2,median))
  lower_dens=data_frame(xgrid,y=apply(MICDens,2,function(x) quantile(x,probs=c(.025))))
  upper_dens=data_frame(xgrid,y=apply(MICDens,2,function(x) quantile(x,probs=c(.975))))

  pltDens=ggplot(densDat,aes(x=xgrid,y))+geom_line()+
    geom_line(data=lower_dens,aes(x=xgrid,y=y,label=NULL),linetype=2,color='deepskyblue4')+
    geom_line(data=upper_dens,aes(x=xgrid,y=y,label=NULL),linetype=2,color='deepskyblue4')+
    scale_x_continuous(breaks = seq(min(xgrid),max(xgrid),by=1))+
    theme_fivethirtyeight()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          plot.title=element_text(size=15))+
    labs(title='',y='',x=expression(MIC~(log["2"]~ug/mL)))
  
  ### MIC/DIA Relationship
  fitMAPDAT=data_frame(xgrid,gx_median=apply(fitMat,2,median))
  fit025DAT=data_frame(xgrid,gx_lower=apply(fitMat,2,function(x) quantile(x,probs=c(.025))))
  fit975DAT=data_frame(xgrid,gx_upper=apply(fitMat,2,function(x) quantile(x,probs=c(.975))))
  
  pltRel=ggplot(a1,aes(x=xobs,y=yobs,label=Freq))+geom_text(size=3.2,color='black')+
    geom_line(data=fitMAPDAT,aes(x=xgrid,y=gx_median,label=NULL),color='deepskyblue4')+
    geom_line(data=fit025DAT,aes(x=xgrid,y=gx_lower,label=NULL),linetype=2,color='deepskyblue4')+
    geom_line(data=fit975DAT,aes(x=xgrid,y=gx_upper,label=NULL),linetype=2,color='deepskyblue4')+
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

  ### Combine Fits
  plt1 <- ggplot_gtable(ggplot_build(pltRel))
  plt2 <- ggplot_gtable(ggplot_build(pltDens))
  maxWidth = unit.pmax(plt1$widths[2:3], plt2$widths[2:3])
  plt1$widths[2:3] <- maxWidth
  plt2$widths[2:3] <- maxWidth
  plot(grid.arrange(plt1, plt2, ncol=1, heights=3:2))
  
  
  par(mfrow=c(2,1))
  xtrue_mean=apply(xtrue_sav,2,mean)
  hist(xtrue_mean,freq=FALSE)
  hist(xobs,freq=FALSE)

  return(list(plt1=plt1,plt2=plt2))
  
}
