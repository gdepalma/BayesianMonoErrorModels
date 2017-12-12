output_graph_one_model_twoMIC_test=function(list_of_draws,xgrid,dat_sav,MICBrkptL,MICBrkptU){


  xobs=dat_sav$xobs
  yobs=dat_sav$yobs
  xcensl=dat_sav$xcensl
  xcensu=dat_sav$xcensu
  ycensl=dat_sav$ycensl
  ycensu=dat_sav$ycensu

  xobs1=xobs
  yobs1=yobs
  xobs[xcensu==1 & xobs==max(xobs)]=max(xobs)+1
  xobs[xcensl==1 & xobs==min(xobs)]=min(xobs)-1
  yobs[ycensu==1 & yobs==max(yobs)]=max(yobs)+1
  yobs[ycensl==1 & yobs==min(yobs)]=min(yobs)-1
  a1=data.frame(table(xobs,yobs))
  a1$xobs=as.numeric(as.character(a1$xobs))
  a1$yobs=as.numeric(as.character(a1$yobs))
  a1=a1[a1$Freq>0,]



  ### MIC Density
  MIC_Dens=list_of_draws$MIC_Dens
  densDat=data_frame(xgrid,y=apply(MIC_Dens,2,mean))
  densDat$lower=apply(MIC_Dens,2,function(x) quantile(x,probs=c(.05)))
  densDat$upper=apply(MIC_Dens,2,function(x) quantile(x,probs=c(.95)))


  pltDens=ggplot(densDat,aes(x=xgrid,y))+geom_line(color='darkblue')+
    geom_ribbon(aes(ymin=lower, ymax=upper),alpha=.5,fill='cyan4')+
    geom_vline(xintercept=MICBrkptL+.5,lty=2,alpha=.5)+
    geom_vline(xintercept=MICBrkptU-.5,lty=2,alpha=.5)+
    scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
           labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
           limits = c(min(xobs1)-1,max(xobs1)+1))+
    theme_fivethirtyeight()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(color='gray90'),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          plot.title=element_text(size=15))+
    labs(title='',y='',x=expression(MIC~(log["2"]~ug/mL)))

  ### MIC/DIA Relationship
  gx=list_of_draws$gx
  gxDat=data_frame(xgrid,y=apply(gx,2,mean))
  gxDat$lower=apply(gx,2,function(x) quantile(x,probs=c(.05)))
  gxDat$upper=apply(gx,2,function(x) quantile(x,probs=c(.95)))


  pltRel=ggplot(data=a1,aes(x=xobs,y=yobs,label=Freq))+geom_text(size=3.2,color='black')+
    geom_line(data=gxDat,aes(x=xgrid,y=y),color='darkblue',inherit.aes = FALSE)+
    geom_ribbon(data=gxDat,aes(x=xgrid,ymin=lower, ymax=upper),fill='cyan4',alpha=.5,inherit.aes = FALSE)+
    geom_vline(xintercept=MICBrkptL+.5,lty=2,alpha=.5)+
    geom_vline(xintercept=MICBrkptU-.5,lty=2,alpha=.5)+
    scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
            labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
            limits = c(min(xobs1)-1,max(xobs1)+1))+
    scale_y_continuous(breaks = seq(min(yobs1)-1,max(yobs1)+1,by=1),
            labels = c(paste("<",min(yobs1),sep=''),seq(min(yobs1),max(yobs1),by=1), paste(">",max(yobs1),sep='')),
            limits = c(min(yobs1)-1,max(yobs1)+1))+
    theme_fivethirtyeight()+
    theme(axis.line = element_line(colour = "black"),
          legend.position='none',
          panel.grid.major = element_line(color='gray90'),
          panel.grid.minor = element_blank(),
          axis.text.x=element_text(size=11),
          axis.text.y=element_text(size=8),
          axis.title=element_text(size=11),
          plot.title=element_text(size=15))+
    labs(title='Logistic Model',y='DIA (mm)',x="")


  plt1 <- ggplot_gtable(ggplot_build(pltRel))
  plt2 <- ggplot_gtable(ggplot_build(pltDens))
  maxWidth = unit.pmax(plt1$widths[2:3], plt2$widths[2:3])
  plt1$widths[2:3] <- maxWidth
  plt2$widths[2:3] <- maxWidth
  plot(grid.arrange(plt1, plt2, ncol=1, heights=c(5,2)))

  return(invisible())

}

output_graph_one_model_oneMIC_test=function(list_of_draws,xgrid,dat_sav,MICBrkpt){


  xobs=dat_sav$xobs
  yobs=dat_sav$yobs
  xcensl=dat_sav$xcensl
  xcensu=dat_sav$xcensu
  ycensl=dat_sav$ycensl
  ycensu=dat_sav$ycensu

  xobs1=xobs
  yobs1=yobs
  xobs[xcensu==1 & xobs==max(xobs)]=max(xobs)+1
  xobs[xcensl==1 & xobs==min(xobs)]=min(xobs)-1
  yobs[ycensu==1 & yobs==max(yobs)]=max(yobs)+1
  yobs[ycensl==1 & yobs==min(yobs)]=min(yobs)-1
  a1=data.frame(table(xobs,yobs))
  a1$xobs=as.numeric(as.character(a1$xobs))
  a1$yobs=as.numeric(as.character(a1$yobs))
  a1=a1[a1$Freq>0,]



  ### MIC Density
  MIC_Dens=list_of_draws$MIC_Dens
  densDat=data_frame(xgrid,y=apply(MIC_Dens,2,mean))
  densDat$lower=apply(MIC_Dens,2,function(x) quantile(x,probs=c(.05)))
  densDat$upper=apply(MIC_Dens,2,function(x) quantile(x,probs=c(.95)))


  pltDens=ggplot(densDat,aes(x=xgrid,y))+geom_line(color='darkblue')+
    geom_ribbon(aes(ymin=lower, ymax=upper),alpha=.5,fill='cyan4')+
    geom_vline(xintercept=MICBrkpt+.5,lty=2,alpha=.5)+
    scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                       labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                       limits = c(min(xobs1)-1,max(xobs1)+1))+
    theme_fivethirtyeight()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(color='gray90'),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          plot.title=element_text(size=15))+
    labs(title='',y='',x=expression(MIC~(log["2"]~ug/mL)))

  ### MIC/DIA Relationship
  gx=list_of_draws$gx
  gxDat=data_frame(xgrid,y=apply(gx,2,mean))
  gxDat$lower=apply(gx,2,function(x) quantile(x,probs=c(.05)))
  gxDat$upper=apply(gx,2,function(x) quantile(x,probs=c(.95)))


  pltRel=ggplot(data=a1,aes(x=xobs,y=yobs,label=Freq))+geom_text(size=3.2,color='black')+
    geom_line(data=gxDat,aes(x=xgrid,y=y),color='darkblue',inherit.aes = FALSE)+
    geom_ribbon(data=gxDat,aes(x=xgrid,ymin=lower, ymax=upper),fill='cyan4',alpha=.5,inherit.aes = FALSE)+
    geom_vline(xintercept=MICBrkpt+.5,lty=2,alpha=.5)+
    scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                       labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                       limits = c(min(xobs1)-1,max(xobs1)+1))+
    scale_y_continuous(breaks = seq(min(yobs1)-1,max(yobs1)+1,by=1),
                       labels = c(paste("<",min(yobs1),sep=''),seq(min(yobs1),max(yobs1),by=1), paste(">",max(yobs1),sep='')),
                       limits = c(min(yobs1)-1,max(yobs1)+1))+
    theme_fivethirtyeight()+
    theme(axis.line = element_line(colour = "black"),
          legend.position='none',
          panel.grid.major = element_line(color='gray90'),
          panel.grid.minor = element_blank(),
          axis.text.x=element_text(size=11),
          axis.text.y=element_text(size=8),
          axis.title=element_text(size=11),
          plot.title=element_text(size=15))+
    labs(title='Logistic Model',y='DIA (mm)',x="")


  plt1 <- ggplot_gtable(ggplot_build(pltRel))
  plt2 <- ggplot_gtable(ggplot_build(pltDens))
  maxWidth = unit.pmax(plt1$widths[2:3], plt2$widths[2:3])
  plt1$widths[2:3] <- maxWidth
  plt2$widths[2:3] <- maxWidth
  plot(grid.arrange(plt1, plt2, ncol=1, heights=c(5,2)))

  return(invisible())

}
