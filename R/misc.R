theme_dbets <- function () {
  theme_fivethirtyeight() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(color='gray90'),
          panel.grid.minor = element_blank(),
          axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          plot.title=element_text(size=15,hjust = 0))
}

parse_file=function(a1,logMIC=FALSE){

  xcens=rep(0,nrow(a1))
  ycens=rep(0,nrow(a1))
  xobs=rep(NA,nrow(a1))
  yobs=rep(NA,nrow(a1))

  for(i in 1:nrow(a1)){
    if(grepl('>',a1[i,1])==TRUE){
      xobs[i]=substr(a1[i,1],2,nchar(a1[i,1]))
      xcens[i]=1
    }else if(grepl('<',a1[i,1])==TRUE){
      xobs[i]=substr(a1[i,1],2,nchar(a1[i,1]))
      xcens[i]=-1
    }else
      xobs[i]=a1[i,1]

    if(grepl('>',a1[i,2])==TRUE){
      yobs[i]=substr(a1[i,2],2,nchar(a1[i,2]))
      ycens[i]=1
    }else if(grepl('<',a1[i,2])==TRUE){
      yobs[i]=substr(a1[i,2],2,nchar(a1[i,2]))
      ycens[i]=-1
    }else
      yobs[i]=a1[i,2]
  }

  xobs=as.numeric(xobs)
  yobs=as.numeric(yobs)
  if(logMIC==TRUE) xobs = log(xobs,2)

  return(list(xobs=xobs,yobs=yobs,xcens=xcens,ycens=ycens))

}
