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
  
  xobs=ifelse(logMIC==TRUE,log(as.numeric(xobs),2),as.numeric(xobs))
  yobs=as.numeric(yobs)
  
  return(list(xobs=xobs,yobs=yobs,xcens=xcens,ycens=ycens))
  
}