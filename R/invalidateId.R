InvalidateId<-function(ped,id=NULL,markerId=1,double=FALSE,propagate=FALSE){
 if(is.null(id)) warning("se debe especificar id de individuo a mutar\n")
 
  
  # paso a variables internas...(orig ids)
  iid <- which(ped$orig.ids==id)
  
  parentIds     <- ped$pedigree[iid,c("FID","MID")]
  parentAlleles <- ped$markerdata[[markerId]][parentIds,]
  parentAlleles <- unique(as.vector(parentAlleles))
  
  u <- (1:attr(ped$markerdata[[markerId]],"nalleles"))
  u <- u[!u %in% parentAlleles]
  if(!double){
   ped$markerdata[[markerId]][iid,sample(1:2,1)]<-sample(u,1)
  }else{
   ped$markerdata[[markerId]][iid,]<-sample(u,2)
  }
 
  #propago el nuevo valor aguas abajo
  if(propagate){
    df      <- as.data.frame(ped)[,1:5]
    descIds <- df[which(df$FID==id | df$MID==id),"ID"]
          
    #simulo la rama que cambia
    a<-branch(ped,id)
    a<-modifyMarker(a,markerId,descIds,c(0,0))
    a<-markerSim(a,available=descIds,N=1,partialmarker=markerId)
    
    adf <- as.data.frame(a)
    iaDescIds <- match(descIds,adf$ID)
    iDescIds  <- match(descIds,df$ID)    
    ped$markerdata[[markerId]][iDescIds,] <- a$markerdata[[markerId]][iaDescIds,]
  
  }
  return(ped)
}

