simTestIDMarkers<-function(ped,testID,numSim=10,seed=123457){ 
  set.seed(seed)
  markerNames <- unlist(lapply(ped$MARKERS,function(x){attr(x,'name')}))
  # .... lsimulation 
  ipeople <- seq_along(testID)#1:2
  lsimulation<-list()
  for(imarker in seq_along(markerNames)){
    # (a<-markerSim(ped,N=numSim,partialmarker=imarker,
    #               available = testID[ipeople],verbose = FALSE,seed=seed))
    (a<-forrel::markerSim(ped,N=numSim,partialmarker=imarker,ids = testID[ipeople],verbose = FALSE))
    laux <- lapply(a$MARKERS,function(x){
      xx<-attr(x,'alleles')[as.vector(t(x[testID,]))]})
    #xx<-attr(x,'alleles')[x[testID[ipeople],]]})
    
    laux <- lapply(laux,function(x){
      ma <-apply(matrix(x,byrow=TRUE,nrow=length(ipeople)),1,function(y){
        paste(sort(as.numeric(y)),collapse='/')})
      return(ma)
    })
    maa <-matrix(unlist(laux),byrow=TRUE,ncol=length(ipeople))
    colnames(maa)<-testID[ipeople]
    lsimulation[[markerNames[imarker]]]<-maa
  }
  return(lsimulation)
}
if(FALSE){
  # (6.2) Reconstrucciuon de distribuciones IT para simulacion 
  # Funcion para obtener on-the-fly los IT values del marcador 'marker' 
  # del ensemble de simulaciones corridas para los posibles genotipos
  # del nuevo contribuyente 'cdi'
  getMarkerITsimValues <-function(marker,cdi,lsimu=lsimulation,
                                  ITtab=ITtable,newp=testID){
    if(is.character(marker)){
      imrkr<-which(names(lsimu)==marker)
    }else{
      imrkr<-marker 
    }
    if(is.character(cdi)){
      icdi <- which(colnames(lsimu[[1]])%in%cdi)
    }else{
      icdi <- cdi
    }    
    tt  <- ITtab[ITtab$marker==names(lsimu)[imrkr],]
    itt <- match(lsimu[[imrkr]][,icdi],tt[tt$cdi==newp[icdi],'allele'])
    return(tt[tt$cdi==newp[icdi],][itt,])
  }
  a  <- getMarkerITsimValues(marker=3,cdi=1)
  aa <- melt(a,id.vars=c('cdi','marker','allele'))
  ggplot(aa,aes(x=value))
}
