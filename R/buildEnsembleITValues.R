buildEnsembleITValues <-function(lsimu=lsimulation,ITtab=sim$ITtable,bFullIT=FALSE){
  # if(is.character(marker)){
  #   imrkr<-which(names(lsimu)==marker)
  # }else{
  #   imrkr<-marker 
  # }
  # 
  # if(is.character(cdi)){
  #   icdi <- which(colnames(lsimu[[1]])%in%cdi)
  # }else{
  #   icdi <- cdi
  # }    
  
  #res <- 
  fullIT<-ensembleIT<-c()
  for(icdi in 1:ncol(lsimu[[1]])){
    newp <- colnames(lsimu[[1]])[icdi]
    #resb <- c()
    for(imrkr in seq_along(lsimu)){
      tt  <- ITtab[ITtab$marker==names(lsimu)[imrkr],]
      itt <- match(lsimu[[imrkr]][,icdi],tt[tt$cdi==newp,'allele'])
      # a   <- t(apply(tt[tt$cdi==newp,][itt,4:8],2,function(x){
      #   c(sum(x),sum(x^2))
      # }))
      # colnames(a) <- paste0(names(lsimu)[imrkr],c('','__2'))
      # resa <- cbind(resa,a)
      
      n       <- nrow(lsimu[[imrkr]])
      aa      <- data.frame(sample=seq_along(itt),id=newp,marker=names(lsimu)[imrkr],tt[tt$cdi==newp,][itt,4:ncol(tt)])
      fullIT  <- rbind(fullIT,aa)
      
      #b <- apply(tt[tt$cdi==newp,][itt,4:ncol(tt)],2,function(x){c(sum(x)/n)})
      #resb <- rbind(resb,b)
    }
    #resb[is.nan(resb)]<-0
    #rownames(resb) <- names(lsimu)

    aa <- aggregate(fullIT[,4:ncol(fullIT)],
                    by = list(sample=fullIT$sample,id=fullIT$id),
                    function(x){sum(x)})
    #Exclusion power
    ep <- aggregate(fullIT[,'pexclusion'],
                    by = list(sample=fullIT$sample,id=fullIT$id),
                    function(x){1-prod(1-(x))})[,3]
    
    aa               <- cbind(aa[,c(1:8,11:12)],exclusionPower=ep)
    colnames(aa)[9]  <- 'numX'
    aa[,10]          <- aa[,10]/length(lsimu)
    ensembleIT       <- rbind(ensembleIT,aa)
  }
  if(bFullIT){
    return(list(ensembleIT=ensembleIT,fullIT=fullIT))
  }else{
    return(list(ensembleIT=ensembleIT))
  }
}
