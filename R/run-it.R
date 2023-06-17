runIT <-function(lped=NULL,freqs,QP,dbg,numCores,bOnlyIT=FALSE,lprobg_ped=NULL,bsigma=FALSE,blog=FALSE){
  if(TRUE){
    library(foreach)
    library(doParallel)
    registerDoParallel(cores=numCores)  
  }else{
    library(foreach)
    library(doSNOW)
  }
  if(bOnlyIT & is.null(lped)) warning('bOnlyIT is TRUE but lprobG is NULL.')
  if(blog) writeLines(c(""), "log.txt")
  t1 <- Sys.time()
  if(!bOnlyIT){
    a <- foreach(irun=1:length(lped))%dopar%{
      #sink("log.txt", append=TRUE)
      if(blog) cat(paste("Starting irun",irun,"\n"),file='log.txt',append=TRUE)
      ped_fbnet <- convertPed(lped[[irun]])
      pbn  <- initBN(ped_fbnet)
      bnet <- buildBN(pbn,QP=QP)
      bn1  <- buildCPTs(bnet) 
      resQ <- velim.bn(bn1,ordMethod="min_fill",verbose=FALSE)
      lprobG <- genotypeProbTable_FM(resQ, freq = freqs)
      lprobG <- lprobG$lprobG
      bnet_pop  <- compareBnetPopGenoPDFs(lprobG,dbg=dbg,bsigma=bsigma)
      if(blog) cat(paste("Ending irun",irun,"\n"),file='log.txt',append=TRUE)
      
      list(info_out=bnet_pop,lprobG_out=lprobG)
    }
    #sink(file=NULL)
    #Rearreglo listas
    # paso de [sample][info|probG][marker] a [info|probG][sample][marker]
    info_ped <- lprobg_ped <- list()
    for(isample in seq_along(a)){
      info_ped[[isample]]  <- a[[isample]][['info_out']]  
      lprobg_ped[[isample]] <- a[[isample]][['lprobG_out']] 
    }
  }else{
    info_ped
    a<-lapply(lprobg_ped[1:11],function(x){
      info_out  <- compareBnetPopGenoPDFs(x,dbg=dbg,bsigma=TRUE)
    })
  }
  
  t2 <- Sys.time()
  cat(length(info_ped),'runs in',t2-t1,ifelse(bOnlyIT,'(onlyIT)\n','(fbnet & IT)\n'))
  return(list(info=info_ped,lprobg=lprobg_ped))
}
