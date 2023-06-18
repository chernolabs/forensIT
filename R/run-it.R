#' @title runIT
#' @description run information theory (IT) metrics
#' @param lped list of pedigree objects
#' @param freqs list of allele frequencies
#' @param QP QP
#' @param dbg debug
#' @param numCores number of cores
#' @param bOnlyIT boolean to only run IT
#' @param lprobg_ped list of probG
#' @param bsigma boolean to compute sigma
#' @param blog boolean to write log
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @import iterators
#' @return runIT
runIT <-function(lped=NULL,freqs,QP,dbg,numCores,bOnlyIT=FALSE,lprobg_ped=NULL,bsigma=FALSE,blog=FALSE){
  irun <- NULL
  if(TRUE){
    library(foreach)
    library(doParallel)
    registerDoParallel(cores=numCores)  
  }else{
    library(foreach)
    library(doSNOW)
  }
  if(bOnlyIT && is.null(lped)) warning('bOnlyIT is TRUE but lprobG is NULL.')
  if(blog) writeLines(c(""), "log.txt")
  t1 <- Sys.time()
  if(!bOnlyIT){
    a <- foreach(irun=1:length(lped))%dopar%{ #nolint
      #sink("log.txt", append=TRUE)
      if(blog) cat(paste("Starting irun",irun,"\n"),file='log.txt',append=TRUE)
      ped_fbnet <- convertPed(lped[[irun]])
      pbn  <- fbnet::initBN(ped_fbnet)
      bnet <- fnet::buildBN(pbn,QP=QP)
      bn1  <- fnet::buildCPTs(bnet) 
      resQ <- fnet::velim.bn(bn1,ordMethod="min_fill",verbose=FALSE)
      lprobG <- genotypeProbTable_FM(bn1,resQ, freq = freqs)
      lprobG <- lprobG$lprobG
      bnet_pop  <- compareBnetPopGenoPDFs(lprobG)
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
    a<-lapply(lprobg_ped[1:11],function(x){ #info_out  <- compareBnetPopGenoPDFs(x)
    })
  }
  t2 <- Sys.time()
  cat(length(info_ped),'runs in',t2-t1,ifelse(bOnlyIT,'(onlyIT)\n','(fbnet & IT)\n'))
  return(list(info=info_ped,lprobg=lprobg_ped))
}
