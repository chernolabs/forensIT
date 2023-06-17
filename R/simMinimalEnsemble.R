simMinimalEnsemble <- function(ped,QP,testID,freqs,numCores=1,seed=123457,bVerbose=TRUE,bJustGetNumber=FALSE,bdbg=FALSE){
  
  # LangeGoradia para ve cuantos genotipos posibles
  lLangeGoradia <- list()
  numGeno<-c()
  markerNames <- unlist(lapply(ped$MARKERS,function(x){attr(x,'name')}))
  
  for(imarker in seq_along(markerNames)){
    alleles <- attr(ped$MARKERS[[imarker]],'alleles')
    
    a<-elimLangeGoradia(ped,iMarker=imarker,bitera=TRUE,bverbose = FALSE)[testID]
    #desarmo el alelo 666
    a<-lapply(a,function(x){
      #desarmo 666's
      ma <- ma0 <- t(apply(matrix(as.numeric(strsplit2(x,'/')),byrow=FALSE,ncol=2),1,sort))
      if(any(ma[,1]=='666')){  #si esto es verdad tengo 666/666 asi que no hay ninguna restriccion
        maux <- expand.grid(alleles,alleles)
        return(unique(apply(maux,1,function(x){paste(sort(x),collapse='/')})))
      }else{
        i666 <- grep('666',ma[,2])
        if(length(i666)>0){
          for(ii in i666){
            mnew <- cbind(rep(ma0[ii,1],length(alleles)),as.numeric(alleles))
            mnew <- t(apply(mnew,1,sort))
            ma <- rbind(ma,mnew)
          }
          ma <- ma[-i666,]
        }
      }
      return(unique(apply(ma,1,paste,collapse='/')))
    })
    
    
    lLangeGoradia[[markerNames[imarker]]]<-a
    numGeno <- rbind(numGeno,unlist(lapply(a,length)))
  }
  rownames(numGeno) <- markerNames
  
  nruns <- apply(numGeno,2,max)
  if(bVerbose | bJustGetNumber){
    cat(paste('ID:',testID,'Number of fbnet runs:',nruns,'  (',rownames(numGeno)[which.max(numGeno[,1])],')\n'))
    cat('\nGenotypes per marker to explore:\n')
    print(numGeno)
  }
  if(bJustGetNumber){
    return()
  }
  
  
  
  # convierto la lista de LangeGoradia a matriz de genotipos para correr
  maxNumGenotypes <- apply(numGeno,2,max)
  lMatrixGenotype <- list()
  for(ip in seq_along(maxNumGenotypes)){
    maux <- matrix(NA,ncol=maxNumGenotypes[ip],nrow=length(lLangeGoradia))
    rownames(maux) <- names(lLangeGoradia)
    for(i in 1:nrow(maux)){
      jmax  <- numGeno[i,ip]
      x     <- lLangeGoradia[[i]][[ip]]
      maux[i,] <- c(x,rep(x[jmax],ncol(maux)-jmax))  
    }
    lMatrixGenotype[[ip]] <- maux
  }
  dim(lMatrixGenotype[[1]])
  
  
  
  # calculo de probs para diferentes genotipos 
  
  if(numCores>1){
    library(doParallel)
    registerDoParallel(cores=numCores)
  }
  
  lprobG0 <- list()
  for(inew in seq_along(testID)){
    cat('running testID:',testID[inew],'\n')
    pednew <- ped
    lprobG0[[as.character(testID[inew])]] <- list()
    
    if(bdbg) cat('', file=paste0("mylog.",inew,".txt"), append=FALSE)
    if(numCores==1){
      #set genotype for i-newopeople  
      for(irun in 1:ncol(lMatrixGenotype[[inew]])){
        cat(paste0('\n contributor: ',inew,'/',length(testID),' - ',irun,'/',nruns[inew],' runs.\n'))
        genos <- lMatrixGenotype[[inew]][,irun]
        #cargo genotipos en markerdata
        pedMarkers <- unlist(lapply(pednew$MARKERS,function(x){attr(x,'name')}))
        for(i in seq_along(genos)){
          imarker  <- match(names(genos)[i],pedMarkers)
          moi      <- pednew$MARKERS[[imarker]]
          alleles  <- attr(moi,'alleles')
          ialleles <- match(unlist(strsplit(genos[pedMarkers[imarker]],'/')),alleles)
          moi[testID[inew],] <- ialleles
          pednew$MARKERS[[imarker]] <-moi
        }
        pednew$available <- sort(c(pednew$available,as.numeric(testID[inew])))
        
        ped_fbnet <- convertPed(pednew)
        pbn  <- initBN(ped_fbnet)
        bnet <- buildBN(pbn,QP=QP)
        bn1  <- buildCPTs(bnet) 
        resQ <- velim.bn(bn1,ordMethod="min_fill",verbose=FALSE)
        lprobG0[[as.character(testID[inew])]][[paste0('sample_',irun)]] <- genotypeProbTable_FM(resQ,freq = freqs)[[1]]
      }
    }else{
      a <- foreach(irun=1:ncol(lMatrixGenotype[[inew]]))%dopar%{
        
        genos <- lMatrixGenotype[[inew]][,irun]
        #cargo genotipos en markerdata
        pedMarkers <- unlist(lapply(pednew$MARKERS,function(x){attr(x,'name')}))
        for(i in seq_along(genos)){
          imarker  <- match(names(genos)[i],pedMarkers)
          moi      <- pednew$MARKERS[[imarker]]
          alleles  <- attr(moi,'alleles')
          ialleles <- match(unlist(strsplit(genos[pedMarkers[imarker]],'/')),alleles)
          moi[testID[inew],] <- ialleles
          pednew$MARKERS[[imarker]] <-moi
        }
        pednew$available <- sort(c(pednew$available,as.numeric(testID[inew])))
        
        pped_fbnet <- convertPed(pednew)
        ppbn  <- initBN(pped_fbnet)
        pbnet <- buildBN(ppbn,QP=QP)
        pbn1  <- buildCPTs(pbnet) 
        presQ <- velim.bn(pbn1,ordMethod="min_fill",verbose=FALSE)
        
        
        res <- genotypeProbTable_FM(presQ,freq=freqs)[[1]]
        if(FALSE) cat(res[[10]],file=paste0('run_',irun,'.csv'))
        if(bdbg)cat(irun,res[[10]][1,2],'\n',file=paste0("mylog.",inew,".txt"),append=TRUE)
        res
      }
      names(a)<-paste0('sample_',seq_along(a))
      lprobG0[[as.character(testID[inew])]]  <- a
    }
    #close(pb)
  }
  
    # IT metrics   
    lprobGenoMOI <- list()
    ITtable<-c()
    for(iCDI in seq_along(lprobG0)){
      lprobGenoMOI[[iCDI]] <- list()
      for(moi in names(lprobG0[[1]][[1]])){
        #solo recupero las primeras  numGeno[moi,iCDI] samples
        a<-lapply(lprobG0[[iCDI]][1:numGeno[moi,iCDI]],function(x){
          x[[moi]]
        })
        saux     <-lLangeGoradia[[moi]][[iCDI]]
        saux     <- unlist(lapply(strsplit(saux,'/'),function(y){paste(sort(as.numeric(y)),collapse='/')}))
        names(a) <- saux 
        lprobGenoMOI[[iCDI]][[moi]]<-a
        compa <- compareBnetPopGenoPDFs(a)
        ITtable <- rbind(ITtable,cbind(cdi=names(lprobG0)[iCDI],marker=moi,allele=rownames(compa),compa))
      }
    }
    #saux <- names(lprobG0)
    #saux <- unlist(lapply(strsplit(saux,'/'),function(y){paste(sort(as.numeric(y)),collapse='/')}))
    names(lprobGenoMOI)<-names(lprobG0)

  
  # lprobGenoMOI[ID][marker][genotype_realization]: data.frame with fbnet results
  # lprobG[ID][sample][marker]: data.frame with fbnet results
  # lMatrixGenotype[[ID]]: genotype dataframe
  # ITtable: IT metrics data.drame
  return(list(lprobGenoMOI=lprobGenoMOI,lprobG=lprobG0,lMatrixGenotype=lMatrixGenotype,ITtable=ITtable))
}

