buildEnsembleCPTs <- function(lsimu,lminimalProbGenoMOI){
  aa <- simplify2array(lsimu,higher = FALSE)
  aa <- data.frame(sample=1:nrow(lsimu[[1]]),id=rep(colnames(lsimu[[1]]),each=nrow(lsimu[[1]])),aa)

  lz        <- vector("list",length(unique(aa$id)))
  names(lz) <- unique(aa$id)
  for(i in seq_along(lz)) lz[[i]]<-list()
    
  for(irow in 1:nrow(aa)){
    bb <- list()
    for(j in 3:ncol(aa)){
      marker <- colnames(aa)[j]
      geno   <- aa[irow,j]
      
      if(is.null(lminimalProbGenoMOI[[aa[irow,'id']]][[marker]][[geno]])){cat(paste(irow),'\n')}
      
      bb     <- c(bb,list(lminimalProbGenoMOI[[aa[irow,'id']]][[marker]][[geno]]))
    }
    names(bb)  <- colnames(aa)[-c(1,2)]
    bb2        <- list(bb)
    names(bb2) <- paste0('sample_',aa[irow,'sample']) 
    lz[[aa[irow,'id']]] <- c(lz[[aa[irow,'id']]],bb2)
  }

  
  return(lz)
}
