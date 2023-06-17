compareBnetPopGenoPDFs <- function(lprobTable){
  klDiv <- resCrossH<- dH <- dHnorm <- exP1 <- exP2 <- c() 
  bsigma <- FALSE
  dbg    <- TRUE
  for(i in seq_along(lprobTable)){
    resCrossH <- rbind(resCrossH,c(crossH(lprobTable[[i]][,"pop"],lprobTable[[i]][,"bnet"]),
                                  crossH(lprobTable[[i]][,"bnet"],lprobTable[[i]][,"pop"])))
    dH        <- rbind(dH,H(lprobTable[[i]][,'pop'])-H(lprobTable[[i]]['bnet']))
    dHnorm    <- rbind(dHnorm,H(lprobTable[[i]][,'pop'],normalized = TRUE)-H(lprobTable[[i]]['bnet'],normalized=TRUE))
  
    if(dbg){
      aux  <- KLde(lprobTable[[i]][,"pop"],lprobTable[[i]][,"bnet"])
      if(bsigma){
       names(aux) <- c("KL_pop-bnet","KL_pop-bnet_sigma","KL_pop-bnet_g","H_gexclusion","pexclusion","epsilon")
      }else{
        #names(aux) <- c("KL_pop-bnet","KL_pop-bnet_g","H_gexclusion","pexclusion","epsilon")
        names(aux) <- c("KL_pop-bnet","KL_pop-bnet_g","H_gexclusion","pexclusion","epsilon")
      }       
    }else{
      aux        <- KLd(lprobTable[[i]][,"pop"],lprobTable[[i]][,"bnet"])
      #names(aux) <- "KL_pop-bnet"
      names(aux) <- c("KL_pop-bnet","pexclusion","epsilon")
      
    }
    klDiv     <- rbind(klDiv,c('KL_bnet-pop'=KLd(lprobTable[[i]][,"bnet"],lprobTable[[i]][,"pop"]),aux))
        #exP1      <- rbind(exP1,Px(p1=lprobTable[[i]][,"bnet"],p0=lprobTable[[i]][,"pop"],dbg=dbg))
  } 
  colnames(dH)        <- 'deltaH'
  colnames(dHnorm)    <- 'deltaHnorm'
  colnames(resCrossH) <- c("crossH_pop-bnet","crossH_bnet-pop") 
  #colnames(klDiv)     <- c("KL_bnet-pop","KL_pop-bnet")
  #colnames(exP)       <- c('exP')
  rownames(resCrossH)<-rownames(klDiv)<-names(lprobTable)
  df <- data.frame(dH,dHnorm,resCrossH,klDiv)
  
  # Puede ser que la funcion haya sido llamada con una lista nombrada con genotipos
  # indicando una realizacion hecha con un genotipo dado.
  # En ese caso arreglo el orden alelico en genotipo
  if(length(grep('/',rownames(df)))>0){
   aux <- rownames(df)
   aux <- unlist(lapply(lapply(strsplit(aux,"/"),function(x){sort(as.numeric(x))}),paste0,collapse="/"))
   rownames(df) <- aux
  }
  return(df) 
}
