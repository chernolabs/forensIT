Px<-function(p1,p0,dbg=FALSE){
  px1 <- minusHp <- NA
  i0 <- which(p1==0)
  if(length(i0)>0){
    p  <- p0[i0]
    px1 <- minusHp <- 0
    for(i in seq_along(p)){
#      accum1 <- accum1 + p[i] * prod(1-p[-i])
      minusHp <- accum1 + p[i]*log10(p[i])
      px1     <- px1 + p[i] 
    }
  }
  if(!dbg){
   return(c(px_1=px1))
  }else{
    return(c(px_1=px1,
             KL_nonxalleles=sum(p1[-i0]*log10(p1[-i0]/p0[-i0])),
             H_xalleles = -minusHp))
  }
}
