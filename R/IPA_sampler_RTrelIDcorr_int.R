
"IPA.sampler.RT.relID.corr.int" <- function(B, P, Int, RT=NULL, rel.id=NULL,
                                        Corr.matrix=NULL,  RT.win=3,
                                        corr.thr=0.80,no.its=1100, 
                                        burn=100, delta=1,allsamp=F,
                                        unknown.pen=NULL, ratio.toll=0.8){
    if(is.null(RT) & is.null(rel.id) & is.null(Corr.matrix)){
      cat("\n Missing RT, Corr.matrix or rel.id")
      stop()
    }
    
    if(!is.null(RT)){
      counter <- 0
      ind.rem <- lapply(RT, function(x){
        counter <<- counter + 1
        checking.RT(RT,counter, RT.win)
      })
    }else if(!is.null(rel.id)){
      counter <- 0
      ind.rem <- lapply(rel.id, function(x){
        counter <<- counter + 1
        checking.rel.id(rel.id,counter)
      })
    }else{
      counter <- 0
      ind.rem <- lapply(RT, function(x){
        counter <<- counter + 1
        checking.corr(Corr.matrix[,counter],counter, corr.thr)
      })
      
    }
     
   
    Nc <- nrow(B)
    M <- nrow(P)

    sampcomp <- apply(P,1,multsample)
    
    allsampcomp <- matrix(0,no.its, M) # matrix of all the samples
    
    for (it in 1:no.its){
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine){
          tmp <- matrix(B[sampcomp[-ind.rem[[thism]]],],ncol = Nc)*(Int[thism]/Int[-ind.rem[[thism]]])
          ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll))) 
          tmp[ind.ones]<-1
          tmp[tmp!=1] <-0
          p<-colSums(tmp)  
          if(!is.null(unknown.pen)){
            p[length(p)] <- unknown.pen
          } 
          p <- (p + delta)/sum(p+delta)
          po <- p*P[thism,]
          po <- po/sum(po)
          oldval <- sampcomp[thism]
          sampcomp[thism] <- multsample(po)
          
        }
            allsampcomp[it,] <- sampcomp
    }
    
    
    post <- t(apply(allsampcomp,2,compute.post, burn=burn, no.its=no.its, Nc=Nc))
    if(allsamp){
      out <- list(Post=post, allsampcomp=allsampcomp)
      return(out)
    }else{
      return(post)  
    }
}
