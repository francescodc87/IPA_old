### this function breaks if you only have two masses


"IPA.sampler.int" <- function(B, P, Int, no.its=1100, burn=100, delta=1,allsamp=F, unknown.pen=NULL, ratio.toll=0.8){
    Nc <- nrow(B)
    M <- nrow(P)
    
    
    sampcomp <- apply(P,1,multsample)
    allsampcomp <- matrix(0,no.its, M) # matrix of all the samples
    for (it in 1:no.its){
        ordine <- sample(M)  # randomising the order used to check all assignments
        for (m in 1:M){
            thism <- ordine[m]
            tmp <- matrix(B[sampcomp[-thism],],ncol = Nc)*(Int[thism]/Int[-thism])
            ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll))) 
            tmp[ind.ones]<-1
            tmp[tmp!=1] <-0
            p<-colSums(tmp)
            
            
            if(!is.null(unknown.pen)){
              p[length(p)] <- unknown.pen
            } 
            p <- (p + delta)/sum(p+delta)
            #p <- p/sum(p)
            po <- p*P[thism,]
            po <- po/sum(po)
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
