
"IPA.sampler" <- function(B, P, no.its=1100, burn=100, delta=1,allsamp=F, unknown.pen=NULL){
    Nc <- nrow(B)
    M <- nrow(P)
    
    sampcomp <- apply(P,1,multsample)
    pot <-colSums(B[sampcomp,])## all possible annotations to each assignment
    
    allsampcomp <- matrix(0,no.its, M) # matrix of all the samples
    
    for (it in 1:no.its){
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine){
            p <- pot - B[sampcomp[thism],]
            if(!is.null(unknown.pen)){
              p[length(p)] <- unknown.pen
            } 
            p <- (p + delta)/sum(p+delta)
            po <- p*P[thism,]
            po <- po/sum(po)
            oldval <- sampcomp[thism]
            sampcomp[thism] <- multsample(po)
            if(oldval!=sampcomp[thism]){
                 pot <- pot - as.numeric(B[,oldval]>0) 
                 pot <- pot + as.numeric(B[,sampcomp[thism]]>0)
                 }
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
