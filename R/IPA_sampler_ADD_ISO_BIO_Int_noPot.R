
"IPA.sampler.Add.Iso.Bio.int_noPot" <- function(P, Add, Iso, Bio, Int,
                                          RT=NULL, rel.id=NULL,
                                          Corr.matrix=NULL,  RT.win=3,
                                          corr.thr=.80,no.its=1100, 
                                          burn=100, delta.add=.5, delta.iso=.5,
                                          delta.bio=1, allsamp=F,
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
     
   
    Nc <- nrow(Add)
    M <- nrow(P)
    
    sampcomp <- apply(P,1,multsample)
    pot.bio <-apply(Bio[sampcomp,],2,sum)
    
    allsampcomp <- matrix(0,no.its, M) # matrix of all the samples
    
    for (it in 1:no.its){
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine){
            
            #counting adducts
            p.add <- colSums(matrix(Add[sampcomp[-ind.rem[[thism]]],], ncol=Nc))
            ###counting isotopes
            tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]],],ncol = Nc)*(Int[thism]/Int[-ind.rem[[thism]]])
            ind.ones <- which((tmp>=ratio.toll) & (tmp <=(1/ratio.toll))) 
            tmp[ind.ones]<-1
            tmp[tmp!=1] <-0
            p.iso<-colSums(tmp)
            
            ##counting biotransformations
            p.bio <- pot.bio - colSums(matrix(Bio[sampcomp[thism],], ncol=Nc))
            
            ### adding penalities
            if(!is.null(unknown.pen)){
              p.add[length(p.add)] <- unknown.pen
              p.iso[length(p.iso)] <- unknown.pen
             
            } 
            
            ## normalising with deltas
            p.add <- (p.add + delta.add)/sum(p.add+delta.add)
            p.iso <- (p.iso + delta.iso)/sum(p.iso+delta.iso)
            p.bio <- (p.bio + delta.bio)/sum(p.bio+delta.bio)
            
            ## merging scores
            po <- p.add*p.iso*p.bio*P[thism,]
            po <- po/sum(po)
            oldval <- sampcomp[thism]
            sampcomp[thism] <- multsample(po)
            if(oldval!=sampcomp[thism]){
                 pot.bio <- pot.bio - Bio[,oldval] 
                 pot.bio <- pot.bio + Bio[,sampcomp[thism]]
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
