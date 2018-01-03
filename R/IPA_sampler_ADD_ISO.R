
"IPA.sampler.Add.Iso" <- function(P, Add, Iso, RT=NULL, rel.id=NULL,
                                        Corr.matrix=NULL,  RT.win=3,
                                        corr.thr=0.80,no.its=1100, 
                                        burn=100, delta.add=1, delta.iso=1,
                                        allsamp=F,
                                        unknown.pen=NULL){
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
    pot.add <-apply(Add[sampcomp,],2,sum)## all possible annotations to each assignment
    pot.iso <-apply(Iso[sampcomp,],2,sum)
    
    allsampcomp <- matrix(0,no.its, M) # matrix of all the samples
    
    for (it in 1:no.its){
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine){
            
            #counting adducts
            p.add <- pot.add - colSums(matrix(Add[sampcomp[ind.rem[[thism]]],], ncol=Nc))
          
            ###counting isotopes
            p.iso <- pot.iso - colSums(matrix(Iso[sampcomp[ind.rem[[thism]]],], ncol=Nc))
            
            ### adding penalities
            if(!is.null(unknown.pen)){
              p.add[length(p.add)] <- unknown.pen
              p.iso[length(p.iso)] <- unknown.pen
            } 
            
            ## normalising with deltas
            p.add <- (p.add + delta.add)/sum(p.add+delta.add)
            p.iso <- (p.iso + delta.iso)/sum(p.iso+delta.iso)
            
            ## merging scores
            po <- p.add*p.iso*P[thism,]
            po <- po/sum(po)
            oldval <- sampcomp[thism]
            sampcomp[thism] <- multsample(po)
            if(oldval!=sampcomp[thism]){
                 pot.add <- pot.add - Add[,oldval] 
                 pot.add <- pot.add + Add[,sampcomp[thism]]
                 pot.iso <- pot.iso - Iso[,oldval] 
                 pot.iso <- pot.iso + Iso[,sampcomp[thism]]
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
