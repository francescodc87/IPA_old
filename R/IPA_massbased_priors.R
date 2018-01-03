#IPA code - R vesion - ONLY MASSES - COMPUTING ONLY THE PROBABILITIES ONLY BASED ON MASS
#mass:  the measured mass
#compounds.mass: a vector containing the C theoretical masses 
# ppm: accuracy, this parameter is used to compute the precision parameter (gamma in the paper). It should be a multiple of the
#######averege accuracy of the instrument


#list of optional parameters
###precision: mass noise precision (gamma in the paper), default 1  ### gotta change this with accuracy 
###noisetype: 0 for'absolute' or 1 for 'relative' noise (default 0)
### pr.limit: the minimum probability accepted (default = 1e-06)

#out: list of outputs
### pr: a vectorcontaining C prior probabilities


"IPA.massbased.prior" <- function(mass, compounds.mass, ppm, pr.limit =1e-02, unknown.prob=0.05,
                                  prior.knowledge=rep(1,length(compounds.mass)),compounds.id=NULL,
                                  RTs=NULL, RT.ranges=NA, RT.prior.penality=.5,
                                  v=FALSE, it=100){
  
  
  "solve.RT.range" <- function(RT.range,min.max){
      RT.range <- unlist(strsplit(RT.range, split=","))
      if(min.max=="min"){
        RT.range <- min(RT.range)
      }else{
        RT.range <- max(RT.range)
      }
      return(as.numeric(RT.range))
  }
  
  #loading packages..no need in the future
  library(Matrix)
  compounds.mass <- as.numeric(compounds.mass)
  #checking the parameters
  if(ppm <=0){
    cat('\n not allowed parameters values')
    stop()
  }
  if(length(compounds.mass)!=length(prior.knowledge)){
    cat('\n compounds.mass and prior.knowledge are not compatible')
    stop()
  }
  
  # defing the number of masses and the number of the compounds
  Nc <- length(compounds.mass)
  M <- length(mass)
  ###evaluating precision
  deltaMs <- ppm*mass*(1e-06)
  sigma <- deltaMs/2
  precision <- 1/(sigma^2)
  rm(sigma, deltaMs)
  # evaluation of prior probabilities (likelihood based only on mass)
  # initialize some variables
  pr <-Matrix(0,M,(Nc+1))
  
  
  unknown.probs <- NULL
  
  #considering RTs and RT ranges for the 
  if(length(RT.ranges)==Nc){
    RT.min <- apply(matrix(RT.ranges,ncol=1),1,solve.RT.range,min.max="min")
    RT.max <-  apply(matrix(RT.ranges,ncol=1),1,solve.RT.range,min.max="max")
  }
  
  
  
  for(i in 1:M){
    RT.prior <-rep(1,Nc)
    
    if(length(RT.ranges)==Nc){
    ind.RT.out <- which(!is.na(RT.ranges))
    ind.RT.out <- ind.RT.out[which(RTs[i]<RT.min[ind.RT.out] | RTs[i]>RT.max[ind.RT.out])]
    RT.prior[ind.RT.out] <- RT.prior.penality
    }

    pr[i,1:Nc] <- (exp((-0.5*precision[i])*((compounds.mass-mass[i])^2)))*prior.knowledge*RT.prior
    if(!is.na(pr.limit)){
      ind.zeros <- which(pr[i,]<(pr.limit*sum(pr[i,])))
      if(length(ind.zeros)>0){
        pr[i,ind.zeros] <- 0  
      }
    }
    unknown.probs  <- (unknown.prob/(1-unknown.prob))*sum(pr[i,])
    if(unknown.probs==0){unknown.probs <- 1}
    pr[i,Nc+1] <- unknown.probs
    pr[i,] <- pr[i,]/sum(pr[i,])
  
    if(v){
      if(i %% it==0) {
        # Print on the screen some message
        cat(paste0((i*100)/M, "%", "\n"))
      }
    }  
  }
  if(!is.null(compounds.id)){
    colnames(pr) <- c(compounds.id,"unknown")
  }
  
  
  
  
    return(pr)
}
