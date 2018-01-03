#load("/home/mqbpqfd2/Desktop/new_IPA_2017/ecmdb_database/ECMDB_input_for_structured_DB.Rdata")
#load("/home/mqbpqfd2/Desktop/new_IPA_2017/functions/adducts_updated.Rdata") # here I  have to substitute data("adducts)

"formula_isopattern" <- function(isotable, isotopes){
  library(stringr)
  NM <- dim(isotable)
  rnames <- rep(NA,NM[1])
  new.names <- NULL
  for(i in 1:NM[1]){
    tmp <- isotable[i,3:NM[2]]
    tmp.elements <- names(tmp)
    tmp.elements.char <- (str_extract(tmp.elements, "[aA-zZ]+"))
    tmp.elements.num <- (str_extract(tmp.elements, "[0-9]+"))
    tmp.element.char.uni <- unique(tmp.elements.char)
    for(j in 1:length(tmp.element.char.uni)){
      ind.tmp <- which(tmp.elements.char%in%tmp.element.char.uni[j])
      tmp.elements[ind.tmp[1]] <- tmp.elements.char[ind.tmp[1]]
      if(length(ind.tmp)>1){
        ind.tmp <- ind.tmp[2:length(ind.tmp)]
        tmp.elements[ind.tmp] <- paste("[", tmp.elements.num[ind.tmp], "]", tmp.elements.char[ind.tmp], sep="")
      }
    }
    keep <- which(tmp!=0)
    tmp <- paste(tmp.elements[keep], tmp[keep], sep="")
    tmp <- paste(tmp, collapse="")
    new.names <- c(new.names, tmp)
  }
  rownames(isotable) <- new.names
  return(isotable)
}





### I want a function for the structured DB creator



"structured.DB.creator" <- function(input.table, ionisation, adducts.fully.connected=F, adducts.matrix,
                                    iso.prob.threshold= 5){ # change to 5%
  ##I should add some checks in order to be sure that the input table is in the correct format
  library(enviPat)
  library(Matrix)
  data("isotopes")
  
  if(ionisation=="positive"){
    input.table <- input.table[,c(1:6,9:12)]
    colnames(input.table) <- c("ID","Name","Formula","monoisotopic mass", "Main adduct",
                               "adducts and fragments", "possible reactions", "reported connections",
                               "prior knowledge", "RT range") 
  }else if(ionisation=="negative"){
    input.table <- input.table[,c(1:4,7:12)]
    colnames(input.table) <- c("ID","Name","Formula","monoisotopic mass", "Main adduct",
                               "adducts and fragments", "possible reactions", "reported connections",
                               "prior knowledge", "RT range") 
    
  }else{
    cat("\n not valid ionisation")
    stop()
  }
  
  exact.masses.table <- matrix("",0,10)
  colnames(exact.masses.table) <- c("ID", "name", "Formula", "m/z", "adduct type",
                                    "main", "isotope", "abundance", "prior knowledge", "RT range")
  for(i in 1:nrow(input.table)){
    adducts <- input.table[i,"adducts and fragments"]
    adducts <- unlist(strsplit(adducts, split=","))
    
    for (add in 1:length(adducts)){
      ind <- which(adducts.matrix[,1]==adducts[add])
      if(length(ind)==0){
        cat("\n ERROR:", adducts[add], "is not present in the adducts matrix")
        stop()
      }
      flag.isopattern <- 0
      #First, multiply the chemical formula of the molecule by the times it appears in the final adduct; multiform
      formula.tmp <- multiform(input.table[i,"Formula"], adducts.matrix[ind,4])
      
      #Second, add the chemical formula of any adduct to that of the molecule; mergeform
      if(adducts.matrix[ind,7]!="FALSE"){ ### if there is something to aDD
        formula.tmp <- mergeform(formula.tmp, adducts.matrix[ind,7])
        flag.isopattern <- 1
      }
      #Third, subtract the chemical formula of any deduct from that of the molecule; check_ded & subform.
      if(adducts.matrix[ind,8]!="FALSE"){ ##if there is something to subtract
        if(check_ded(formula.tmp,adducts.matrix[ind,8])=="FALSE"){ ##skip this iteration if I can't do this
          formula.tmp <- subform(formula.tmp,adducts.matrix[ind,8])
          if(formula.tmp!="NANA") {flag.isopattern <- 1}
        }
      }
       
        #Finally, calculate the isotopic fine structure using the correct charge argument in isopattern.
        # if the flag is ==1
        if(flag.isopattern==1){
          
        
        isotable <- isopattern(isotopes = isotopes, formula.tmp, threshold = iso.prob.threshold,
                             charge = adducts.matrix[ind,3],verbose=FALSE)[[1]]
      
        isotable <- formula_isopattern(isotable, isotopes)
        for(iso in 1:nrow(isotable)){
          isotope.f <-1
          if(iso==1){isotope.f <- 0}
          main <- 0
          if(input.table[i,"Main adduct"]==adducts[add] & isotope.f==0){main <- 1}
          exact.masses.table <- rbind(exact.masses.table, c(input.table[i,"ID"], input.table[i,"Name"],
                                                          rownames(isotable)[iso], isotable[iso,"m/z"],
                                                          adducts[add], main, isotope.f, isotable[iso,"abundance"],
                                                          input.table[i,"prior knowledge"], input.table[i,"RT range"]))
      }
    }    
    
  
 
    }# end add loop
  
    
    
  } ## end i loop
  ### in the input.table the information regarding the biochemical connections can be provided in two forms
  ### 1) list of "possible reaction", if for example this section contains H2O the main adduct of this compound will
  ### be connected to all compounds that show a formula which is exactly the initial compound +/- H2O
  ### 2) reported connections just contains a list of comma separated IDs of molecules connected to it
  entry.names <- apply(exact.masses.table[,1:3],MARGIN = 1,paste, collapse=",")
  
  ###
  ###a biotransformation matrix
  ###
  Bio.M <- Matrix(0,nrow(exact.masses.table),nrow(exact.masses.table))
  colnames(Bio.M) <- entry.names
  rownames(Bio.M) <- entry.names
  ind.main <- which(exact.masses.table[,"main"]==1 & exact.masses.table[,"isotope"]==0)
  
  #reported connections
  ind.rep.conn <- which(!is.na(input.table[,"reported connections"]))
  for(b in ind.rep.conn){
    rep.connections <- unlist(strsplit(input.table[b,"reported connections"], split=","))
    ind.ones <- which(as.vector(exact.masses.table[ind.main, "ID"]) %in% rep.connections)
    ind.ones <- ind.main[ind.ones]
    Bio.M[b,ind.ones] <- 1
    Bio.M[ind.ones,b] <- 1
    rm(ind.ones)
  }
  rm(ind.rep.conn)
  ##possible reactions
  ind.pos.react <- which(!is.na(input.table[,"possible reactions"]))
  for(i in ind.pos.react){
    possible.reactions <- unique(unlist(strsplit(input.table[i,"possible reactions"], split=",")))
    possible.reactions <- check_chemform(isotopes, possible.reactions)
    possible.reactions <- possible.reactions[possible.reactions[,1]==FALSE,2]
    for(r in possible.reactions){
     tmp <- subform(input.table[i,"Formula"], r)  
     ind.curr <- which(exact.masses.table[,"ID"]==input.table[i,"ID"] & exact.masses.table[,"main"]==1)
     ind.ones <- which(input.table[,"Formula"]==tmp)
     ind.ones <- which(exact.masses.table[,"ID"]%in%input.table[ind.ones,"ID"] & exact.masses.table[,"main"]==1)
     if(length(ind.ones)>0){
       Bio.M[ind.curr, ind.ones] <- 1
       Bio.M[ind.ones, ind.curr] <- 1  
     }
    }
  }
  
  #rm(tmp, possible.reactions, ind.curr, ind.ones)
  
  ###
  ### a adduct matrix (fully connected or not)
  ###
  Iso.M <- Matrix(0,nrow(exact.masses.table),nrow(exact.masses.table))
  colnames(Iso.M) <- entry.names
  rownames(Iso.M) <- entry.names
  
  Add.M <- Matrix(0,nrow(exact.masses.table),nrow(exact.masses.table))
  colnames(Add.M) <- entry.names
  rownames(Add.M) <- entry.names
  for(i in 1:nrow(input.table)){
    if(adducts.fully.connected==TRUE){
      ind <- which(exact.masses.table[,"ID"]==input.table[i,"ID"] & exact.masses.table[,"isotope"]==0)
      if(length(ind)>0){
        ind <- t(combn(ind,2))
        Add.M[ind] <- 1
        Add.M[cbind(ind[,2], ind[,1])] <- 1  
      }
      rm(ind)
    }else{
      ind.central <- which(exact.masses.table[,"ID"]==input.table[i,"ID"] & exact.masses.table[,"main"]==1)
      ind.other <- which(exact.masses.table[,"ID"]==input.table[i,"ID"] & exact.masses.table[,"main"]==0 & exact.masses.table[,"isotope"]==0)
      if(length(ind.other)>0){
        Add.M[cbind(ind.central, ind.other)] <- 1
        Add.M[cbind(ind.other, ind.central)] <- 1  
      }
      
      rm(ind.central, ind.other)
    }
    
    
    
    #### THIS IS THE SECTION I HAVE TO CHANGE
    ##### WHAT TO I WANT?
    ###
    ###isotopes matrix (which contains intensity ratios and it's not symmetric!)
    ###
    ind <- which(exact.masses.table[,"ID"]==input.table[i,"ID"])
    adducts <- exact.masses.table[ind,"adduct type"]
    adducts <- unique(adducts)
    for(j in 1:length(adducts)){
      ## now I should consider all the entries with the same ID and the same
      ## adduct type (they are isotopes)
      ind.iso <- intersect(ind, which(exact.masses.table[,"adduct type"]==adducts[j]))
      ### now I order them according to abundance
      ind.iso <- ind.iso[order(as.numeric(exact.masses.table[ind.iso,"abundance"]),decreasing = T)]
      if(length(ind.iso)>1){
        for(is in 1:(length(ind.iso)-1)){
          ratio <- as.numeric(exact.masses.table[ind.iso[is],"abundance"])/as.numeric(exact.masses.table[ind.iso[(is+1)],"abundance"])
          Iso.M[ind.iso[is],ind.iso[(is+1)]] <- ratio
          Iso.M[ind.iso[(is+1)],ind.iso[is]] <- 1/ratio
        }
      }
    }
    rm(ind,ind.iso)
  
  }
  

  
  out <- list(exact.masses.table=exact.masses.table, Bio.M=Bio.M, Add.M=Add.M, Iso.M=Iso.M)
  return(out)
}