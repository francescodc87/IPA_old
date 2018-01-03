rm(list = ls())
library(mzmatch.R)
mzmatch.init(version.1 = F)




PeakML.Data <- PeakML.Read("~/allpeaks_filtered.peakml")

PeakML.Data.Table <- PeakML.Methods.getCompleteTable(PeakML.Data)
dim(PeakML.Data.Table$Intensities)

Masses <- apply(PeakML.Data.Table$Masses, 2, mean, na.rm=TRUE)
RTs <- apply(PeakML.Data.Table$Retentiontimes,2, mean, na.rm=TRUE)
Intensities <- t(PeakML.Data.Table$Intensities)
colnames(Intensities) <- PeakML.Data$sampleNames
ids <- PeakML.Data$GroupAnnotations$id
rel.id <- PeakML.Data$GroupAnnotations$relation.id
Pos.All.data <- cbind(Masses, RTs, Intensities, ids,rel.id)

tmp <- PeakML.Data$phenoData
tmp <- paste(tmp, c(1:7, 1:7,1:7,1:7,1:7,1:3), sep = "")
colnames(Pos.All.data) <- c("m/z", "RTs", tmp, "ids", "rel.ids")




####Priors only one entry per formula
load("~/data/DBs_all_compounds_only1entry.Rdata")
source('~/R/IPA_massbased_priors.R')

compounds.mass = DB2.structured.POS$exact.masses.table[,4]
### for the prior I only need the Masses
Prior.Pos <- IPA.massbased.prior(Masses,
                                 compounds.mass = compounds.mass,ppm = 3,
                                 unknown.prob=0.05, compounds.id=1:length(compounds.mass),v=T)
rownames(Prior.Pos) <- ids

Prior.List <- list()
for(i in 1:nrow(Prior.Pos)){
  ind <- which(Prior.Pos[i,]>0)
  if(length(ind)>0){
    Prior.List[[i]] <- Prior.Pos[i,ind]
    names(Prior.List[[i]]) <- colnames(Prior.Pos)[ind]
  }else{
    Prior.List[[i]] <- NA
  }
}

names(Prior.List) <- rownames(Prior.Pos)

#save(Prior.Pos,Prior.List, file="~/data/POSPriors_oneEntry.Rdata")


