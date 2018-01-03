rm(list = ls())

## posteriors for synthetic experiment

####################Only  Real compounds
####################Positive
###################ADD+ISO


#### loading functions
source("~/R/checking_corr.R")
source("~/R/checking_relID.R")
source("~/R/checking_RT.R")
source("~/R/compute_post.R")
source("~/R/multsample.R")
source("~/R/IPA_sampler_ADD_ISO_BIO_Int_noPot.R")
source("~/R/IPA_sampler_ADD_ISO_BIO.R")
source("~/R/IPA_sampler_ADD_ISO_Int_noPot.R")
source("~/R/IPA_sampler_ADD_ISO.R")
source("~/R/IPA_sampler_int.R")
source("~/R/IPA_sampler.R")
source("~/R/IPA_sampler_RTrelIDcorr_int.R")
source("~/R/IPA_sampler_RTrelIDcorr.R")

### loading priors

###dataset
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
ind.kept<- which(Masses>=90 & Masses<=476.18 & apply(Intensities,1, max, na.rm=T) > 8.59e+05 & RTs>=38.5 & RTs<=440.31)
Pos.All.data <- Pos.All.data[ind.kept,]
rm(Masses, PeakML.Data, PeakML.Data.Table, RTs, Intensities, ids, rel.id, tmp)
###database
load("~/data/DBs_all_compounds_only1entry.Rdata") ##database
###Priors
load("~/data/POSPriors_oneEntry.Rdata") # Priors
Prior.Pos <- Prior.Pos[ind.kept, ]  ### I have to remove the masses from here as well of course
#save(Prior.Pos,file="~/data/POSPriors_oneEntry_filtered.Rdata")


#### adding unknowns in the connectivity matrices
Adducts.matrix <- DB2.structured.POS$Add.M
Adducts.matrix <- cbind(Adducts.matrix, rep(0,nrow(Adducts.matrix)))
Adducts.matrix <- rbind(Adducts.matrix, rep(0,ncol(Adducts.matrix)))

Isotopes.matrix <- DB2.structured.POS$Iso.M
Isotopes.matrix <- cbind(Isotopes.matrix, rep(0,nrow(Isotopes.matrix)))
Isotopes.matrix <- rbind(Isotopes.matrix, rep(0,ncol(Isotopes.matrix)))

Bio.matrix <- DB2.structured.POS$Bio.M
Bio.matrix <- cbind(Bio.matrix, rep(0,nrow(Bio.matrix)))
Bio.matrix <- rbind(Bio.matrix, rep(0,ncol(Bio.matrix)))

### removing masses and formulas I don't need to consider
ind.masses.kept <- which(Prior.Pos[,ncol(Prior.Pos)]!=1)  
ind.formulas.kept <- which(colSums(as.matrix(Prior.Pos))>0)

Prior.red <- Prior.Pos[ind.masses.kept, ind.formulas.kept]

Bio.matrix <- Bio.matrix[ind.formulas.kept, ind.formulas.kept]

Isotopes.matrix <- Isotopes.matrix[ind.formulas.kept, ind.formulas.kept]

Adducts.matrix <- Adducts.matrix[ind.formulas.kept, ind.formulas.kept] 

### making binary Isotopes matrix
Isotopes.matrix.bin <- Isotopes.matrix
Isotopes.matrix.bin[Isotopes.matrix.bin>0] <- 1


RT <- as.numeric(Pos.All.data[,"RTs"])[ind.masses.kept]
### FOR THE INTENSITIES I AM TAKING THE MAXIMUM
Int <- as.numeric(apply(Pos.All.data[,3:40],1, max, na.rm=T)[ind.masses.kept])





Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen <- IPA.sampler.Add.Iso.Bio.int_noPot(P=Prior.red,Add = Adducts.matrix,
                                                                            Bio = Bio.matrix,
                                                                            Iso = Isotopes.matrix, RT = RT, Int = Int,
                                                                            RT.win = 5, no.its =2000, burn = 1000,
                                                                            allsamp = T, delta.add =0.4,
                                                                            delta.iso = 0.2, delta.bio =1, ratio.toll = 0.5)


tmp.Post <- Prior.Pos
tmp.Post[ind.masses.kept, ind.formulas.kept] <-Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$Post
Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$Post <- tmp.Post
rm(tmp.Post)

tmp.allsampcomp <- matrix(0,2000, nrow(Prior.Pos))
tmp.allsampcomp[,ind.masses.kept] <- Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$allsampcomp
Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen$allsampcomp <- tmp.allsampcomp
rm(tmp.allsampcomp)

#save(Post.real.Pos.oneEntry.ADD.ISO.BIO.Int.noPen, 
#     file="~/data/Post_real_oneEntry_ADD_ISO_BIO_Int_noPen2.Rdata")



