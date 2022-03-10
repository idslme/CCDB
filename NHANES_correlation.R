## Version March 10th 2022
## The script will work only in a linux machine. 
## Donwnload the data dictionary and matrix from zenodo repository.

Link 1 : https://zenodo.org/record/6345265/files/master_cdf.sb.RData?download=1
Link 2 : https://zenodo.org/record/6345265/files/ndf.list.sb.RData?download=1

------------------------ CODE START below----------
setwd("path to the working directory")

studyid <- "nhanes"

## Load the data dictionary and the list of data matrices. 
load(paste0("master_cdf.sb.RData"))
load(paste0("ndf_list.sb.RData"))

ndf.mat <- do.call(cbind, ndf.list)
ndf.list <- NA
gc()

ndf.mat[which(ndf.mat==0)] <- NA
ndf.mat <- log(ndf.mat, base = 2)

library(parallel)
resl <- mclapply(1:nrow(master_cdf), function(k) {
  vec1 <- as.numeric(ndf.mat[,k])
  qindex <- which(is.na(ndf.mat[,k])==F)

  na_len_overlap <- sapply(1:ncol(ndf.mat), function(x){
    length(which(qindex %in%which(is.na(ndf.mat[,x])==F)==T))
  })
  gc()
  corvec <- rep(0,ncol(ndf.mat))
  ndf.mat.sb <- ndf.mat[,which(na_len_overlap > 100)]
  corvec1 <- corvec[which(na_len_overlap > 100)]
  tryCatch(corvec1 <- WGCNA::cor(x = ndf.mat.sb, y = vec1, use = "pairwise.complete.obs"), error=function(e) {})
  corvec1 <- round(corvec1, digits = 2)
  corvec[which(na_len_overlap > 100)] <- corvec1
  save(corvec, file=paste0(master_cdf$compound_id[k],".RData"))
  ndf.mat.sb <- NA
  #ndf.mat <- NA
  na_len_overlap <- NA
  corvec <- NA
  corvec1 <- NA
  gc()
  return(1)
},mc.cores = 5)