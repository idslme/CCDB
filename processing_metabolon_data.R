## Version March 10th 2022
## The script will work only in a linux machine with at least 12 CPUs 
## Donwnload the input format from the zenodo repository

Link 1 : https://zenodo.org/record/6345265/files/ST002020_input.xlsx?download=1

setwd("path to the working directory")

datadir <- getwd()
studyid <- "ST002020"

## Always make sure that metabolon's dataset are properly transformred
cdf <- as.data.frame(readxl::read_xlsx(paste0(datadir,studyid,"_input.xlsx"),sheet = "data_dictionary"))
ndf <- as.data.frame(readxl::read_xlsx(paste0(datadir,studyid,"_input.xlsx"),sheet = "data_matrix"))

master_cdf <- cdf[,c("compound_id","BIOCHEMICAL","KEGG","PUBCHEM")]
master_ndf <- ndf[,-1]

master_ndf <- master_ndf[order(master_cdf$BIOCHEMICAL),]
master_cdf <- master_cdf[order(master_cdf$BIOCHEMICAL),]

save(master_ndf, file=paste0(studyid,"_master_ndf.RData"))
save(master_cdf, file=paste0(studyid,"_master_cdf.RData"))

ndf.list <- lapply(1:nrow(master_ndf), function(x){
  as.numeric(master_ndf[x,])
})

save(ndf.list, file=paste0(studyid,"_ndf_list.RData"))

ndf.mat <- do.call(cbind, ndf.list)
na.mat <- is.na(ndf.mat)

na_len_vec <- sapply(1:ncol(ndf.mat), function(x) {
  length(which(!is.na(ndf.mat[,x])))
})
save(na_len_vec, file= paste0(studyid,"_na_len.RData"))

library(parallel)
resl <- mclapply(1:ncol(ndf.mat), function(k) {
  vec1 <- ndf.mat[,k]
  na_len_overlap <- sapply(1:ncol(ndf.mat), function(x){
    length(which(which(!na.mat[,x])%in%which(!na.mat[,k])))
  })
  corvec <- rep(0,ncol(ndf.mat))
  ndf.mat.sb <- ndf.mat[,which(na_len_overlap > nrow(na.mat)/10)]
  corvec1 <- corvec[which(na_len_overlap > nrow(na.mat)/10)]
  tryCatch(corvec1 <- WGCNA::cor(x = ndf.mat.sb, y = vec1, use = "pairwise.complete.obs"), error=function(e) {})
  corvec1 <- round(corvec1, digits = 2)
  corvec[which(na_len_overlap > nrow(na.mat)/10)] <- corvec1
  save(corvec, file=paste0(master_cdf$compound_id[k],".RData"))
  ndf.mat.sb <- NA
  na_len_overlap <- NA
  corvec <- NA
  corvec1 <- NA
  return(1)
},mc.cores = 24)

