# 05/23/2016
# Yered H.
#
# The purpose of this script is to get the following
# - Number of samples per GSE series
# - Overlap coefficient between pathway pairs
# - Names of all pathway pairs
#
# ODYSSEY
# pcxn_pathprint_estimates04.sh 

rm(list=ls())
gc()
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# ==== GSE annotation ====
gse_annot = readRDS(paste0(pcxn_dir,"data/GSE_annotation.RDS"))
# get sample size per GSE series
gse_count = table(gse_annot$GSE)
gse_count = sort(gse_count,decreasing=T)
# keep series with at least 15 samples
gse_ids = names(gse_count[gse_count >= 15])

# ==== Pathway Annotation ====
gs_lst = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/pathprint_gs.RDS"))
n_pairs = choose(length(gs_lst),2)

library(parallel)
nc = detectCores()

# ==== Sample Size ====
# helper function to get the number of samples for each experiment (GSE)
GetN = function(ic){
  gse=gse_ids[ic]
  gse_lst = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/",gse,"_pathprint_pathcor.RDS"))
  setTxtProgressBar(pb,ic)
  return(gse_lst[[1]]$n)
}

cat("\n Sample Size \n")
pb = txtProgressBar(min=0,max=length(gse_ids),style=3,initial=0)
cat("\n")
n_vec = mclapply(1:length(gse_ids),GetN,mc.cores=nc)
n_vec = simplify2array(n_vec, higher = TRUE)
close(pb)

# save results
saveRDS(n_vec,paste0(pcxn_dir,"output/pathprint_v1.2.3/res/n_vec.RDS"))

# load results for a single experiment 
# (since the overlap coefficient and the pathway names are all the same)
gse_lst = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/",gse_ids[1],"_pathprint_pathcor.RDS"))
# ==== Overlap Coefficient ====
# helper function to get the overlap coefficient
over_coef = mclapply(gse_lst,function(x){
  if(length(x) == 7){
    return(x$Overlap.Coeff)
  }else{
    return(NA)
  }},mc.cores=nc)
over_coef = unlist(over_coef)

# save results
saveRDS(over_coef,paste0(pcxn_dir,"output/pathprint_v1.2.3/res/over_coef.RDS"))

# ==== Pathway Names ====
# helper function to get the pathway names for each pair
pathway_names = mclapply(gse_lst,function(x){
  if(length(x) == 7){
    return(c(x$Pathway.A,x$Pathway.B))
  }else{
    return(c(NA,NA))
  }},mc.cores=nc)

pathway_names = t(simplify2array(pathway_names, higher = TRUE))

# save results
saveRDS(pathway_names,paste0(pcxn_dir,"output/pathprint_v1.2.3/res/pathway_names.RDS"))

# clean workspace
rm(list=ls())
gc()
