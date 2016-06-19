# 05/23/2016
#
# Combine the experiment-level p-values
# using Liptak's method with experiment sample
# sizes as weights
#
# ODYSSEY
# pcxn_pathprint_estimates06.sh

rm(list=ls())
gc()
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# adjust p-values and correlation estimates,
# otherwise the functions to combine the p-values
# cannot handle values near 0 and values near 1
AdjustPmat = function(p_mat,eps=1E-16){
  res = t(apply(p_mat,1,function(pval){
    pval[pval <= eps] = eps
    pval[pval >= 1-eps] = 1 - eps
    return(pval)
  }))
  return(res)
}

library(metap)
library(pbapply)
# ==== Experiment-level p-values ====
p_mat = cbind(readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat1.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat2.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat3.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat4.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat5.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat6.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat7.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat8.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat9.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat10.RDS")))

# adjust p-values
p_mat = AdjustPmat(p_mat)


# ==== Liptak's Method (SS) ====
# sample size per experiment
n_vec = unlist(readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/n_vec.RDS")))

# Liptak Method (sample size as weights)
p_combined = pbapply(p_mat,1,function(x){
  if(any(is.na(x))){
    return(NA)
  }else{
    return(sumz(p=x,weights=n_vec)$p)
  }})

# save results
saveRDS(p_combined,paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_combined.RDS"))

rm(list=ls())
gc()