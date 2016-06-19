# 06/19/2016
#
# The purpose of this script is to estimate the weighted average of the 
# experiment-level correlation coefficients that have been previously
# stored as a matrix
# 
# ODYSSEY
# pcxn_pathprint_estimates05.sh

rm(list=ls())
gc()
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"


library(pbapply)
# ==== Experiment-level correlation estimates ====
r_mat = cbind(readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat1.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat2.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat3.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat4.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat5.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat6.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat7.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat8.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat9.RDS")),
              readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_mat10.RDS")))

dim(r_mat)

# sample size per experiment
n_vec = unlist(readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/n_vec.RDS")))


# weighted average for the correlation estimates
n_mult = n_vec/sum(n_vec)
rm(n_vec)
r_bar = r_mat%*%n_mult

# save results
saveRDS(r_bar,paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_bar.RDS"))


rm(list=ls())
gc()