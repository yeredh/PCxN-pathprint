# 06/19/2016
#
# Script to check for missing experiment-level estimates
#
# Once the missing estimates are determined, the script will submit
# the appropriate job to get them
#
# ODYSSEY (INTERACTIVE)
#
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem-per-cpu=2000 -t 0-3:30 --pty R


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

# set working directory to output directory
setwd(paste0(pcxn_dir,"output/pathprint_v1.2.3/"))
# get RDS files with the experiment-level estimates
tmp = list.files(pattern="_pathprint_pathcor.RDS")
# Get GSE series with available results
gse_done = gsub("(GSE[0-9]+).*","\\1",tmp)
# determine GSE with missing results
gse_miss = setdiff(gse_ids,gse_done)
cat("Missing:",length(gse_miss),"\n")
miss_ind = which(gse_ids %in% gse_miss)

if(length(miss_ind) > 0){
  # split missing results by consecutive indices
  res = list()
  ic = 1
  res[[ic]] = miss_ind[1]
  if(length(miss_ind) >= 2){
    for(k in 1:(length(miss_ind)-1)){
      if(miss_ind[k+1] - miss_ind[k] == 1){
        res[[ic]] = c(res[[ic]],miss_ind[k+1])
      }else{
        res = c(res,list(miss_ind[k+1]))
        ic = ic + 1
      }
    }
  }
  
  # get upper and lower bounds for job array indices
  job_lims = lapply(res,function(x){
    if(length(x) > 1){
      return(c(head(x,1),tail(x,1)))
    }else{
      return(x)
    }
  })

  
  # submit SLURM scripts
  setwd(paste0(pcxn_dir,"src"))
  for(i in 1:length(job_lims)){
    # command line to submit job
    if(length(job_lims[[i]]) == 1){
      cmd_str = paste0("sbatch --array=",job_lims[[i]]," pcxn_pathprint_estimates01.sh")
    }else{
      cmd_str = paste0("sbatch --array=",job_lims[[i]][1],"-",job_lims[[i]][2]," pcxn_pathprint_estimates01.sh")
    }
    cat(cmd_str,"\n")
    # submit job
    system(cmd_str)
    Sys.sleep(0.5)
  }
}
  
