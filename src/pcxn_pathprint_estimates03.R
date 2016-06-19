# 05/23/2016
#
# The purpose of this script is to retrieve the experiment-level
# p-values  and store them into a matrix.
#
# I decided to split the experiments into 10 groups to speed up the process
# the job array index is used to select the group
# 
# RUN ON ODYSSEY (pcxn_pathprint_estimates03.sh) [1-10]

rm(list=ls())
gc()
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# ==== Pathway Annotation ====
gs_lst = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/pathprint_gs.RDS"))
n_pairs = choose(length(gs_lst),2)

# ==== GSE annotation ====
gse_annot = readRDS(paste0(pcxn_dir,"data/GSE_annotation.RDS"))
# get sample size per GSE series
gse_count = table(gse_annot$GSE)
gse_count = sort(gse_count,decreasing=T)
# keep series with at least 15 samples
gse_ids = names(gse_count[gse_count >= 15])


# split experiments into ten groups
n_part = ceiling(length(gse_ids)/10)
gse_split = split(gse_ids, ceiling(seq_along(gse_ids)/n_part))

# read command line arguments
args <- as.numeric(commandArgs(trailingOnly = T))


# matrix to store correlation coefficients
p_mat = matrix(NA,ncol=length(gse_split[[args]]),nrow=n_pairs,dimnames=list(NULL,gse_split[[args]]))

library(parallel)
nc = detectCores()

pb = txtProgressBar(min=0,max=length(gse_split[[args]]),style=3,initial=0)
for(k in 1:length(gse_split[[args]])){
  gse = gse_split[[args]][k]
  gse_lst = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/",gse,"_pathprint_pathcor.RDS"))
  # get experiment-level correlation coefficients
  tmp = mclapply(gse_lst,function(x){
    if(length(x) == 7){
      return(x$p.value)
    }else{
      return(NA)
    }
  },mc.cores=nc)
  tmp = simplify2array(tmp, higher = TRUE)
  # store experiment-level estimates in matrix
  p_mat[,k] = tmp
  rm(tmp)
  gc()
  setTxtProgressBar(pb,k)
}
cat("\n")
close(pb)
cat("\n")

saveRDS(p_mat,paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_mat",args,".RDS"))
rm(list=ls())
gc()