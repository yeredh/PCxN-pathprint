# 06/18/2016
#
# Subset and a copy the matrix with the gene expression background to include
# only the genes in pathprint in order to save memory while estimating
# the gene set correlations
#
# ODYSSEY (INTERACTIVE SESSION)
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem-per-cpu=15000 -t 0-2:30 --pty R

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

# ===== Pathprint Pathways =====
pathprint_gs = readRDS(paste0(pcxn_dir,"data/pathprint_v1_2_3.RDS"))


# ==== Expression Background ====
load(paste0(pcxn_dir,"data/GPL570.R.mat.RData"))

# genes in microarray background
entrez_ids = rownames(GPL570.R.mat)
# keep only genes in microarray background
pathprint_r_gs = lapply(pathprint_gs,function(x){x[x %in% entrez_ids]})
# filter empty gene sets (if applicable)
pathprint_rlen = sapply(pathprint_r_gs,length)
pathprint_r_gs = pathprint_r_gs[pathprint_rlen > 0]

# save filtered gene set annotation 
saveRDS(pathprint_r_gs,paste0(pcxn_dir,"output/pathprint_v1.2.3/pathprint_gs.RDS"))


# subset and save expression rank matrix to save memory 
# (i.e. only keep ranks for genes in pathprint gene sets)
saveRDS(GPL570.R.mat[rownames(GPL570.R.mat) %in% unlist(unique(pathprint_r_gs)),],
        paste0(pcxn_dir,"output/pathprint_v1.2.3/GPL570_pathprint_mat.RDS"))

rm(list=ls())
gc()

