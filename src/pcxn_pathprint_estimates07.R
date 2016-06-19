# 05/23/2016
#
# The purpose of this script to combine all the results into a data frame
#
# ODYSSEY (INTERACTIVE SESSION)
# 
# source new-modules.sh
# module load R/3.2.2-fasrc01
# export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE
# srun -p irizarry,serial_requeue --mem=24000 -t 0-2:00 --pty R

rm(list=ls())
gc()
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# ==== Results ====
over_coeff = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/over_coef.RDS"))
pathway_names = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/pathway_names.RDS"))
r_bar = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/r_bar.RDS"))
p_combined = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/res/p_combined.RDS"))


# ==== PCxN results ====
pathCor = data.frame(Pathway.A = pathway_names[,1],Pathway.B = pathway_names[,2],
                     PathCor = r_bar, p.value = p_combined, Overlap.Coefficient=over_coeff)


# remove NAs
pathCor = pathCor[!is.na(pathCor$PathCor),]
# adjust p-values for multiple comparison
pathCor$p.Adjust = p.adjust(pathCor$p.value,"fdr")

saveRDS(pathCor,paste0(pcxn_dir,"/output/pathprint_v1.2.3/res/pathCor_pathprint_v1.2.3_dframe.RDS"))
rm(list=ls())
gc()
        
        