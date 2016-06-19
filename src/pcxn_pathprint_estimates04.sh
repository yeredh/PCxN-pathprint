#!/bin/sh
#SBATCH -J pcxn_pathprint_res # A single job name for the array
#SBATCH -n 16 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue,irizarry # Partition
#SBATCH --mem 12000 # Memory request
#SBATCH -t 0-00:07 # (D-HH:MM)
#SBATCH -o /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_pathprint_estimates04.out # Standard output
#SBATCH -e /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_pathprint_estimates04.err # Standard error
source new-modules.sh
module load R/3.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE

R CMD BATCH /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/src/pcxn_pathprint_estimates04.R /net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/log/pcxn_pathprint_estimates04.Rout
