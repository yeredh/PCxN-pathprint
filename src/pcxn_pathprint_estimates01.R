# 06/18/2016
#
# Get experiment-level estimates for pathprint 
#
# The job array index is used to select the GSE series with the script
#
# ODYSSEY
#
# (pcxn_pathprint_estimates01.sh) [1-863]



rm(list=ls())
gc()
library(parallel)
library(matrixStats)
library(corpcor)
options(stringsAsFactors = F)
pcxn_dir = "/net/hsphfs1/srv/export/hsphfs1/share_root/hide_lab/PCxN/"

# ==== PCxN Functions ====
OverlapCoefficient <- function(x,y){
  # function to calculate the overlap coefficient between x and y
  # which is defined as the size of the intersection divided by the
  # size of the smaller set
  #
  # Args
  #   x: a vector
  #   y: a vector
  #
  # Returns
  #   the overlap coefficient, a number between 0 and 1
  
  length(intersect(x,y))/min(length(unique(x)),length(unique(y)))
}

GetSummary = function(dat,gs,sum_fun){
  # function to calculate the summary statistic for the pathway
  #
  # Args.
  #   dat: genes by samples matrix
  #   gs: vector with the names of the genes in the gene set
  #   sum_fun: function to calculate the summary
  #
  # Returns
  #   a 1 by samples vector with the summary statistic for the pathway
  
  if(length(gs) > 1){
    # calculate summary for pathways with more than 1 element
    return(sum_fun(dat[rownames(dat) %in% gs,]))
  }else{
    # return actual value for pathways with a single element
    return(dat[rownames(dat) %in% gs,])
  }
}

ShrinkCor = function(x,y,method="pearson"){
  # wrapper to estimate the correlation coefficient between x and y using the 
  # shrinkage estimator from Schaffer & Strimmer (2005) [corpcor package]
  # and the corresponding t-statistic and p-value
  #
  # Args
  #   x: a vector with n observations
  #   y: a vector with n observations
  #   method: character to pick either the Pearson or Spearman correlation coefficient
  #
  # Returns
  #   a named vector with the correlation estimate, the sample size n, the t-statistic
  #   and its corresponding p-value
  
  # function to get the t-statistic
  GetStatistic <- function(r,n){r*sqrt((n-2)/(1-r^2))}
  # get sample size
  if(length(x) == length(y)){
    n <- length(x)
  }else{
    cat("\n x and y have different lengths! >=( \n")
    return(NA)
  }
  # determine method
  selected_method <- match(method,c("pearson","spearman"))
  # Pearson correlation
  if(selected_method == 1){
    estimate <- cor.shrink(cbind(x,y),verbose=F)
    statistic <- GetStatistic(estimate[2,1],n)
    p.value <- 2*pt(-abs(statistic),n-2)
  }else if(selected_method == 2){
    estimate <- cor.shrink(cbind(rank(x),rank(y)),verbose=F)
    statistic <- GetStatistic(estimate[2,1],n)
    p.value <- 2*pt(-abs(statistic),n-2)
  }else{
    cat("invalid method! >=( \n")
    return(NA)
  }
  # prepare results
  res <- c(estimate[2,1],n,statistic,p.value)
  names(res) <- c("estimate","n","statistic","p.value")
  return(res)
}


# ==== GSE annotation ====
gse_annot = readRDS(paste0(pcxn_dir,"data/GSE_annotation.RDS"))
# get sample size per GSE series
gse_count = table(gse_annot$GSE)
gse_count = sort(gse_count,decreasing=T)
# keep series with at least 15 samples
gse_ids = names(gse_count[gse_count >= 15])

# ==== Pathway Annotation ====
gs_lst = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/pathprint_gs.RDS"))


# ==== Expression Background ====
exprs_rnk = readRDS(paste0(pcxn_dir,"output/pathprint_v1.2.3/GPL570_pathprint_mat.RDS"))

# read command line arguments
args <- as.numeric(commandArgs(trailingOnly = T))
# argument to pick GSE series
cat("GSE Series:", gse_ids[args],"\n")

# select GSE series
gse = gse_ids[args]
gsm_targets = gse_annot$GSM[gse_annot$GSE == gse]
gsm_ind = colnames(exprs_rnk) %in% gsm_targets
# subset expression ranks
exprs_rnk = exprs_rnk[,gsm_ind]
# R garbage collection
gc()

# helper function to get the experiment-level estimates for a 
# gene-set pair
ProcessElement = function(ic,METHOD,sum_fun){
  i = ceiling((sqrt(8*(ic+1)-7)+1)/2)
  j = ic-choose(floor(1/2+sqrt(2*ic)),2)
  
  # pathway gene sets
  gsA=gs_lst[[i]]
  gsB=gs_lst[[j]]
  
  # get genes unique to each gene set
  # 1. genes unique to pathway A
  gsAu <- setdiff(gsA,gsB)
  # 2. genes unique to pathway B
  gsBu <- setdiff(gsB,gsA)
  
  if(length(gsAu) == 0 | length(gsBu) == 0){
    
    cat("\n-------------------------------------\n")
    cat("\n",ic," no unique genes \n")
    cat(names(gs_lst)[i],"\n")
    cat(length(gsAu),"\n")
    cat(names(gs_lst)[j],"\n")
    cat(length(gsBu),"\n")
    cat("\n-------------------------------------\n")
    
    fu = data.frame(Pathway.A=names(gs_lst)[i],Pathway.B=names(gs_lst)[j])
    tmp = rep(NA,4)
    names(tmp) <- c("estimate","n","statistic","p.value")
    fu = c(fu,tmp)
    # overlap coefficient
    fu$Overlap.Coeff= OverlapCoefficient(gs_lst[[i]],gs_lst[[j]])
    return(fu)
  }
  
  # get pathway summaries for unique genes in each pathway
  summaryAu = GetSummary(dat=exprs_rnk,gs=gsAu,sum_fun=sum_fun)
  summaryBu = GetSummary(dat=exprs_rnk,gs=gsBu,sum_fun=sum_fun)
  
  # get correlation between the summaries for the unique genes
  tmp = data.frame(Pathway.A=names(gs_lst)[i],Pathway.B=names(gs_lst)[j])
  tmp = c(tmp,ShrinkCor(summaryAu,summaryBu,method=METHOD))
  # overlap coefficient
  tmp$Overlap.Coeff= OverlapCoefficient(gs_lst[[i]],gs_lst[[j]])
  
  setTxtProgressBar(pb,ic)
  return(tmp)
}


(nc = detectCores())

# loop through pathways  
number_of_pathways = choose(length(gs_lst),2)
input = 1:number_of_pathways
# ==== PCxN Estimates: Mean/Spearman ====
ProcessElementMnSp = function(ic){ProcessElement(ic,METHOD="spearman",sum_fun=colMeans)}


pb = txtProgressBar(min=0,max=number_of_pathways,style=3,initial=0)
cat("\n")
res = mclapply(input,ProcessElementMnSp,mc.cores=nc)
close(pb)

saveRDS(res,paste0(pcxn_dir,"output/pathprint_v1.2.3/",gse,"_pathprint_pathcor.RDS"))

rm(list=ls())
gc()
