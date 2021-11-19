# Standard DE Analysis
library(dplyr)
library(edgeR)
# Set work Directory
setwd('C:/Users/harri/Desktop/Bulk_samples')
# Read expression matrix
exp_hm <- read.csv('hm_design.csv')
gene <- exp_hm[2]
exp_hm <- as.matrix(exp_hm[,-c(1,2)])
exp_hs <- read.csv('hs_design.csv')
exp_hs <- as.matrix(exp_hs[,-c(1,2)])
exp_ms <- read.csv('ms_design.csv')
exp_ms <- as.matrix(exp_ms[,-c(1,2)])
# Set groups
hm_group <- factor(c(rep(x='Healthy',23), 
                  rep(x='Moderate',10)))
hs_group <- factor(c(rep(x='Healthy',23), 
                     rep(x='Severe',12,)))
ms_group <- factor(c(rep(x='Moderate',10), 
                     rep(x='Severe',12)))
# DE Analysis procedure
DEAnalysis <- function(counts, group){
  # Create DGElist object
  dgelist <- DGEList(counts = counts, group = group, genes = gene)
  # Preliminary gene filtering
  keep.exprs<- edgeR::filterByExpr(counts, group = group)
  dgelist.filtered <- dgelist[keep.exprs, ,keep.lib.sizes = FALSE]
  # TMM nomalization
  dgelist.filtered.normalized <- calcNormFactors(dgelist.filtered, 
                                                 method = 'TMM')
  dgelist_analysis <- dgelist.filtered.normalized
  # Estimate dispersons
  dgelist_analysis <- estimateDisp(dgelist_analysis)
  # Calculate DE Analysis results
  dgelist_analysis <- edgeR::exactTest(dgelist_analysis)
  dgelist.dtest <- decideTests(dgelist_analysis, adjust.method="BH",
                               p.value=0.05, lfc = 1)
  dgelist.tested <- dgelist_analysis[as.logical(dgelist.dtest),]
  dgelist_analysis$tested <- dgelist_analysis[as.logical(dgelist.dtest),]
  return(dgelist_analysis)
  
}
dgelist_hm <- DEAnalysis(exp_hm, hm_group)
dgelist_hs <- DEAnalysis(exp_hs, hs_group)
dgelist_ms <- DEAnalysis(exp_ms, ms_group)

# View and store DE genes
StoreDE <- function(dgelist, filename){
  write.csv(dgelist, filename)
}
StoreDE(dgelist_hm$tested, 'dgelist_hm.csv')
StoreDE(dgelist_hs$tested, 'dgelist_hs.csv')
StoreDE(dgelist_ms$tested, 'dgelist_ms.csv')

degenes_aggregate = append(dgelist_hm$tested$genes$Gene,
                           dgelist_hs$tested$genes$Gene)
degenes_aggregate = unique(append(degenes_aggregate,
                                  dgelist_ms$tested$genes$Gene))
write(degenes_aggregate, 'degenes.txt')