# Script: Plotting the expression of different gene sets (PAM50 and top DEGs)
#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB or METABRIC

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/2_transcriptomic/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
#library(Hmisc)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)

#######################################################################
# load data for boxplot comparisons
#######################################################################

# for SCANB
if (cohort=="SCANB") {
  
  # load annotation data and select subgroup data
  anno <- loadRData(file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>%
    dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50)
  
  # load gex data
  gex.data <- scanb_gex_load(gex.path = "data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata", geneanno.path = "data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata", ID.type = "Gene.Name") %>% 
    dplyr::select(any_of(anno$sampleID)) %>% # select subgroup gex 
    select_if(~ !any(is.na(.))) # otherwise error when scaling
  
  # log transform FPKM data
  gex.data <- as.data.frame(log2(gex.data + 1))
  # z-transform
  gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases
  
  #---------------------------------------------------------------------#
  
} else if (cohort=="METABRIC") {
  
  # load annotation data
  load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
  
  # extract relevant variables
  anno <- loadRData("data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
    dplyr::rename(sampleID=METABRIC_ID) # rename to match SCANB variables
  
  # load and select subgroup data
  gex.data <- metabric_gex_load("./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt",ID.type = "Hugo_Symbol") %>% 
    dplyr::select(any_of(anno$sampleID)) %>% 
    mutate_all(function(x) as.numeric(x))
  
  # exclude samples from anno without associated gex data
  anno <- anno %>% 
    filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# GENESET 1: PAM50 DEGs 
#######################################################################

# load data
DE.res <- loadRData(paste(data.path,"DE_results.RData",sep=""))
pam50.genes <- as.data.frame(read_table("data/SCANB/1_clinical/raw/Spiral_SRIQ_PAM50_Gex_Clusters_6_Info_ann.txt")) %>% pull(HGNC)

# DEGs
DEGs <- DE.res %>% 
    filter(Her2.LumA.padj <= 0.05) %>% 
    filter(Her2.LumB.padj <= 0.05) %>% 
    dplyr::select(Her2.LumA.de, Her2.LumA.padj, Her2.LumB.de, Her2.LumB.padj) %>% rownames_to_column("Gene")

pam50.DEGs <- DEGs %>% 
    filter(Gene %in% pam50.genes) %>% 
    pull(Gene)

# plot expression
plot.list <- list()

for (i in 1:length(pam50.DEGs)) {
    geneID <- pam50.DEGs[i]
    gex <- get_gex(geneID,gex.data,anno)
    
    # base statistics
    #get_stats(gex,"PAM50",geneID)
    
    # test
    res <- pair_ttest(gex, 
                      group.var = "PAM50",
                      test.var = geneID, 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    
    # plot
    plot.list <- append(plot.list, list(
      three_boxplot(gex,
                    group.var = "PAM50",
                    test.var = geneID,
                    g1="Her2",g2="LumA",g3="LumB",
                    colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("Her2","LumA","LumB")),
                    ylab = "Expression (log2)", 
                    title = paste(
                      geneID," expression (LumA_pval=",
                      as.numeric(gsub('^.{2}', '', sub('.+p-value (.+)', '\\1', res[6]))),
                      "; LumB_pval= ",
                      as.numeric(gsub('^.{2}', '', sub('.+p-value (.+)', '\\1', res[19]))),
                      "; ERpHER2n)",sep = ""))))
}

# save plots
pdf(file = paste(output.path,cohort,"_HER2n_PAM50DEGs_expression.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()

#######################################################################
# GENESET 1: PAM50 DEGs 
#######################################################################
# quick check FGF19
# head(DE.res)
# DE.res <- loadRData(paste(data.path,"DE_results.RData",sep=""))
# DE.res["FGF19",]

# top DEGs (high difference)
top.DEGs <- DE.res %>% 
    filter(Her2.LumA.padj <= 0.05) %>% 
    filter(Her2.LumB.padj <= 0.05) %>% 
    filter(abs(Her2.LumA.diff) >= 1.5) %>% 
    filter(abs(Her2.LumB.diff) >= 1.5) %>% 
    dplyr::select(Her2.LumA.de, Her2.LumA.padj, Her2.LumB.de, Her2.LumB.padj) %>% rownames_to_column("Gene") %>% pull(Gene)

# plot expression
plot.list <- list()

for (i in 1:length(top.DEGs)) {
  geneID <- top.DEGs[i]
  gex <- get_gex(geneID,gex.data,anno)
  
  # base statistics
  #get_stats(gex,"PAM50",geneID)
    
  # test
  res <- pair_ttest(gex, 
                    group.var = "PAM50",
                    test.var = geneID, 
                    g1 = "Her2", g2 = "LumA", g3 = "LumB")
  
  # plot
  plot.list <- append(plot.list, list(
      three_boxplot(gex,
                    group.var = "PAM50",
                    test.var = geneID,
                    g1="Her2",g2="LumA",g3="LumB",
                    colors=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("Her2","LumA","LumB")),
                    ylab = "Expression (log2)", 
                    title = paste(
                      geneID," expression (LumA_pval=",
                      as.numeric(gsub('^.{2}', '', sub('.+p-value (.+)', '\\1', res[6]))),
                      "; LumB_pval= ",
                      as.numeric(gsub('^.{2}', '', sub('.+p-value (.+)', '\\1', res[19]))),
                      "; ERpHER2n)",sep = ""))))
}

# save plots
pdf(file = paste(output.path,cohort,"_HER2n_topDEGs_expression.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()
