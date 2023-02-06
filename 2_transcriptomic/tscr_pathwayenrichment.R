# Script: Pathway enrichment analyses for different gene sets
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
  anno <- as.data.frame(
    read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx")) %>%
    filter(Follow.up.cohort==TRUE) %>% 
    filter(NCN.PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
    filter(ER=="Positive" & HER2=="Negative") %>% 
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
  anno <- anno %>% 
    filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
    filter(grepl('ERpHER2n', ClinGroup)) %>% 
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

#-----------------------plot expression-------------------------------#

# plot expression
plot.list <- list()

if (cohort == "SCANB") {
    for (i in 1:length(pam50.DEGs)) {
        geneID <- mart.res %>% 
            filter(hgnc_symbol == pam50.DEGs[i]) %>% 
            pull(ensembl_gene_id)
        gex <- get_gex(geneID,gex.data,anno)
        # base statistics
        get_stats(gex,"PAM50",geneID)
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
                          g3.pos = 5.5, g3.sign = res[2,3],
                          g2.pos = 6.5, g2.sign = res[1,3],
                          ylim = c(-7,7),
                          ylab = "Expression (log2)", 
                          title = paste(pam50.DEGs[i]," (",geneID,") expression in PAM50 subtypes (ERpHER2n;",cohort,")",sep = ""))))
    }
} else if (cohort == "METABRIC") {
    geneID <- pam50.DEGs[i]
    for (i in 1:length(pam50.DEGs)) {
        gex <- get_gex(geneID,gex.data,anno)
        # base statistics
        #get_stats(gex,"PAM50",pam50.DEGs[i])
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
                          g3.pos = 6, g3.sign = res[2,3],
                          g2.pos = 7, g2.sign = res[1,3],
                          ylim = c(-8,8),
                          ylab = "Expression (log2)", 
                          title = substitute(paste(italic(geneID)," expression in PAM50 subtypes (ERpHER2n)")))))
}
    
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

# top DEGs (high difference in )
top.DEGs <- DE.res %>% 
    filter(Her2.LumA.padj <= 0.05) %>% 
    filter(Her2.LumB.padj <= 0.05) %>% 
    filter(abs(Her2.LumA.diff) >= 1.5) %>% 
    filter(abs(Her2.LumB.diff) >= 1.5) %>% 
    dplyr::select(Her2.LumA.de, Her2.LumA.padj, Her2.LumB.de, Her2.LumB.padj) %>% rownames_to_column("Gene") %>% pull(Gene)

for (i in 1:length(top.DEGs)) {
    geneID <- mart.res %>% 
        filter(ensembl_gene_id == top.DEGs[i]) %>% 
        pull(hgnc_symbol)
    print(geneID)
}

#-----------------------plot expression-------------------------------#

# plot expression
plot.list <- list()

for (i in 1:length(top.DEGs)) {
    hgncID <- mart.res %>% 
        filter(ensembl_gene_id == top.DEGs[i]) %>% 
        pull(hgnc_symbol)
    geneID <- top.DEGs[i]
    gex <- get_gex(geneID,gex.data,anno)
    # base statistics
    get_stats(gex,"PAM50",geneID)
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
                      g3.pos = 5.5, g3.sign = res[2,3],
                      g2.pos = 6.5, g2.sign = res[1,3],
                      ylim = c(-7,7),
                      ylab = "Expression (log2)", 
                      title = paste(hgncID," (",geneID,") expression in PAM50 subtypes (ERpHER2n;",cohort,")",sep = ""))))
}

# save plots
pdf(file = paste(output.path,cohort,"_HER2n_topDEGs_expression.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()
