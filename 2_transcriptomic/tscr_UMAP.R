# Script: Generate UMAP

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

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_heatmaps.pdf",sep = "")

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#library(Hmisc)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)
library(amap)
library(umap)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
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
  anno <- loadRData("data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
    dplyr::rename(sampleID=METABRIC_ID,NHG=Grade) # rename to match SCANB variables
  
  # load and select subgroup data
  gex.data <- metabric_gex_load("./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt",ID.type = "Hugo_Symbol") %>% 
    dplyr::select(any_of(anno$sampleID)) %>% 
    mutate_all(function(x) as.numeric(x))
  
  # exclude samples from anno without associated gex data
  anno <- anno %>% 
    filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# Load DE results 
#######################################################################

# load DE results
load(file = "./data/SCANB/2_transcriptomic/processed/DE_results.RData")
DE.res.scanb <- DE.res
load(file = "./data/METABRIC/2_transcriptomic/processed/DE_results.RData")
DE.res.metabric <- DE.res

#######################################################################
# UMAP with top DEGs
#######################################################################

# define DEG set
top.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  filter(abs(Her2.LumA.diff) >= 1) %>% 
  filter(abs(Her2.LumB.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% pull(Gene)

top.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  filter(abs(Her2.LumA.diff) >= 1) %>% 
  filter(abs(Her2.LumB.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 

# core top gex
top.DEGs.core <- intersect(top.DEGs.scanb,top.DEGs.metabric)

# get gex of selected genes
top.gex <- gex.data %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% top.DEGs.core) %>% 
  column_to_rownames(var="Gene") 

#---------------------------------------------------------------------#

# merge to one df to have pam50 annotation for each sample in the right order
top.gex <- as.data.frame(t(top.gex)) %>% rownames_to_column(var="sampleID")
top.umap.gex <- merge(top.gex,anno[c("PAM50","sampleID")],by="sampleID") %>% column_to_rownames(var="sampleID") 

# run umap
umap <- umap(top.umap.gex[1:ncol(top.umap.gex)-1])
df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 PAM50 = top.umap.gex["PAM50"])

# plot
ggplot(df, aes(x, y, colour = PAM50)) +
  geom_point()
