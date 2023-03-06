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
plot.file <- paste(output.path,cohort,"_HER2n_umap.pdf",sep = "")

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
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

top.umap.gex$PAM50 <- factor(top.umap.gex$PAM50,levels = c("LumA","LumB","Her2"))

# run umap
umap <- umap(top.umap.gex[1:ncol(top.umap.gex)-1],n_neighbors=50)
df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 PAM50 = top.umap.gex["PAM50"])

# plot
plot <- ggplot(df %>% arrange(PAM50), aes(x, y, colour = PAM50)) + #,order=
  geom_point(alpha=0.7, size=2.5) +
  theme(axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 35),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 35),
        plot.title = element_text(size=25)) +
  scale_color_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("Her2","LumA","LumB"))) +
  ggtitle(paste("UMAP based on top core DEGs (n genes=",ncol(top.umap.gex)-1,"; ",cohort,")",sep=""))

# plot
plot.list <- append(plot.list, list(plot))

#######################################################################
# UMAP with all core DEGs
#######################################################################

# define DEG set
DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene)

DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 

# core gex
DEGs.core <- intersect(DEGs.scanb,DEGs.metabric)

# get gex of selected genes
deg.gex <- gex.data %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% DEGs.core) %>% 
  column_to_rownames(var="Gene") 

#---------------------------------------------------------------------#

# merge to one df to have pam50 annotation for each sample in the right order
deg.gex <- as.data.frame(t(deg.gex)) %>% rownames_to_column(var="sampleID")
umap.gex <- merge(deg.gex,anno[c("PAM50","sampleID")],by="sampleID") %>% column_to_rownames(var="sampleID") 

umap.gex$PAM50 <- factor(umap.gex$PAM50,levels = c("LumA","LumB","Her2"))

# run umap
umap <- umap(umap.gex[1:ncol(umap.gex)-1],n_neighbors=50)
df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 PAM50 = umap.gex["PAM50"])

# plot
plot <- ggplot(df %>% arrange(PAM50), aes(x, y, colour = PAM50)) + #,order=
  geom_point(alpha=0.7, size=2.5) +
  theme(axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 35),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 35),
        plot.title = element_text(size=25)) +
  scale_color_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                     c("Her2","LumA","LumB"))) +
  ggtitle(paste("UMAP based on core DEGs (n genes=",ncol(umap.gex)-1,"; ",cohort,")",sep=""))

# plot
plot.list <- append(plot.list, list(plot))

#######################################################################
# UMAP with all genes
#######################################################################

# merge to one df to have pam50 annotation for each sample in the right order
all.gex <- as.data.frame(t(gex.data)) %>% rownames_to_column(var="sampleID")
umap.gex <- merge(all.gex,anno[c("PAM50","sampleID")],by="sampleID") %>% column_to_rownames(var="sampleID") 

umap.gex$PAM50 <- factor(umap.gex$PAM50,levels = c("LumA","LumB","Her2"))

# run umap
umap <- umap(umap.gex[1:ncol(umap.gex)-1],n_neighbors=50)
df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 PAM50 = umap.gex["PAM50"])

# plot
plot <- ggplot(df %>% arrange(PAM50), aes(x, y, colour = PAM50)) + #,order=
  geom_point(alpha=0.7, size=2.5) +
  theme(axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 35),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 35),
        plot.title = element_text(size=25)) +
  scale_color_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                     c("Her2","LumA","LumB"))) +
  ggtitle(paste("UMAP based on all genes (n genes=",ncol(umap.gex)-1,"; ",cohort,")",sep=""))

# plot
plot.list <- append(plot.list, list(plot))

#######################################################################

# save plots
pdf(file = plot.file, 
    onefile = TRUE, width = 14.8, height = 10.5) 
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()
