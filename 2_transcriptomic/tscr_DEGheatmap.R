# Script: Generate heatmaps for DEGs
#TODO: 
# exclude other group in pairwise heatmap?
# adapt group colors

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
library(amap)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for SCANB
if (cohort=="SCANB") {
    
    # load annotation data
    clin.rel4 <- as.data.frame(
        read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx")) #coltypes='text'
    
    # load gene anno data to convert IDs
    load("./data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata")
    # load gex data
    load("data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata")
    
    # select subgroup data
    anno <- clin.rel4 %>% 
        filter(Follow.up.cohort==TRUE) %>% 
        filter(NCN.PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
        filter(ER=="Positive" & HER2=="Negative") %>% 
        dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50)
    
    # filter to select subgroup gex data
    # modfiy ensembl ids to remove version annotation
    gex.data <- as.data.frame(genematrix_noNeg[,colnames(genematrix_noNeg) %in% anno$sampleID]) %>% 
        rownames_to_column("ensembl_gene_id") %>% 
        mutate(geneID = Gene.ID.ann[which(Gene.ID.ann$Gene.ID==ensembl_gene_id),]$Gene.Name) %>% 
        dplyr::select(-c(ensembl_gene_id)) %>% 
        drop_na(geneID) %>% 
        distinct(geneID,.keep_all = TRUE) %>% 
        column_to_rownames("geneID") %>% 
        select_if(~ !any(is.na(.))) # need this here because i scale/row-center
    
    # exclude samples from anno without associated gex data
    anno <- anno %>% 
        filter(sampleID %in% colnames(gex.data))
    
    # log transformed FPKM data
    gex.data <- as.data.frame(log2(gex.data + 1))
    # filter based on stdev cutoff before z-transform
    #gex.data <- gex.data %>% mutate(stdev=rowSds(as.matrix(.[colnames(gex.data)]))) %>% filter(stdev >= 0.5) %>% dplyr::select(-c(stdev)) 
    # z-transform
    gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases
    
#-----------------------------------------------------------------------#
    
} else if (cohort=="METABRIC") {
    
    # load annotation data
    load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
    
    # load gex data
    gex.data <- as.data.frame(read.table("data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt", sep="\t")) %>%
        row_to_names(row_number = 1) %>% 
        mutate_all(na_if,"") %>% 
        drop_na(Hugo_Symbol) %>% 
        distinct(Hugo_Symbol,.keep_all = TRUE) %>% 
        column_to_rownames(var="Hugo_Symbol") %>% 
        dplyr::select(-c(Entrez_Gene_Id)) 
    
    # extract relevant variables
    anno <- anno %>% 
        filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
        filter(grepl('ERpHER2n', ClinGroup)) %>% 
        dplyr::rename(sampleID=METABRIC_ID, NHG = Grade) # rename to match SCANB variables
    
    # filter to select subgroup gex data HERE
    gex.data <- gex.data[,colnames(gex.data) %in% anno$sampleID] %>% 
        mutate_all(function(x) as.numeric(x))
    
    
    # exclude samples from anno without associated gex data
    anno <- anno %>% 
        filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# Load DE results and convert ids
#######################################################################

# load DE results
load(file = "./data/SCANB/2_transcriptomic/processed/DE_results.RData")
DE.res.scanb <- DE.res
load(file = "./data/METABRIC/2_transcriptomic/processed/DE_results.RData")
DE.res.metabric <- DE.res

#######################################################################
# Heatmap parameters
#######################################################################

n.tree <- 10  # nbr of branches/clusters 
nbr.c <- 2 # color palette size = max value
breaksList = seq(-2, nbr.c, by = 0.1) # this defines the color palette range
my.link.method <- "ward.D" 
my.dist.method <- "euclidean" #"correlation" "euclidean"

mg_colors <- c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00")
# colors
my_colour = list(
    Cluster=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A"),
    PAM50 = c("Her2" = "#d334eb" , "LumA" = "#2176d5" , "LumB" = "#34c6eb"),
    NHG = c("missing"="white","1"="#fecc5c","2"="#fd8d3c","3"="#e31a1c"),
    Basal = mg_colors,
    Early_response = mg_colors,
    IR = mg_colors,
    Lipid = mg_colors,
    Mitotic_checkpoint = mg_colors,
    Mitotic_progression = mg_colors,
    SR = mg_colors,
    Stroma = mg_colors)

#######################################################################
# Heatmap 1: Top DEGs
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

# load metagene scores for hm annotation
load(paste("./data/",cohort,"/2_transcriptomic/processed/mg_anno.RData",sep = ""))

# create hm annotaiton object
hm.anno <- merge(
    mg.anno.list[[1]],anno[c("sampleID","NHG")],by="sampleID") %>% 
    relocate(NHG, .after = PAM50) %>% 
    #dplyr::select(-c(ER, HER2)) %>% 
    column_to_rownames(var = "sampleID")
    
# handle missing NHG values
hm.anno$NHG <- as.character(hm.anno$NHG)
hm.anno$NHG[hm.anno$NHG == ""] <- "missing"

# cluster samples (prev. with dist())
c1 <- cutree(hclust(
    Dist(t(top.gex),method=my.dist.method),method=my.link.method),n.tree)

# add clusters to hm annotation object 
hm.anno <- merge(hm.anno,as.data.frame(c1),by=0) %>% 
    column_to_rownames(var = "Row.names") %>% 
    dplyr::rename(Cluster=c1)

# final hm annotation 
final.hm.anno <- data.frame(
    Cluster=factor(hm.anno$Cluster),
    PAM50=factor(hm.anno$PAM50),
    NHG=factor(hm.anno$NHG),
    Basal=factor(mg_converter(hm.anno$Basal)),
    IR=factor(mg_converter(hm.anno$IR)),
    Lipid=factor(mg_converter(hm.anno$Lipid)),
    Mitotic_checkpoint=factor(mg_converter(hm.anno$Mitotic_checkpoint)),
    SR=factor(mg_converter(hm.anno$SR)),
    Stroma=factor(mg_converter(hm.anno$Stroma)))
rownames(final.hm.anno) <- rownames(hm.anno)

# create the heatmap
pheatmap(top.gex, cluster_rows=T, treeheight_row = 0,
         clustering_distance_rows = my.dist.method,
         clustering_distance_cols = my.dist.method, 
         clustering_method = my.link.method,
         show_rownames=F, show_colnames=F, 
         main=paste(
             "Top HER2E core DEGs (ERpHER2n; ",cohort,")",sep = ""), 
         cutree_cols=n.tree, legend = T,
         color=colorRampPalette(
             c("#998ec3", "#f7f7f7", "#f1a340"))(length(breaksList)),
         annotation_col=
             final.hm.anno[,
                           c("Stroma","SR","Mitotic_checkpoint",
                             "Lipid","IR","Basal","NHG","PAM50")], 
         annotation_colors=my_colour,breaks=breaksList,
         filename= paste(output.path,cohort,"_coretop_heatmap.pdf",sep = ""))


#######################################################################
# Heatmap 2: HER2E vs. LumA
#######################################################################

# load metagene scores for hm annotation
load(paste("./data/",cohort,"/2_transcriptomic/processed/mg_anno.RData",sep = ""))

# create hm annotaiton object
hm.anno <- merge(
    mg.anno.list[[1]],anno[c("sampleID","NHG")],by="sampleID") %>% 
    relocate(NHG, .after = PAM50) %>% 
    filter(PAM50 %in% c("Her2","LumA")) %>% 
    column_to_rownames(var = "sampleID")

# handle missing NHG values
hm.anno$NHG <- as.character(hm.anno$NHG)
hm.anno$NHG[hm.anno$NHG == ""] <- "missing"

# define DEG set
LumA.DEGs.scanb <- DE.res.scanb %>% 
    filter(Her2.LumA.padj <= 0.05) %>% 
    filter(abs(Her2.LumA.diff) >= 1) %>% 
    rownames_to_column("Gene") %>% 
    pull(Gene)

LumA.DEGs.metabric <- DE.res.metabric %>% 
    filter(Her2.LumA.padj <= 0.05) %>% 
    filter(abs(Her2.LumA.diff) >= 1) %>% 
    rownames_to_column("Gene") %>% 
    pull(Gene) 

# core top gex
LumA.DEGs.core <- intersect(LumA.DEGs.scanb,LumA.DEGs.metabric)

# get gex of selected genes
top.gex <- gex.data %>% 
    dplyr::select(any_of(rownames(hm.anno))) %>% 
    rownames_to_column("Gene") %>% 
    filter(Gene %in% LumA.DEGs.core) %>% 
    column_to_rownames(var="Gene") 

# cluster samples (prev. with dist())
c1 <- cutree(hclust(
    Dist(t(top.gex),method=my.dist.method),method=my.link.method),n.tree)

# add clusters to hm annotation object 
hm.anno <- merge(hm.anno,as.data.frame(c1),by=0) %>% 
    column_to_rownames(var = "Row.names") %>% 
    dplyr::rename(Cluster=c1)

# final hm annotation 
final.hm.anno <- data.frame(
    Cluster=factor(hm.anno$Cluster),
    PAM50=factor(hm.anno$PAM50),
    NHG=factor(hm.anno$NHG),
    Basal=factor(mg_converter(hm.anno$Basal)),
    IR=factor(mg_converter(hm.anno$IR)),
    Lipid=factor(mg_converter(hm.anno$Lipid)),
    Mitotic_checkpoint=factor(mg_converter(hm.anno$Mitotic_checkpoint)),
    SR=factor(mg_converter(hm.anno$SR)),
    Stroma=factor(mg_converter(hm.anno$Stroma)))
rownames(final.hm.anno) <- rownames(hm.anno)

# create the heatmap
pheatmap(top.gex, cluster_rows=T, treeheight_row = 0,
         clustering_distance_rows = my.dist.method,
         clustering_distance_cols = my.dist.method, 
         clustering_method = my.link.method,
         show_rownames=F, show_colnames=F, 
         main=paste(
             "Core HER2E-LUMA DEGs (ERpHER2n; ",cohort,")",sep = ""), 
         cutree_cols=n.tree, legend = T,
         color=colorRampPalette(
             c("#998ec3", "#f7f7f7", "#f1a340"))(length(breaksList)),
         annotation_col=
             final.hm.anno[,
                           c("Stroma","SR","Mitotic_checkpoint",
                             "Lipid","IR","Basal","NHG","PAM50")], 
         annotation_colors=my_colour,breaks=breaksList,
         filename= paste(output.path,cohort,"_LUMA_heatmap.pdf",sep = ""))

#######################################################################
# Heatmap 3: HER2E vs. LumB
#######################################################################

# load metagene scores for hm annotation
load(paste("./data/",cohort,"/2_transcriptomic/processed/mg_anno.RData",sep = ""))

# create hm annotaiton object
hm.anno <- merge(
    mg.anno.list[[1]],anno[c("sampleID","NHG")],by="sampleID") %>% 
    relocate(NHG, .after = PAM50) %>% 
    filter(PAM50 %in% c("Her2","LumB")) %>% 
    column_to_rownames(var = "sampleID")

# handle missing NHG values
hm.anno$NHG <- as.character(hm.anno$NHG)
hm.anno$NHG[hm.anno$NHG == ""] <- "missing"

# define DEG set
LumB.DEGs.scanb <- DE.res.scanb %>% 
    filter(Her2.LumB.padj <= 0.05) %>% 
    filter(abs(Her2.LumB.diff) >= 1) %>% 
    rownames_to_column("Gene") %>% 
    pull(Gene)

LumB.DEGs.metabric <- DE.res.metabric %>% 
    filter(Her2.LumB.padj <= 0.05) %>% 
    filter(abs(Her2.LumB.diff) >= 1) %>% 
    rownames_to_column("Gene") %>% 
    pull(Gene) 

# core top gex
LumB.DEGs.core <- intersect(LumB.DEGs.scanb,LumB.DEGs.metabric)

# get gex of selected genes
top.gex <- gex.data %>% 
    dplyr::select(any_of(rownames(hm.anno))) %>% 
    rownames_to_column("Gene") %>% 
    filter(Gene %in% LumB.DEGs.core) %>% 
    column_to_rownames(var="Gene") 

# cluster samples (prev. with dist())
c1 <- cutree(hclust(
    Dist(t(top.gex),method=my.dist.method),method=my.link.method),n.tree)

# add clusters to hm annotation object 
hm.anno <- merge(hm.anno,as.data.frame(c1),by=0) %>% 
    column_to_rownames(var = "Row.names") %>% 
    dplyr::rename(Cluster=c1)

# final hm annotation 
final.hm.anno <- data.frame(
    Cluster=factor(hm.anno$Cluster),
    PAM50=factor(hm.anno$PAM50),
    NHG=factor(hm.anno$NHG),
    Basal=factor(mg_converter(hm.anno$Basal)),
    IR=factor(mg_converter(hm.anno$IR)),
    Lipid=factor(mg_converter(hm.anno$Lipid)),
    Mitotic_checkpoint=factor(mg_converter(hm.anno$Mitotic_checkpoint)),
    SR=factor(mg_converter(hm.anno$SR)),
    Stroma=factor(mg_converter(hm.anno$Stroma)))
rownames(final.hm.anno) <- rownames(hm.anno)

# create the heatmap
pheatmap(top.gex, cluster_rows=T, treeheight_row = 0,
         clustering_distance_rows = my.dist.method,
         clustering_distance_cols = my.dist.method, 
         clustering_method = my.link.method,
         show_rownames=F, show_colnames=F, 
         main=paste(
             "Core HER2E-LUMB DEGs (ERpHER2n; ",cohort,")",sep = ""), 
         cutree_cols=n.tree, legend = T,
         color=colorRampPalette(
             c("#998ec3", "#f7f7f7", "#f1a340"))(length(breaksList)),
         annotation_col=
             final.hm.anno[,
                           c("Stroma","SR","Mitotic_checkpoint",
                             "Lipid","IR","Basal","NHG","PAM50")], 
         annotation_colors=my_colour,breaks=breaksList,
         filename= paste(output.path,cohort,"_LUMB_heatmap.pdf",sep = ""))
