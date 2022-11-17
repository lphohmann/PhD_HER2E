# Script: Generate heatmaps for DEGs
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
        mutate(ensembl_gene_id = gsub("\\..*","",ensembl_gene_id))  %>% # remove characters after dot 
        drop_na(ensembl_gene_id) %>% 
        distinct(ensembl_gene_id,.keep_all = TRUE) %>% 
        column_to_rownames("ensembl_gene_id") %>% 
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
# Heatmap parameters
#######################################################################

n.tree <- 10  # nbr of branches/clusters 
nbr.c <- 2 # color palette size = max value
breaksList = seq(-2, nbr.c, by = 0.1) # this defines the color palette range
my.link.method <- "ward.D" 
my.dist.method <- "euclidean" #"correlation" "euclidean"

# colors
my_colour = list(
    Cluster=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A"),
    PAM50 = c("Her2" = "red" , "LumA" = "blue" , "LumB" = "yellow"),
    NHG = c("missing"="white","1"="#E41A1C","2"="lightblue","3"="#4DAF4A"),
    Basal = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
    Early_response = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
    IR = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
    Lipid = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
    Mitotic_checkpoint = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
    Mitotic_progression = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
    SR = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
    Stroma = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"))

#######################################################################
# Heatmap data
#######################################################################

# load DE results
load(file = paste(data.path,"DE_results.RData",sep=""))

# define DEG set
DEGs <- DE.res %>% 
    filter(Her2.LumA.padj <= 0.05) %>% 
    filter(Her2.LumB.padj <= 0.05) %>% 
    filter(abs(Her2.LumA.diff) >= 1) %>% 
    filter(abs(Her2.LumB.diff) >= 1) %>% 
    dplyr::select(Her2.LumA.de, Her2.LumA.padj, Her2.LumB.de, Her2.LumB.padj) %>% rownames_to_column("Gene") %>% pull(Gene)

# get gex of selected genes
top.gex <- gex.data %>% 
    rownames_to_column("Gene") %>% 
    filter(Gene %in% DEGs) %>% 
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

#######################################################################
# Heatmap plotting
#######################################################################

# create the heatmap
pheatmap(top.gex, cluster_rows=T, treeheight_row = 0,
         clustering_distance_rows = my.dist.method,
         clustering_distance_cols = my.dist.method, 
         clustering_method = my.link.method,
         show_rownames=F, show_colnames=F, 
         main=paste(
             "Top HER2E DEGs (ERpHER2n; ",cohort,")",sep = ""), 
         cutree_cols=n.tree, legend = T,
         color=colorRampPalette(
             c("blue", "white", "yellow"))(length(breaksList)),
         annotation_col=
             final.hm.anno[,
                           c("Stroma","SR","Mitotic_checkpoint",
                             "Lipid","IR","Basal","NHG","PAM50")], 
         annotation_colors=my_colour,breaks=breaksList,
         filename= paste(output.path,cohort,"_top_heatmap.pdf",sep = ""))