# Script: differential expression analysis; cohorts=SCANB, METABRIC
#TODO: 
# Filter based on diff in expression absolute (add column)

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "METABRIC" # SCANB or METABRIC

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
        dplyr::rename(sampleID=METABRIC_ID) # rename to match SCANB variables
    
    # filter to select subgroup gex data HERE
    gex.data <- gex.data[,colnames(gex.data) %in% anno$sampleID] %>% 
        mutate_all(function(x) as.numeric(x))
    
    # exclude samples from anno without associated gex data
    anno <- anno %>% 
        filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# 4. Perform t-test for each gene
#######################################################################
#source("./scripts/2_transcriptomic/src/tscr_functions.R")

# sampleIDs per pam50 subtype
Her2.cols <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA.cols <- anno %>% filter(PAM50=="LumA") %>% pull(sampleID)
LumB.cols <- anno %>% filter(PAM50=="LumB") %>% pull(sampleID)

# initialize vector of stored p-values and expression differences (based on mean comparison)
Her2.LumA.pval <- rep(NA,nrow(gex.data))
Her2.LumB.pval <- rep(NA,nrow(gex.data))
Her2.LumA.diff <- rep(NA,nrow(gex.data))
Her2.LumB.diff <- rep(NA,nrow(gex.data))

# loop through the genes
pb = txtProgressBar(min = 0, max = nrow(gex.data), initial = 0, style = 3)
for (i in 1:nrow(gex.data)) { #nrow(gex.data)
    setTxtProgressBar(pb,i)
    
    # set vars
    her2.data <- as.numeric(gex.data[i,Her2.cols])
    luma.data <- as.numeric(gex.data[i,LumA.cols])
    lumb.data <- as.numeric(gex.data[i,LumB.cols])
    
    # for Her2 vs LumA
    res.A <- var_ttest(her2.data,luma.data)
    
    # save results
    Her2.LumA.pval[i] <- res.A$p.value
    Her2.LumA.diff[i] <- res.A$estimate[1]-res.A$estimate[2]
    
    # for Her2 vs LumB
    res.B <- var_ttest(her2.data,lumb.data)
    
    # save results
    Her2.LumB.pval[i] <- res.B$p.value
    Her2.LumB.diff[i] <- res.B$estimate[1]-res.B$estimate[2]
    close(pb)
}

# Process the output

# create the final output
DE.res <- gex.data %>% add_column(
    Her2.LumA.pval = Her2.LumA.pval,
    Her2.LumB.pval = Her2.LumB.pval,
    Her2.LumA.diff = Her2.LumA.diff, 
    Her2.LumB.diff = Her2.LumB.diff) %>% 
    dplyr::select(Her2.LumA.pval,Her2.LumB.pval,
                  Her2.LumA.diff,Her2.LumB.diff)

# adjust p value 
DE.res$Her2.LumA.padj <- p.adjust(DE.res$Her2.LumA.pval, "bonferroni")
DE.res$Her2.LumB.padj <- p.adjust(DE.res$Her2.LumB.pval, "bonferroni")

# up or down (refers to the situation in lumHer2 subtype; e.g. "up" indicates a gene whose expression is upregulated in lumHer2 compared to lumB or lumA)
DE.res <- DE.res %>% 
    mutate(Her2.LumA.de = case_when(
        Her2.LumA.diff <= 0 ~ "down", 
        Her2.LumA.diff >= 0 ~ "up")) %>% 
    mutate(Her2.LumB.de = case_when(
        Her2.LumB.diff <= 0 ~ "down", 
        Her2.LumB.diff >= 0 ~ "up"))

# save the file
save(DE.res,file = paste(data.path,"DE_results.RData",sep="") )
load(file = paste(data.path,"DE_results.RData",sep=""))