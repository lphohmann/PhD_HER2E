# Script: differential expression analysis; cohorts=SCANB, METABRIC
#TODO

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
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# TODO: FORGOT I HAVE TO FILTER ROWS WITH LOW STDev

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
    res.B <- var_ttest(her2.data,luma.data)
    
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
DE.res$Her2.LumA.padj <- p.adjust(DE.res$Her2.LumA.pval, "fdr")
DE.res$Her2.LumB.padj <- p.adjust(DE.res$Her2.LumB.pval, "fdr")

# up or down (refers to the situation in lumHer2 subtype; e.g. "up" indicates a gene whose expression is upregulated in lumHer2 compared to lumB or lumA)
DE.res <- DE.res %>% mutate(Her2.LumA.de =
                                  case_when(Her2.LumA.diff <= 0 ~ "down",
                                            Her2.LumA.diff >= 0 ~ "up"))
DE.res <- DE.res %>% mutate(Her2.LumB.de =
                                  case_when(Her2.LumB.diff <= 0 ~ "down",
                                            Her2.LumB.diff >= 0 ~ "up"))

# save the file
save(DE.res,file = paste(data.path,"DE_results.RData",sep="") )

# final result 
# # Her2 vs LumA
# H2vsLA_DEGs <- results %>% filter(H2vsLA_padj <= 0.05) %>% dplyr::select(H2vsLA_de, H2vsLA_padj) %>% rownames_to_column("Gene") #%>% pull(Gene)
# H2vsLA_DEGs_up <- results %>% filter(H2vsLA_padj <= 0.05) %>% filter(H2vsLA_de == "up") %>% rownames_to_column("Gene") #%>% pull(Gene)
# H2vsLA_DEGs_down <- results %>% filter(H2vsLA_padj <= 0.05) %>% filter(H2vsLA_de == "down") %>% rownames_to_column("Gene") #%>% pull(Gene)
# 
# # Her2 vs LumB
# H2vsLB_DEGs <- results %>% filter(H2vsLB_padj <= 0.05) %>% dplyr::select(H2vsLB_de, H2vsLB_padj) %>% rownames_to_column("Gene") #%>% pull(Gene)
# H2vsLB_DEGs_up <- results %>% filter(H2vsLB_padj <= 0.05) %>% filter(H2vsLB_de == "up") %>% rownames_to_column("Gene") #%>% pull(Gene)
# H2vsLB_DEGs_down <- results %>% filter(H2vsLB_padj <= 0.05) %>% filter(H2vsLB_de == "down") %>% rownames_to_column("Gene") #%>% pull(Gene)
# 
# # compare overlap
# venn.diagram(
#     x = list(H2vsLB_DEGs$Gene, H2vsLA_DEGs$Gene),
#     category.names = c("H2vsLB" , "H2vsLA"),
#     filename = 'venn_diagramm.png',
#     output=TRUE)

# #######################################################################
# # 6. Prepare export for the functional enrichment analysis
# #######################################################################
# 
# # i think i have to strip the endings of the ids because the webtools dont recognize these identifiers
# H2vsLA_DEGs_up$Gene <- gsub("\\..*","",H2vsLA_DEGs_up$Gene)
# H2vsLA_DEGs_down$Gene <- gsub("\\..*","",H2vsLA_DEGs_down$Gene)
# H2vsLB_DEGs_up$Gene <- gsub("\\..*","",H2vsLB_DEGs_up$Gene)
# H2vsLB_DEGs_down$Gene <- gsub("\\..*","",H2vsLB_DEGs_down$Gene)
# 
# # export for functional enrichment analysis
# #write.table(H2vsLB_DEGs_up$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLB_DEGs_up.txt",sep = ""), sep="\t", row.names = FALSE)
# #write.table(H2vsLB_DEGs_down$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLB_DEGs_down.txt",sep = ""), sep="\t", row.names = FALSE)
# #write.table(H2vsLA_DEGs_up$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLA_DEGs_up.txt",sep = ""), sep="\t", row.names = FALSE)
# #write.table(H2vsLA_DEGs_down$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLA_DEGs_down.txt",sep = ""), sep="\t", row.names = FALSE)
# 
# #######################################################################
# # 7. PCA
# #######################################################################
# 
# pca_pam50anno <- sample_info %>% column_to_rownames(var="sampleID") %>% dplyr::rename(PAM50_subtype = PAM50)
# pca_data_analysis <- t(gex_data_log)
# pca_data_plotting <- merge(pca_data_analysis,pca_pam50anno,by=0) %>% column_to_rownames(var = "Row.names")
# 
# pca_res <- prcomp(pca_data_analysis) #, scale=T
# pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/pca_plot.pdf",sep =""))
# autoplot(pca_res, data=pca_data_plotting, colour="PAM50_subtype", main= "Principal component analysis of gene expression data") +
#     theme(legend.position=c(0.9,0.9))
# dev.off()
# 
# #######################################################################
# # 8. compare bonferroni vs. fdr DEGs
# #######################################################################
# 
# # load scan-b fdr DEGs
# output_dir <- "SCANB"
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/complete_DE_results.RData",sep = ""))
# scanb_DEGs <- results
# scanb_DEGs <- scanb_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,
#                                     H2vsLB_padj_fdr = H2vsLB_padj)
# 
# # load tcga fdr DEGs
# output_dir <- "TCGA"
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/complete_DE_results.RData",sep = ""))
# tcga_DEGs <- results
# tcga_DEGs <- tcga_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,
#                                   H2vsLB_padj_fdr = H2vsLB_padj)
# 
# # add bf adjustments
# # tcga
# tcga_DEGs$H2vsLA_padj_bf <- p.adjust(tcga_DEGs$H2vsLA_pvalue,"bonferroni")
# tcga_DEGs$H2vsLB_padj_bf <- p.adjust(tcga_DEGs$H2vsLB_pvalue, "bonferroni")
# # scanb
# scanb_DEGs$H2vsLA_padj_bf <- p.adjust(scanb_DEGs$H2vsLA_pvalue,"bonferroni")
# scanb_DEGs$H2vsLB_padj_bf <- p.adjust(scanb_DEGs$H2vsLB_pvalue, "bonferroni")
# 
# # DEGs
# H2vsLA_DEGs_tcga_fdr <- tcga_DEGs %>% dplyr::filter(H2vsLA_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_tcga_fdr <- tcga_DEGs %>% dplyr::filter(H2vsLB_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLA_DEGs_tcga_bf <- tcga_DEGs %>% dplyr::filter(H2vsLA_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_tcga_bf <- tcga_DEGs %>% dplyr::filter(H2vsLB_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# 
# H2vsLA_DEGs_scanb_fdr <- scanb_DEGs %>% dplyr::filter(H2vsLA_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_scanb_fdr <- scanb_DEGs %>% dplyr::filter(H2vsLB_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLA_DEGs_scanb_bf <- scanb_DEGs %>% dplyr::filter(H2vsLA_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_scanb_bf <- scanb_DEGs %>% dplyr::filter(H2vsLB_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# 
# # compare overlap
# # her2 vs lumA
# venn.diagram(
#     x = list(H2vsLA_DEGs_scanb_fdr, H2vsLA_DEGs_scanb_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_H2vsLA_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumA - SCANB)")
# 
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_fdr, H2vsLA_DEGs_tcga_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/tcga_H2vsLA_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumA - TCGA)")
# 
# # her2 vs lumB
# venn.diagram(
#     x = list(H2vsLB_DEGs_scanb_fdr, H2vsLB_DEGs_scanb_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_H2vsLB_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumB - SCANB)")
# 
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_fdr, H2vsLB_DEGs_tcga_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/tcga_H2vsLB_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumB - SCANB)")
# 
# # scanb vs tcga DEGs
# # FDR Her2-LumB
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_fdr, H2vsLB_DEGs_scanb_fdr),
#     category.names = c("TCGA DEGs" , "SCAN-B DEGs"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLB_DEGs_fdr.png",
#     main = "Shared Her2/LumB-DEGs between TCGA and SCAN-B (core set)")
# 
# # BF Her2-LumB
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_bf, H2vsLB_DEGs_scanb_bf),
#     category.names = c("TCGA" , "SCANB"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLB_DEGs_bf.png",
#     main = "Comparison Bonferroni-DEGs (Her2/LumB)")
# 
# # FDR Her2-LumA
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_fdr, H2vsLA_DEGs_scanb_fdr),
#     category.names = c("TCGA DEGs" , "SCAN-B DEGs"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLA_DEGs_fdr.png",
#     main = "Shared Her2/LumA-DEGs between TCGA and SCAN-B (core set)")
# 
# # BF Her2-LumA
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_bf, H2vsLA_DEGs_scanb_bf),
#     category.names = c("TCGA" , "SCANB"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLA_DEGs_bf.png",
#     main = "Comparison Bonferroni-DEGs (Her2/LumA)")
# 
# #######################################################################
# # 9. compare the DEGs identified using count data to the ones 
# # identified using FPKM data
# #######################################################################
# 
# # load tcga fdr DEGs
# output_dir <- "TCGA"
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/complete_DE_results.RData",sep = ""))
# tcga_fpkm_DEGs <- results
# tcga_fpkm_DEGs <- tcga_fpkm_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,
#                                   H2vsLB_padj_fdr = H2vsLB_padj)
# 
# # DEGs
# H2vsLA_DEGs_tcga_fpkm_fdr <- tcga_fpkm_DEGs %>% filter(H2vsLA_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% pull(Gene)
# H2vsLB_DEGs_tcga_fpkm_fdr <- tcga_fpkm_DEGs %>% filter(H2vsLB_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% pull(Gene)
# 
# 
# # load count DEGs
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/gex_count_luma_diffExp.RData", sep =""))
# 
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/gex_count_lumb_diffExp.RData", sep =""))
# 
# H2vsLA_DEGs_tcga_counts <- res_table_luma %>% filter(padj <= 0.05) %>% dplyr::rename(Gene = genename_luma) %>% pull(Gene)
# H2vsLB_DEGs_tcga_counts <- res_table_lumb %>% filter(padj <= 0.05) %>% dplyr::rename(Gene = genename_lumb) %>% pull(Gene)
# 
# # compare
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_counts, H2vsLB_DEGs_tcga_fpkm_fdr),
#     category.names = c("counts" , "FPKM"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/TCGA/count_vs_fpkm_H2vsLB_DEGs.png",
#     main = "TCGA: Comparison DEGs identified using count and FPKM data (Her2/LumB)")
# 
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_counts, H2vsLA_DEGs_tcga_fpkm_fdr),
#     category.names = c("counts" , "FPKM"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/TCGA/count_vs_fpkm_H2vsLA_DEGs.png",
#     main = "TCGA: Comparison DEGs identified using count and FPKM data (Her2/LumA)")
