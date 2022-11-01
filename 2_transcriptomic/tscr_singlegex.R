# Script: Single gene expression assessment (incl. ERBB2, ESR1, ...)

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
#library(matrixStats)
library(readxl)
library(biomaRt)
library(gridExtra)
library(ggsignif)
library(janitor)

#######################################################################
# Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for SCANB
if (cohort=="SCANB") {
    
    # load annotation data
    clin.rel4 <- as.data.frame(
        read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))
    
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
# 4. for selected single genes (ERBB2/GRB7/ESR1)
#######################################################################
# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant
source("scripts/2_transcriptomic/src/tscr_functions.R")

# list to save plots
plot.list <- list()

# erbb2 #
erbb2.gex <- get_gex("ENSG00000141736",gex.data,anno)
# base statistics
get_stats(erbb2.gex,"PAM50","ENSG00000141736")
# test
res <- pair_ttest(erbb2.gex, 
            group.var = "PAM50",
            test.var = "ENSG00000141736", 
            g1 = "Her2", g2 = "LumA", g3 = "LumB")
# plot
plot.list <- append(plot.list, list(
    uni_quickplot(erbb2.gex,
              group.var = "PAM50",
              test.var = "ENSG00000141736",
              pair_ttest.res = res,
              lumb.pos = 7.4, luma.pos = 9, ylim = c(-5,10),
              ylab = "Expression (log2)", 
              title = substitute(paste(italic("ERBB2")," (ENSG00000141736) expression in PAM50 subtypes (ERpHER2n)")))))



# ESR1 gene (as part of elucidating the steroid response metagene results) #

esr1_score <- gene_score("ENSG00000091831",sg_gex_data_log,method="scale")
#View(esr1_score)
esr1_score <- merge(esr1_score,sample_info,by=0) %>% column_to_rownames(var = "Row.names")
esr1_res <- sg_ttest("ENSG00000091831",sg_gex_data_log,sample_info)
esr1_res$H2vsLA_pvalue # ****
esr1_res$H2vsLB_pvalue # ****

#plot
ggplot(esr1_score, aes(x=as.factor(PAM50),y=ENSG00000091831,fill=as.factor(PAM50))) +
    geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 subtype") +
    ylab("scaled log2 expression") +
    ylim(c(-4,4.3))+
    ggtitle("ESR1 expression (ENSG00000091831)") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="****", tip_length = 0.02, vjust=0.01, y_position = 3.5, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="****", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")
    #+ geom_text(data=as.data.frame(dplyr::count(x=sample_info, group)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -4,nudge_x = 0.3,size=5)

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/ESR1_expression.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm")

} else if (cohort=="Metabric") {
    # erbb2 #
    # HGNC symbol: ERBB2 
    erbb2_score <- gene_score("ERBB2",gex_data,method="simple")
    erbb2_res <- sg_ttest("ERBB2",gex_data,anno)
    erbb2_res$H2vsLA_pvalue # ns
    erbb2_res$H2vsLB_pvalue # ****
    
    erbb2_score <- merge(erbb2_score,anno,by=0) %>% column_to_rownames(var = "Row.names")
    
    #plot
    ggplot(erbb2_score, aes(x=as.factor(PAM50),y=as.numeric(ERBB2),fill=as.factor(PAM50))) + 
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) + 
        xlab("PAM50 subtype") + 
        ylab("scaled log2 expression") + 
        ylim(c(-4,3.5)) + 
        ggtitle("ERBB2 expression")  + 
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations="****", tip_length = 0.02, vjust=0.01, y_position = 2.6, size = 2, textsize = 15) + 
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) + 
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              legend.position = "none")
    #+geom_text(data=as.data.frame(dplyr::count(x=sample_info, group)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -4.5,nudge_x = 0.3,size=5) +
    
    ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/ERBB2_expression.pdf", sep =""),
    width = 300,
    height = 300,
    units = "mm")
    
    # ESR1 gene (as part of elucidating the steroid response metagene results) #
    
    esr1_score <- gene_score("ESR1",gex_data,method="simple")
    #View(esr1_score)
    esr1_score <- merge(esr1_score,anno,by=0) %>% column_to_rownames(var = "Row.names")
    esr1_res <- sg_ttest("ESR1",gex_data,anno)
    esr1_res$H2vsLA_pvalue # ****
    esr1_res$H2vsLB_pvalue # ****
    
    #plot
    ggplot(esr1_score, aes(x=as.factor(PAM50),y=as.numeric(ESR1),fill=as.factor(PAM50))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("PAM50 subtype") +
        ylab("scaled log2 expression") +
        ylim(c(-2.5,4.3))+
        ggtitle("ESR1 expression") +
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations="****", tip_length = 0.02, vjust=0.01, y_position = 2.5, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations="****", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              legend.position="none")
    #+ geom_text(data=as.data.frame(dplyr::count(x=sample_info, group)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -4,nudge_x = 0.3,size=5)
    
    ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/ESR1_expression.pdf", sep =""),
           width = 300,
           height = 300,
           units = "mm")
}





# save plots
pdf(file = paste(output.path,cohort,"_HER2n_metagenes.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()


# #######################################################################
# # 5. for all genes (not finished)
# #######################################################################
# 
# # i have to rerun this because it crashed in the middle of it
# 
# for(i in 1:length(gex_data_log$ensembl_gene_id)) {
#     if(i==1) {
#         res <- gene_score(gex_data_log$ensembl_gene_id[i],gex_data_log)
#     } else {
#         temp_res <- gene_score(gex_data_log$ensembl_gene_id[i],gex_data_log)
#         res <- merge(res,temp_res,by=0) %>% column_to_rownames(var = "Row.names")
#     }
# }
# 
# View(res)
# 
# # add the pam50 annotations
# all_gene_scores <- merge(res,sample_info,by=0) %>% column_to_rownames(var = "Row.names")
# 
# # save the file
# save(all_gene_scores, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/all_gene_scores.RData", sep =""))
# 
# # create pdf with all the boxplots
# pdf("gene_plots.pdf", onefile = TRUE)
# for(i in 1:length(gex_data_log$ensembl_gene_id)-1) {
#     plot <- ggplot(all_gene_scores, aes(x=as.factor(group),y=all_gene_scores[,i])) +
#         geom_boxplot(fill="slateblue",alpha=0.2) +
#         xlab("PAM50 subtype") +
#         ylab("scaled log2 expression") +
#         ggtitle(colnames(all_gene_scores)[i])
#     print(plot)
# }
# dev.off()
# 
# 
# 
# ####
# 
# #reson <- some samples in anno but not in gex
# anno <- sample_info
# Her2_cols <- anno %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="Her2") %>% pull(sampleID)
# LumA_cols <- anno %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumA") %>% pull(sampleID)
# LumB_cols <- anno %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumB") %>% pull(sampleID) 
# # extract the gex
# gex <- gex_data %>% filter(ensembl_gene_id == "ERBB2") %>% column_to_rownames(var = "ensembl_gene_id")
# 
# as.numeric(gex[1,Her2_cols])
# gex[1,LumA_cols]
# gex[1,LumB_cols]
# 
# # for Her2 vs LumA
# # equal variance check
# var.test(unlist(gex[1,Her2_cols]),unlist(gex[1,LumA_cols]), alternative = "two.sided")
# 
# H2vsLA_ttest_result <- t.test(as.numeric(gex[1,Her2_cols]),as.numeric(gex[1,LumA_cols]), var.equal = TRUE)
# t.test(gex[1,Her2_cols],gex[1,LumA_cols])
# 
# 
# 
# 
# #
# hdata <- as.numeric(gex[1,Her2_cols])
# adata <- as.numeric(gex[1,LumA_cols])
# bdata <- as.numeric(gex[1,LumB_cols])
# 
# var.test(unlist(hdata),unlist(adata), alternative = "two.sided")
# 
# H2vsLA_ttest_result <- t.test(as.numeric(gex[1,Her2_cols]),as.numeric(gex[1,LumA_cols]), var.equal = TRUE)
# t.test(gex[1,Her2_cols],gex[1,LumA_cols])