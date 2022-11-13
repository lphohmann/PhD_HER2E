# Script: check if pam50 genes are in the DEGs
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
#load data
#######################################################################
load(file = paste(data.path,"DE_results.RData",sep=""))

pam50.genes <- as.data.frame(read_table("data/SCANB/1_clinical/raw/Spiral_SRIQ_PAM50_Gex_Clusters_6_Info_ann.txt")) %>% pull(HGNC)

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
# Check how many PAM50 genes are DE 
#######################################################################

# DEGs
DEGs <- DE.res %>% 
    filter(Her2.LumA.padj <= 0.05) %>% 
    filter(Her2.LumB.padj <= 0.05) %>% 
    dplyr::select(Her2.LumA.de, Her2.LumA.padj, Her2.LumB.de, Her2.LumB.padj) %>% rownames_to_column("Gene")

# convert to entrez ids
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "http://www.ensembl.org")

res <- getBM(filters = "ensembl_gene_id",
             attributes = c("ensembl_gene_id","hgnc_symbol"),
             values = DEGs[,"Gene"], 
             mart = mart)

pam50.DEGs <- res %>% 
    filter(hgnc_symbol %in% pam50.genes) %>% 
    pull(hgnc_symbol)

#######################################################################
# plot expression
#######################################################################

# plot expression
plot.list <- list()


if (cohort == "SCANB") {
    # MUSS HIER DIE ENSEMBL IDS VON RES NEHMEN, DAS ANDERE DINS DIE GENE SYMBOLS 
    geneID <- res %>% 
        filter(hgnc_symbol == pam50.DEGs[i]) %>% 
        pull(ensembl_gene_id)
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
                              test.var = geneID
                              g1="Her2",g2="LumA",g3="LumB",
                              g3.pos = 6, g3.sign = res[2,3],
                              g2.pos = 7, g2.sign = res[1,3],
                              ylim = c(-8,-8),
                              ylab = "Expression (log2)", 
                              title = substitute(paste(italic(geneID)," expression in PAM50 subtypes (ERpHER2n)")))))
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
                          test.var = geneID
                          g1="Her2",g2="LumA",g3="LumB",
                          g3.pos = 6, g3.sign = res[2,3],
                          g2.pos = 7, g2.sign = res[1,3],
                          ylim = c(-8,-8),
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
