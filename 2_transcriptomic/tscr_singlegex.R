# Script: Single gene expression assessment (incl. ERBB2, ESR1, FGFR4)

#TODO

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
# look into the expression of selected genes
#######################################################################

# *<0.05, **<0.01 ***<0.001 ****< 0.0001 ns not significant

# list to store plots
plot.list <- list()

# for SCANB

if (cohort=="SCANB") {
    
    #ERBB2-----------------------------------------------------------------
    
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
                  lumb.pos = 7.4, lumb.sign = res[2,3],
                  luma.pos = 9, luma.sign = res[1,3],
                  ylim = c(-5,10),
                  ylab = "Expression (log2)", 
                  title = substitute(paste(italic("ERBB2")," (ENSG00000141736) expression in PAM50 subtypes (ERpHER2n)")))))
    
    #ESR1------------------------------------------------------------------
    
    esr1.gex <- get_gex("ENSG00000091831",gex.data,anno)
    # base statistics
    get_stats(esr1.gex,"PAM50","ENSG00000091831")
    # test
    res <- pair_ttest(esr1.gex, 
                      group.var = "PAM50",
                      test.var = "ENSG00000091831", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        uni_quickplot(esr1.gex,
                      group.var = "PAM50",
                      test.var = "ENSG00000091831",
                      lumb.pos = 2.6, lumb.sign = res[2,3],
                      luma.pos = 3.7, luma.sign = res[1,3],
                      ylim = c(-5,5),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("ESR1")," (ENSG00000091831) expression in PAM50 subtypes (ERpHER2n)")))))
    
    #FGFR4-----------------------------------------------------------------
    
    fgfr4.gex <- get_gex("ENSG00000160867",gex.data,anno)
    # base statistics
    get_stats(fgfr4.gex,"PAM50","ENSG00000160867")
    # test
    res <- pair_ttest(fgfr4.gex, 
                      group.var = "PAM50",
                      test.var = "ENSG00000160867", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        uni_quickplot(fgfr4.gex,
                      group.var = "PAM50",
                      test.var = "ENSG00000160867",
                      lumb.pos = 6, lumb.sign = res[2,3],
                      luma.pos = 7, luma.sign = res[1,3],
                      ylim = c(-2,8),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("FGFR4")," (ENSG00000160867) expression in PAM50 subtypes (ERpHER2n)")))))

#######################################################################

# for METABRIC (have to use HUGO gene symbold instead of ensembl ids here)
} else if (cohort=="METABRIC") {
    
    #ERBB2-----------------------------------------------------------------
    
    erbb2.gex <- get_gex("ERBB2",gex.data,anno)
    # base statistics
    get_stats(erbb2.gex,"PAM50","ERBB2")
    # test
    res <- pair_ttest(erbb2.gex, 
                      group.var = "PAM50",
                      test.var = "ERBB2", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        uni_quickplot(erbb2.gex,
                      group.var = "PAM50",
                      test.var = "ERBB2",
                      lumb.pos = 2.5, lumb.sign = res[2,3],
                      luma.pos = 3.5, luma.sign = res[1,3],
                      ylim = c(-4,4),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("ERBB2")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #ESR1------------------------------------------------------------------
    
    esr1.gex <- get_gex("ESR1",gex.data,anno)
    # base statistics
    get_stats(esr1.gex,"PAM50","ESR1")
    # test
    res <- pair_ttest(esr1.gex, 
                      group.var = "PAM50",
                      test.var = "ESR1", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        uni_quickplot(esr1.gex,
                      group.var = "PAM50",
                      test.var = "ESR1",
                      lumb.pos = 1.6, lumb.sign = res[2,3],
                      luma.pos = 2.3, luma.sign = res[1,3],
                      ylim = c(-2,3.5),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("ESR1")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #FGFR4-----------------------------------------------------------------
    
    fgfr4.gex <- get_gex("FGFR4",gex.data,anno)
    # base statistics
    get_stats(fgfr4.gex,"PAM50","FGFR4")
    # test
    res <- pair_ttest(fgfr4.gex, 
                      group.var = "PAM50",
                      test.var = "FGFR4", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        uni_quickplot(fgfr4.gex,
                      group.var = "PAM50",
                      test.var = "FGFR4",
                      lumb.pos = 4.2, lumb.sign = res[2,3],
                      luma.pos = 5.2, luma.sign = res[1,3],
                      ylim = c(-3,7),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("FGFR4")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #######################################################################
    
    
    
}
# save plots
pdf(file = paste(output.path,cohort,"_HER2n_selectedgenes.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()
