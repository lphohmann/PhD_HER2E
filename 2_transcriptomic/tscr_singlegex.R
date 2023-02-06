# Script: Single gene expression assessment (incl. ERBB2, ESR1, FGFR4)

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
source("./scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
#library(matrixStats)
library(readxl)
library(biomaRt)
library(gridExtra)
library(ggsignif)
library(janitor)
library(rstatix)

#######################################################################
# Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
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
# look into the expression of selected genes
#######################################################################

# *<0.05, **<0.01 ***<0.001 ****< 0.0001 ns not significant

# list to store plots
plot.list <- list()

# for SCANB

if (cohort=="SCANB") {
    
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
        three_boxplot(erbb2.gex,
                      group.var = "PAM50",
                      test.var = "ERBB2",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 7.4, g3.sign = res[2,3],
                      g2.pos = 9, g2.sign = res[1,3],
                      ylim = c(-5,10),
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
        three_boxplot(esr1.gex,
                      group.var = "PAM50",
                      test.var = "ESR1",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 2.6, g3.sign = res[2,3],
                      g2.pos = 3.7, g2.sign = res[1,3],
                      ylim = c(-5,5),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("ESR1")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #FGFR4-----------------------------------------------------------------
    
    fgfr4.gex <- get_gex("FGFR4",gex.data,anno)
    # save for later analysis
    save(fgfr4.gex,file = paste(data.path,"fgfr4.gex.RData",sep=""))
    
    # base statistics
    get_stats(fgfr4.gex,"PAM50","FGFR4")
    # test
    res <- pair_ttest(fgfr4.gex, 
                      group.var = "PAM50",
                      test.var = "FGFR4", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        three_boxplot(fgfr4.gex,
                      group.var = "PAM50",
                      test.var = "FGFR4",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 6, g3.sign = res[2,3],
                      g2.pos = 7, g2.sign = res[1,3],
                      ylim = c(-2,8),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("FGFR4")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #FGF19-----------------------------------------------------------------
    # 
    # fgf19.gex <- get_gex("FGF19",gex.data,anno)
    # # base statistics
    # get_stats(fgf19.gex,"PAM50","FGF19")
    # # test
    # res <- pair_ttest(fgf19.gex, 
    #                   group.var = "PAM50",
    #                   test.var = "FGF19", 
    #                   g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # # plot
    # plot.list <- append(plot.list, list(
    #   three_boxplot(fgf19.gex,
    #                 group.var = "PAM50",
    #                 test.var = "FGF19",
    #                 g1="Her2",g2="LumA",g3="LumB",
    #                 g3.pos = 6.5, g3.sign = res[2,3],
    #                 g2.pos = 7.5, g2.sign = res[1,3],
    #                 ylim = c(-0.5,1),
    #                 ylab = "Expression (log2)", 
    #                 title = substitute(paste(italic("FGF19")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #FGFR1-----------------------------------------------------------------
    
    fgfr1.gex <- get_gex("FGFR1",gex.data,anno)
    # base statistics
    get_stats(fgfr1.gex,"PAM50","FGFR1")
    # test
    res <- pair_ttest(fgfr1.gex, 
                      group.var = "PAM50",
                      test.var = "FGFR1", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        three_boxplot(fgfr1.gex,
                      group.var = "PAM50",
                      test.var = "FGFR1",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 7, g3.sign = res[2,3],
                      g2.pos = 8.5, g2.sign = res[1,3],
                      ylim = c(-5,10),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("FGFR1")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #FGFR2-----------------------------------------------------------------
    
    fgfr2.gex <- get_gex("FGFR2",gex.data,anno)
    # base statistics
    get_stats(fgfr2.gex,"PAM50","FGFR2")
    # test
    res <- pair_ttest(fgfr2.gex, 
                      group.var = "PAM50",
                      test.var = "FGFR2", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        three_boxplot(fgfr2.gex,
                      group.var = "PAM50",
                      test.var = "FGFR2",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 6.5, g3.sign = res[2,3],
                      g2.pos = 7.5, g2.sign = res[1,3],
                      ylim = c(-3.5,9),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("FGFR2")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    #FGFR3-----------------------------------------------------------------
    
    fgfr3.gex <- get_gex("FGFR3",gex.data,anno)
    # base statistics
    get_stats(fgfr3.gex,"PAM50","FGFR3")
    # test
    res <- pair_ttest(fgfr3.gex, 
                      group.var = "PAM50",
                      test.var = "FGFR3", 
                      g1 = "Her2", g2 = "LumA", g3 = "LumB")
    # plot
    plot.list <- append(plot.list, list(
        three_boxplot(fgfr3.gex,
                      group.var = "PAM50",
                      test.var = "FGFR3",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 5, g3.sign = res[2,3],
                      g2.pos = 6, g2.sign = res[1,3],
                      ylim = c(-2,8),
                      ylab = "Expression (log2)", 
                      title = substitute(paste(italic("FGFR3")," expression in PAM50 subtypes (ERpHER2n)")))))
    
    

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
        three_boxplot(erbb2.gex,
                      group.var = "PAM50",
                      test.var = "ERBB2",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 2.5, g3.sign = res[2,3],
                      g2.pos = 3.5, g2.sign = res[1,3],
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
        three_boxplot(esr1.gex,
                      group.var = "PAM50",
                      test.var = "ESR1",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 1.6, g3.sign = res[2,3],
                      g2.pos = 2.3, g2.sign = res[1,3],
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
        three_boxplot(fgfr4.gex,
                      group.var = "PAM50",
                      test.var = "FGFR4",
                      g1="Her2",g2="LumA",g3="LumB",
                      g3.pos = 4.2, g3.sign = res[2,3],
                      g2.pos = 5.2, g2.sign = res[1,3],
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
