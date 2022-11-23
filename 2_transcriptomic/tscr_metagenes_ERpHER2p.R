# Script: Metagene analysis for HER2E and ERpHER2p contrast groups in SCAN-B 

# TODO:
# no data for metabric (only erpher2n or tnbc in dataset)

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB (or METABRIC - no data)

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
library(forcats)
library(janitor)

#######################################################################
# load metagene definitions
#######################################################################

# metagene definitions
metagene.def <- as.data.frame(read_excel(
    "data/SCANB/2_transcriptomic/raw/metagene_definitions.XLSX")) %>% 
    dplyr::rename(entrezgene_id = `Entrez Gene ID`,
                  module = `Module Name`,
                  gene_symbol = `Gene symbol`)
    
#######################################################################
# 2. Cohort-specific data preprocessing including selection of
# the clinical ER+Her2- subtyped samples
#######################################################################

# for SCANB
if (cohort=="SCANB") {
    
    # load annotation data
    clin.rel4 <- as.data.frame(
        read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))
    
    # load gene anno data to convert IDs
    load("./data/SCANB/1_clinical/raw/Gene.ID.ann.Rdata")
    
    # load gex data
    load("data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata")
    
    # select subgroup data
    anno <- clin.rel4 %>% 
        filter(Follow.up.cohort==TRUE) %>% 
        filter(ER=="Positive") %>% 
        dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50) %>% 
        mutate(Group = case_when(
            HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
            HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
            HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
        filter(Group %in% 
                   c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E"))
    
    # filter to select subgroup gex data
    # modfiy ensembl ids to remove version annotation
    gex.data <- as.data.frame(genematrix_noNeg[,colnames(genematrix_noNeg) %in% anno$sampleID]) %>% 
        rownames_to_column("ensembl_gene_id") %>% 
        mutate(ensembl_gene_id = gsub("\\..*","",ensembl_gene_id)) # remove characters after dot
    
    # convert to entrez ids
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                              dataset = "hsapiens_gene_ensembl",
                              host = "http://www.ensembl.org")
    res <- getBM(filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","entrezgene_id"),
                  values = gex.data$ensembl_gene_id, 
                  mart = mart)
    save(res,file = paste(data.path,"mart_res.RData",sep="") )
    #load(paste(data.path,"mart_res.RData",sep=""))
    
    # get final metagene gex data
    gex.data <- as.data.frame(merge(gex.data, res, by="ensembl_gene_id")) %>% 
        filter(entrezgene_id %in% metagene.def$entrezgene_id) %>% 
        column_to_rownames(var="entrezgene_id") %>% 
        dplyr::select(-c(ensembl_gene_id)) %>% 
        select_if(~ !any(is.na(.)))
    
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
        drop_na(Entrez_Gene_Id) %>% 
        distinct(Entrez_Gene_Id,.keep_all = TRUE) %>% 
        column_to_rownames(var="Entrez_Gene_Id") %>% 
        dplyr::select(-c(Hugo_Symbol)) 
    
    # extract relevant variables
    anno <- anno %>% 
        mutate(ER = if_else("ERp" %in% ClinGroup,"Positive","Negative")) %>% 
        mutate(HER2 = if_else("HER2p" %in% ClinGroup,"Positive","Negative")) %>% 
        filter(ER == "Positive") %>% 
        mutate(Group = case_when(
            HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
            HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
            HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
        filter(Group %in% 
                   c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E")) %>% 
        dplyr::rename(sampleID=METABRIC_ID) # rename to match SCANB variables
    
    # filter to select subgroup gex data HERE
    gex.data <- gex.data[,colnames(gex.data) %in% anno$sampleID] %>% 
        rownames_to_column("entrezgene_id") %>% 
        filter(entrezgene_id %in% metagene.def$entrezgene_id) %>% 
        column_to_rownames(var="entrezgene_id") %>% 
        mutate_all(function(x) as.numeric(x)) %>% 
        select_if(~ !any(is.na(.))) 
    
    # exclude samples from anno without associated gex data
    anno <- anno %>% 
        filter(sampleID %in% colnames(gex.data))
}

#######################################################################
# 4. calc. score for each metagene in each sample
#######################################################################

#basal
basal.scores <- mgscore(metagene = "Basal",
                        metagene.def = metagene.def,
                        gex.data = gex.data) 

#earlyresponse
earlyresponse.scores <- mgscore(metagene = "Early_response",
                                metagene.def = metagene.def,
                                gex.data = gex.data)

metagene.scores <- merge(basal.scores,earlyresponse.scores,by=0) %>% 
    column_to_rownames(var = "Row.names")

#immuneresponse
ir.scores <- mgscore(metagene = "IR",
                     metagene.def = metagene.def,
                     gex.data = gex.data)

metagene.scores <- merge(metagene.scores,ir.scores,by=0) %>% column_to_rownames(var = "Row.names")

#lipid
lipid.scores <- mgscore(metagene = "Lipid",
                        metagene.def = metagene.def,
                        gex.data = gex.data)

metagene.scores <- merge(metagene.scores,lipid.scores,by=0) %>% column_to_rownames(var = "Row.names")

#mitoticcheckpoint
mitoticcheckpoint.scores <- mgscore(metagene = "Mitotic_checkpoint",
                                    metagene.def = metagene.def,
                                    gex.data = gex.data)

metagene.scores <- merge(metagene.scores,mitoticcheckpoint.scores,by=0) %>% column_to_rownames(var = "Row.names")


#mitoticprogression
mitoticprogression.scores <- mgscore(metagene = "Mitotic_progression",
                                     metagene.def = metagene.def,
                                     gex.data = gex.data)

metagene.scores <- merge(metagene.scores,mitoticprogression.scores,by=0) %>% column_to_rownames(var = "Row.names")


#SR
sr.scores <- mgscore(metagene = "SR",
                     metagene.def = metagene.def,
                     gex.data = gex.data)

metagene.scores <- merge(metagene.scores,sr.scores,by=0) %>% column_to_rownames(var = "Row.names")

#stroma
stroma.scores <- mgscore(metagene = "Stroma",
                         metagene.def = metagene.def,
                         gex.data = gex.data)

metagene.scores <- merge(metagene.scores,stroma.scores,by=0) %>% column_to_rownames(var = "Row.names")


#######################################################################
# 5. statistics for metagene scores between groups
#######################################################################



#mutate(Group = fct_relevel(Group,"HER2p_nonHER2E","HER2p_HER2E","HER2n_HER2E")) 

# get pvalues
mg.pvals <- data.frame()
for(i in 1:ncol(metagene.scores)) {
    mg <- colnames(metagene.scores)[i]
    res <- pair_ttest(metagene.scores,
                      anno = anno,
                      group.var = "Group",
                      test.var = mg, 
                      g1 = "HER2n_HER2E", g2 = "HER2p_nonHER2E", g3 = "HER2p_HER2E")
    mg.pvals <- rbind(mg.pvals, c(mg,res$pval[1],res$signif[1],res$pval[2],res$signif[2]))
}

# name column and set rownames
mg.pvals <- mg.pvals %>% data.table::setnames(., old = colnames(mg.pvals), 
        new = c("metagene", "HER2n_HER2E.HER2p_nonHER2E.pval", "HER2n_HER2E.HER2p_nonHER2E.signif", "HER2n_HER2E.HER2p_HER2E.pval", "HER2n_HER2E.HER2p_HER2E.signif")) %>% column_to_rownames(var="metagene")

# create group label column (ERpHER2nHER2E, ERpHER2p, ERpHER2pHER2E)
mg.anno <- merge(metagene.scores %>% rownames_to_column(var="sampleID"),anno[,c("sampleID","Group")],by="sampleID")

#######################################################################
# 5. Boxplots
#######################################################################
source("./scripts/2_transcriptomic/src/tscr_functions.R")
# plot
# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant

# list to save plots
plot.list <- list()

# SCANB
if (cohort=="SCANB") {
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "Basal",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                     g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 2, g3.sign = mg.pvals["Basal",4],
                      g2.pos = 3.5, g2.sign = mg.pvals["Basal",2],
                      ylim = c(-1.5,4),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "Basal metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                  group.var = "Group",
                  test.var = "Early_response",
                  g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                  g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                  g3.pos = 2.4, g3.sign = mg.pvals["Early_response",4],
                  g2.pos = 2.9, g2.sign = mg.pvals["Early_response",2],
                  ylim = c(-2.5,3.5),
                  ylab = "Metagene score", 
                  xlab = "HER2 subtype",
                  title = "Early_response metagene scores in HER2 subtypes (ERp)")))
                                        
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                  group.var = "Group",
                  test.var = "IR",
                  g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                  g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                  g3.pos = 2.8, g3.sign = mg.pvals["IR",4],
                  g2.pos = 3.3, g2.sign = mg.pvals["IR",2],
                  ylim = c(-2,4.1),
                  ylab = "Metagene score", 
                  xlab = "HER2 subtype",
                  title = "IR metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                  group.var = "Group",
                  test.var = "Lipid",
                  g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                  g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                  g3.pos = 2.5, g3.sign = mg.pvals["Lipid",4],
                  g2.pos = 3, g2.sign = mg.pvals["Lipid",2],
                  ylim = c(-1.8,3.5),
                  ylab = "Metagene score", 
                  xlab = "HER2 subtype",
                  title = "Lipid metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                  group.var = "Group",
                  test.var = "Mitotic_checkpoint",
                  g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                  g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                  g3.pos = 2.7, g3.sign = mg.pvals["Mitotic_checkpoint",4],
                  g2.pos = 3.2, g2.sign = mg.pvals["Mitotic_checkpoint",2],
                  ylim = c(-3,3.8),
                  ylab = "Metagene score", 
                  xlab = "HER2 subtype",
                  title = "Mitotic_checkpoint metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                  group.var = "Group",
                  test.var = "Mitotic_progression",
                  g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                  g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                  g3.pos = 2.8, g3.sign = mg.pvals["Mitotic_progression",4],
                  g2.pos = 3.4, g2.sign = mg.pvals["Mitotic_progression",2],
                  ylim = c(-3,4),
                  ylab = "Metagene score", 
                  xlab = "HER2 subtype",
                  title = "Mitotic_progression metagene scores in HER2 subtypes (ERp)")))
    
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                  group.var = "Group",
                  test.var = "SR",
                  g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                  g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                  g3.pos = 1.4, g3.sign = mg.pvals["SR",4],
                  g2.pos = 2.1, g2.sign = mg.pvals["SR",2],
                  ylim = c(-4,3),
                  ylab = "Metagene score", 
                  xlab = "HER2 subtype",
                  title = "SR metagene scores in HER2 subtypes (ERp)")))

    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                  group.var = "Group",
                  test.var = "Stroma",
                  g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                  g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                  g3.pos = 1.6, g3.sign = mg.pvals["Stroma",4],
                  g2.pos = 2.1, g2.sign = mg.pvals["Stroma",2],
                  ylim = c(-2.7,2.7),
                  ylab = "Metagene score", 
                  xlab = "HER2 subtype",
                  title = "Stroma metagene scores in HER2 subtypes (ERp)")))

#-----------------------------------------------------------------------#
    
# METABRIC
} else if(cohort=="METABRIC") {

    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "Basal",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 2, g3.sign = mg.pvals["Basal",4],
                      g2.pos = 3.5, g2.sign = mg.pvals["Basal",2],
                      ylim = c(-1.5,4),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "Basal metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "Early_response",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 2.4, g3.sign = mg.pvals["Early_response",4],
                      g2.pos = 2.9, g2.sign = mg.pvals["Early_response",2],
                      ylim = c(-2.5,3.5),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "Early_response metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "IR",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 2.8, g3.sign = mg.pvals["IR",4],
                      g2.pos = 3.3, g2.sign = mg.pvals["IR",2],
                      ylim = c(-2,4.1),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "IR metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "Lipid",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 2.5, g3.sign = mg.pvals["Lipid",4],
                      g2.pos = 3, g2.sign = mg.pvals["Lipid",2],
                      ylim = c(-1.8,3.5),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "Lipid metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "Mitotic_checkpoint",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 2.7, g3.sign = mg.pvals["Mitotic_checkpoint",4],
                      g2.pos = 3.2, g2.sign = mg.pvals["Mitotic_checkpoint",2],
                      ylim = c(-3,3.8),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "Mitotic_checkpoint metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "Mitotic_progression",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 2.8, g3.sign = mg.pvals["Mitotic_progression",4],
                      g2.pos = 3.4, g2.sign = mg.pvals["Mitotic_progression",2],
                      ylim = c(-3,4),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "Mitotic_progression metagene scores in HER2 subtypes (ERp)")))
    
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "SR",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 1.4, g3.sign = mg.pvals["SR",4],
                      g2.pos = 2.1, g2.sign = mg.pvals["SR",2],
                      ylim = c(-4,3),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "SR metagene scores in HER2 subtypes (ERp)")))
    
    #
    plot.list <- append(plot.list, list(
        three_boxplot(mg.anno,
                      group.var = "Group",
                      test.var = "Stroma",
                      g1="HER2n_HER2E",g2="HER2p_nonHER2E",g3="HER2p_HER2E", 
                      g1.col="#d334eb", g2.col="#d8b365", g3.col="#5ab4ac", 
                      g3.pos = 1.6, g3.sign = mg.pvals["Stroma",4],
                      g2.pos = 2.1, g2.sign = mg.pvals["Stroma",2],
                      ylim = c(-2.7,2.7),
                      ylab = "Metagene score", 
                      xlab = "HER2 subtype",
                      title = "Stroma metagene scores in HER2 subtypes (ERp)")))
    
    
}
#plot
pdf(file = paste(output.path,cohort,"_HER2p_metagenes.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()

###########################################################################
###########################################################################
