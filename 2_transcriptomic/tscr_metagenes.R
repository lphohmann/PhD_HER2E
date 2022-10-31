# Script: Metagene analysis in SCAN-B and METABRIC

# TODO:
# try once with scaling MB data and without and compare output

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
        mutate(ensembl_gene_id = gsub("\\..*","",ensembl_gene_id)) # remove characters after dot
    
    # convert to entrez ids
    # mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
    #                          dataset = "hsapiens_gene_ensembl",
    #                          host = "http://www.ensembl.org")
    # res <- getBM(filters = "ensembl_gene_id",
    #              attributes = c("ensembl_gene_id","entrezgene_id"),
    #              values = gex.data$ensembl_gene_id, 
    #              mart = mart)
    #save(res,file = paste(data.path,"mart_res.RData",sep="") )
    load(paste(data.path,"mart_res.RData",sep=""))
    
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
        filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
        filter(grepl('ERpHER2n', ClinGroup)) %>% 
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


# check if genes were not in the scanb data
setdiff(metagene.def$entrezgene_id,rownames(gex.data)) 

#METABRIC: 125 -> ; 9370 -> ; 1520 -> ; 3507 -> "IGHM" ; 28755 <- "TRAC"
#SCANB: 3507 -> "IGHM" ; 28755 <- "TRAC"

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

metagene.scores <- merge(basal.scores,earlyresponse.scores,by=0) %>% column_to_rownames(var = "Row.names")

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

# get pvalues
mg.pvals <- mgtest(metagene.scores,anno)

# save as a summarized object
obj.anno <- anno %>% dplyr::select(c(sampleID, PAM50))
mg.scores <- metagene.scores %>% rownames_to_column(var="sampleID")
mg.anno <- as.data.frame(merge(obj.anno,mg.scores,by="sampleID"))

#mg.anno.list <- list(mg.anno, mg.pvals)
#save(mg.anno.list,file = paste(data.path,"mg_anno.RData",sep=""))

# round
mg.pvals <- round(mg.pvals, digits = 5)

#######################################################################
# 5. Boxplots
#######################################################################

# SCANB
if (cohort=="SCANB") {
    

# list to save plots
plot.list <- list()

# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant
mg.pvals["Basal",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Basal", 
    lumb.sig = "****", lumb.pos = 1, 
    luma.sig = "****", luma.pos = 2.9, 
    c(-1.5,4))))

mg.pvals["Early_response",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Early_response", 
    lumb.sig = "ns", lumb.pos = 2.5, 
    luma.sig = "****", luma.pos = 3.6, 
    c(-2,4.6))))

mg.pvals["IR",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "IR", 
    lumb.sig = "****", lumb.pos = 4, 
    luma.sig = "****", luma.pos = 4.6, 
    c(-2,5.6))))

mg.pvals["Lipid",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Lipid", 
    lumb.sig = "ns", lumb.pos = 2.2, 
    luma.sig = "****",  luma.pos = 3.1, 
    c(-1.5,4))))

mg.pvals["Mitotic_checkpoint",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Mitotic_checkpoint", 
    lumb.sig = "ns", lumb.pos = 3.7, 
    luma.sig = "****",  luma.pos = 4.2, 
    c(-2,5))))

mg.pvals["Mitotic_progression",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Mitotic_progression", 
    lumb.sig = "ns", lumb.pos = 3.6, 
    luma.sig = "****",  luma.pos = 4.1, 
    c(-2,5))))

mg.pvals["SR",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "SR", 
    lumb.sig = "****", lumb.pos = 1.8, 
    luma.sig = "****", luma.pos = 2.3, 
    c(-1.5,4))))

mg.pvals["Stroma",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Stroma", 
    lumb.sig = "****", lumb.pos = 2, 
    luma.sig = "ns", luma.pos = 2.5, 
    c(-3.5,3.5))))

#plot
pdf(file = paste(output.path,cohort,"_HER2n_metagenes.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()


#-----------------------------------------------------------------------#

# METABRIC
} else if (cohort=="METABRIC") {
    
# list to save plots
plot.list <- list()

# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant
mg.pvals["Basal",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Basal", 
    lumb.sig = "ns", lumb.pos = 0.8, 
    luma.sig = "****", luma.pos = 2.1, 
    c(-1.5,4))))

mg.pvals["Early_response",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Early_response", 
    lumb.sig = "*", lumb.pos = 2.2, 
    luma.sig = "****", luma.pos = 2.7, 
    c(-2,4.6))))

mg.pvals["IR",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "IR", 
    lumb.sig = "*", lumb.pos = 3, 
    luma.sig = "**", luma.pos = 3.6, 
    c(-2,5.6))))

mg.pvals["Lipid",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Lipid", 
    lumb.sig = "*", lumb.pos = 2.5, 
    luma.sig = "****",  luma.pos = 3.3, 
    c(-1.5,4))))

mg.pvals["Mitotic_checkpoint",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Mitotic_checkpoint", 
    lumb.sig = "ns", lumb.pos = 2.2, 
    luma.sig = "****",  luma.pos = 2.7, 
    c(-2,5))))

mg.pvals["Mitotic_progression",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Mitotic_progression", 
    lumb.sig = "ns", lumb.pos = 2.5, 
    luma.sig = "****",  luma.pos = 3, 
    c(-2,5))))

mg.pvals["SR",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "SR", 
    lumb.sig = "****", lumb.pos = 2.1, 
    luma.sig = "****", luma.pos = 2.6, 
    c(-1.5,4))))

mg.pvals["Stroma",]
plot.list <- append(plot.list, list(quickplot(
    mg.anno, "Stroma", 
    lumb.sig = "ns", lumb.pos = 2, 
    luma.sig = "***", luma.pos = 2.5, 
    c(-3.5,3.5))))

#plot
pdf(file = paste(output.path,cohort,"_HER2n_metagenes.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()

}