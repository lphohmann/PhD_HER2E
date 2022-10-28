# Script: Metagene analysis for HER2E and ERpHER2p contrast groups in SCAN-B 

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

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
    
# load annotation data
clin.rel4 <- as.data.frame(
    read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))

# load gex data
load("data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata")

# select subgroup data
anno <- clin.rel4 %>% 
    filter(Follow.up.cohort==TRUE) %>% 
    filter(ER=="Positive") %>% 
    dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50)

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

# get final metagene gex data
gex.data <- as.data.frame(merge(gex.data, res, by="ensembl_gene_id")) %>% 
    filter(entrezgene_id %in% metagene.def$entrezgene_id) %>% 
    column_to_rownames(var="entrezgene_id") %>% 
    dplyr::select(-c(ensembl_gene_id)) %>% 
    select_if(~ !any(is.na(.)))

# exclude samples from anno without associated gex data
anno <- anno %>% 
    filter(sampleID %in% colnames(gex.data))

setdiff(metagene.def$entrezgene_id,rownames(gex.data)) # "IGHM" "TRAC"

#######################################################################
# 3. gex data processing 
#######################################################################

# log transformed FPKM data
gex.data <- as.data.frame(log2(gex.data + 1))
# z-transform
gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases

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
# create summarized object
#######################################################################

# create group label column (ERpHER2nHER2E, ERpHER2p, ERpHER2pHER2E)
obj.anno <- anno %>% 
    dplyr::select(c(sampleID, PAM50, ER, HER2))
mg.scores <- metagene.scores %>% 
    rownames_to_column(var="sampleID")
mg.anno <- as.data.frame(merge(obj.anno,mg.scores,by="sampleID"))

# data prep for plotting
mg.anno <- mg.anno %>% 
    mutate(Group = case_when(
        HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
        HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
        HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
    mutate(Group = fct_relevel(Group,"HER2p_nonHER2E","HER2p_HER2E","HER2n_HER2E")) 

#######################################################################
# 5. statistics for metagene scores between groups
#######################################################################

# get pvalues 
mg.pvals <- her2p_mgtest(metagene.scores,mg.anno)
# round
mg.pvals <- round(mg.pvals, digits = 5)

#######################################################################
# 5. compare scaled metagene scores between groups
#######################################################################

# plot
# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant

plot.list <- list()

#
mg.pvals["Basal",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"Basal","**",2,"****",2.5, c(-1.5,3))))

#
mg.pvals["Early_response",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"Early_response","ns",2,"ns",2.3, c(-2.5,2.6))))
                                    
#
mg.pvals["IR",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"IR","ns",3.3,"*",3.7, c(-2,4.1))))

#
mg.pvals["Lipid",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"Lipid","ns",2,"****",2.5, c(-1.8,3))))

#
mg.pvals["Mitotic_checkpoint",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"Mitotic_checkpoint","ns",3,"****",3.4, c(-2,3.8))))

#
mg.pvals["Mitotic_progression",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"Mitotic_progression","ns",3.2,"****",3.6, c(-2,4))))

#
mg.pvals["SR",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"SR","ns",1.1,"****",1.5, c(-4,1.9))))

#
mg.pvals["Stroma",]
plot.list <- append(plot.list, list(
    her2p_quickplot(mg.anno,"Stroma","ns",1.5,"ns",1.9, c(-2.7,2.3))))

#plot
pdf(file = paste(output.path,cohort,"_HER2p_metagenes.pdf", sep=""), 
    onefile = TRUE, width = 15, height = 15) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()

###########################################################################
###########################################################################
