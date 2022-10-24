# Script: Pathway enrichment analysis of PAM50 gene clusters

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/1_clinical/processed/",sep="")
dir.create(data.path)

#packages
#source("scripts/1_clinical/src/clin_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape)
library(openxlsx)
library(org.Hs.eg.db)
library(enrichR)

#######################################################################
# Load data and convert to entrez id
#######################################################################

# load data
clust.anno <- as.data.frame(read_table("data/SCANB/1_clinical/raw/Core_Clusters_6_Info_ann.txt"))

# convert hgnc to entrez ids 
hgnc.ids <- clust.anno %>% pull(HGNC) 
hs <- org.Hs.eg.db
res <- select(hs, 
       keys = hgnc.ids,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL") %>% dplyr::rename(HGNC=SYMBOL)
head(clust.anno)
head(res)
# merge
clust.anno <- as.data.frame(merge(clust.anno,res,by="HGNC"))

# pwe analysis
setEnrichrSite("Enrichr") 
dbs <- listEnrichrDbs() # see avaiable dbs and select 
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021","KEGG_2021_Human")
# for each cluster
enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs)