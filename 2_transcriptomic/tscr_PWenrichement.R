# Script: Pathway enrichment analysis for DEGs

# TODO

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" 

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# output filenames
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_pwenrichment_plots.pdf", sep="")
xlsx.list <- list() # object to store dataframes -> excel tables 
xlsx.file <- paste("output/supplementary_data/HER2n_pwenrichment.xlsx", sep="")

#packages
source("./scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape)
library(openxlsx)
library(enrichR)
library(janitor)

#######################################################################

# load data
DE.res.scanb <- loadRData("./data/SCANB/2_transcriptomic/processed/DE_results.RData")
DE.res.metabric <- loadRData("./data/METABRIC/2_transcriptomic/processed/DE_results.RData")

# set database parameters
setEnrichrSite("Enrichr") 
dbs <- listEnrichrDbs() # see avaiable dbs and select 
dbs <- c("KEGG_2021_Human") #"GO_Molecular_Function_2021", "GO_Biological_Process_2021

#######################################################################
# Top DEGs
#######################################################################

# define DEG set
top.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  filter(abs(Her2.LumA.diff) >= 1) %>% 
  filter(abs(Her2.LumB.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% pull(Gene)

top.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  filter(abs(Her2.LumA.diff) >= 1) %>% 
  filter(abs(Her2.LumB.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 

# core top gex
top.DEGs.core <- intersect(top.DEGs.scanb,top.DEGs.metabric)

#---------------------------------------------------------------------#

#assign("res", enrichr(top.DEGs.core, dbs)) 

res <- as.data.frame(enrichr(top.DEGs.core, dbs)[[1]])

plot.list <- append(plot.list, list(
  plotEnrich(res, showTerms = 5, numChar = 60, title = "Top core DEGs") +
    theme(text = element_text(size = 20))))

xlsx.list <- append(xlsx.list, list("Top_coreDEGs" = res))

#######################################################################
# HER2E vs. LumA
#######################################################################

# define DEG set
LumA.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(abs(Her2.LumA.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene)

LumA.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(abs(Her2.LumA.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene) 

# core top gex
LumA.DEGs.core <- intersect(LumA.DEGs.scanb,LumA.DEGs.metabric)

#---------------------------------------------------------------------#

res <- as.data.frame(enrichr(LumA.DEGs.core, dbs)[[1]])

plot.list <- append(plot.list, list(
  plotEnrich(res, showTerms = 5, numChar = 60, y = "Count", orderBy = "P.value", title = "LumA core DEGs") +
    theme(text = element_text(size = 20))))

xlsx.list <- append(xlsx.list, list("LumA_coreDEGs" = res))

#######################################################################
# HER2E vs. LumB
#######################################################################

# define DEG set
LumB.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  filter(abs(Her2.LumB.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene)

LumB.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  filter(abs(Her2.LumB.diff) >= 1) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene) 

# core top gex
LumB.DEGs.core <- intersect(LumB.DEGs.scanb,LumB.DEGs.metabric)

#---------------------------------------------------------------------#

res <- as.data.frame(enrichr(LumB.DEGs.core, dbs)[[1]])

plot.list <- append(plot.list, list(
  plotEnrich(res, showTerms = 5, numChar = 60, y = "Count", orderBy = "P.value", title = "LumB core DEGs") +
    theme(text = element_text(size = 20))))

xlsx.list <- append(xlsx.list, list("LumB_coreDEGs" = res))

#######################################################################
# save plots and excel file
#######################################################################

# save excel file
write.xlsx(xlsx.list, file = xlsx.file)

# save plots
pdf(file = plot.file, onefile = TRUE, width = 20, height = 5.5) 

for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}

dev.off()
