# Script: FGFR4 deep dive in HER2E SCAN-B; mechanism of altered mRNA expression
# Author: Lennart Hohmann
# Date: 27.05.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")
# cohort
cohort <- "SCANB"
#-------------------
# packages
source("scripts/4_CN/src/cn_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, dplyr, stringr)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/plots/3_genomic/"
dir.create(output.path)
# for data
data.path <- "data/SCANB/3_genomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/3_genomic/processed/driver_mutations_all.RData"
infile.2 <- "./data/SCANB/3_genomic/processed/ASCAT_genelevel.RData"
infile.3 <- "./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx"

# # output paths
# plot.file <- paste0(output.path,cohort,"_FGFR4_deepdive.pdf")
# txt.file <- paste0(output.path,cohort,"_FGFR4_deepdive.txt")
# #-------------------
# # storing objects
# plot.list <- list() # object to store plots
# plot.parameters <- list() # object to store parameters to plot base R plots again later
# txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

df1 <- as.data.frame(read_excel(infile.3, sheet = "AllCodingASMD_CLPM")) 
df2 <- as.data.frame(read_excel(infile.3, sheet = "AllCodingPindel")) 

# no FGFR4 alterations in these

# scanb mut data
scanb.dmut <- loadRData(infile.1) %>% 
  mutate(PAM50 = "Her2e") %>% 
  dplyr::rename(Sample=sample, Gene=gene, Effect=variant_class) %>% 
  # split the mutation_type column in two columns (add effect column)
  mutate(Mutation_Type = sub("_$","",str_extract(Effect,"^.*?_"))) %>% 
  mutate(Effect = gsub("^.*?_","",Effect)) %>% 
  mutate(Mutation_Type = if_else(Effect=="amplified","CN",Mutation_Type)) 

# no FGFR4 only FGFR2 alterations

genes <- "FGFR4"

# for the following genes
#genes <- unique(scanb.dmut$Gene)
samples <- unique(scanb.dmut$Sample)
cna.genes <- loadRData(infile.2)
#View(cna.genes)
names(cna.genes) <- gsub("\\..*", "", names(cna.genes))
cna.genes <- cna.genes[names(cna.genes) %in% samples]
cna.driver.genes <- lapply(cna.genes, function(x) {
  #x <- cna.genes[[1]]
  #filt.x <- x[x$gene %in% genes & x$Amp==1,]
  filt.x <- x[x$gene %in% genes,]
  filt.x <- filt.x[,c("sample","gene","nMajor","nMinor","nAraw", "nBraw","size","nTot","CNA","LOH","cnnLOH","Amp","HomDel")]
  filt.x$sample <- gsub("\\..*", "", filt.x$sample)
  return(filt.x)
})
driv.amp <- do.call(rbind, cna.driver.genes)
View(driv.amp)
