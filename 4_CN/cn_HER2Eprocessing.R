# Script: Process SCANB CN data to prepare for plotting
# Purpose: 
# 1. add genome position for probes
# 2. map the probes to genes for the driver mutation plots 

# TODO:
# - 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run (for plot naming)
cohort <- "SCANB"

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/",sep="")
dir.create(data.path)

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
#library(readxl)

################################################################################
# load data
################################################################################

# key file for gene to probe mapping
map.key <- loadRData("data/SCANB/4_CN/processed/CN_mapped_probes.RData")

# chr lengths
# get chr lengths to get genome positions of probes (excluding chr X)
chr.lengths <- as.data.frame(read.table(file = "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv", sep = '\t', header = FALSE))[1:22,] %>% 
  dplyr::rename(Chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% 
  mutate(genome = lag(genome,default = 0)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))
chr.lengths$Chr <- as.numeric(gsub('^.{3}','',chr.lengths$Chr))

# files which to process
gainloss.cn.scanb <- loadRData('data/SCANB/4_CN/raw/CN_gainloss.RData')
gainloss.cn.scanb$Chr <- as.numeric(gainloss.cn.scanb$Chr)
amp.cn.scanb <- loadRData('data/SCANB/4_CN/raw/CN_amplification.RData')
amp.cn.scanb$Chr <- as.numeric(amp.cn.scanb$Chr)

# save into file list
cnfile.list <- list("gainloss" = gainloss.cn.scanb, 
                  "amp" = amp.cn.scanb)

################################################################################
# add genome position of probes ans corresponding genes
################################################################################

for(i in c(1:length(cnfile.list))) {
  cnfile <- cnfile.list[[i]]
  # add genome position
  cnfile <- add_genomepos(cnfile,chr.lengths)
  # add genes
  # chr seq need to be in format chr1 etc.
  cnfile$Chr <- sub("^", "chr", cnfile$Chr)
  cnfile <- merge(cnfile,map.key,by=c("Chr","ProbeID","Position")) %>% 
    relocate(Gene_symbol,.after = ProbeID)
  
  # save
  save(cnfile, 
       file = paste(data.path,"CN_",names(cnfile.list)[1],"_genpos_genmap.RData",sep=""))
} 