# Script: Process SCANB CN data to prepare for plotting
# Purpose: This script processes the SCANB HER2E CN file to match the format of the BASIS files for the other subtypes.

# and adding the genome position of probes (in addition to the exisitng chr positions)

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
library(readxl)

################################################################################
# loading data
################################################################################

gainloss.cn.scanb <- loadRData('data/SCANB/4_CN/raw/CN_gainloss.RData')
gainloss.cn.scanb$Chr <- as.numeric(gainloss.cn.scanb$Chr)

# load BASIS files to see format
#basis.A <- loadRData("data/BASIS/4_CN/processed/LumA_CollectedFrequencyData.RData")
#basis.B <- loadRData("data/BASIS/4_CN/processed/LumB_CollectedFrequencyData.RData")
#basis.H <- loadRData("data/BASIS/4_CN/processed/ERpHER2p_CollectedFrequencyData.RData")

################################################################################
# add genome position of probes
################################################################################

# get chr lengths to get genome positions of probes (excluding chr X)
chr.lengths <- as.data.frame(read.table(file = "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv", sep = '\t', header = FALSE))[1:22,] %>% 
  dplyr::rename(Chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% 
  mutate(genome = lag(genome,default = 0)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))
chr.lengths$Chr <- as.numeric(gsub('^.{3}','',chr.lengths$Chr))

gainloss.cn.scanb <- add_genomepos(gainloss.cn.scanb,chr.lengths)

################################################################################
# calc. the gain/loss freq. for HER2E 
################################################################################

# calc loss/gain freqs per group
gainloss.cn.scanb$freqloss.HER2E <- apply(
    gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)], 1, function(x) (
        length(which(x==-1))/ncol(gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)]))*-100) # i add a minus to make it easier for plotting

gainloss.cn.scanb$freqgain.HER2E <- apply(
    gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)], 1, function(x) (
        length(which(x==1))/ncol(gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)]))*100)

#head(gainloss.cn.scanb)

# remove the indiv. sample data
gainloss.cn.scanb <- gainloss.cn.scanb %>% 
    relocate(freqgain.HER2E, .after=Genome_pos) %>% 
    relocate(freqloss.HER2E, .after=Genome_pos) 
gainloss.cn.scanb[,7:ncol(gainloss.cn.scanb)] <- NULL

#head(gainloss.cn.scanb)

# save
save(gainloss.cn.scanb, 
     file = "data/SCANB/4_CN/processed/CN_gainloss_frequencies.RData")


