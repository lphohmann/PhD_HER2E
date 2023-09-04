# Script: 
# test gene level CN matrices for SCANB and BASIS for singificant differences in CN gain/loss alterations

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

#-------------------------------------------------------------------------------
# packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
#library(data.table)
library(purrr)
library(readxl)
library(IRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(Repitools)
#-------------------------------------------------------------------------------
# set/create output directories
# for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)
# for data
data.path <- "data/SCANB/4_CN/processed/"
dir.create(data.path)
#-------------------------------------------------------------------------------
# input & output file paths
# input
scanb.segments <- "data/SCANB/4_CN/processed/Segment_CN_states.RData"
basis.segments <- "data/BASIS/4_CN/raw/ASCAT_CEL_Total/ASCAT_CEL_Total_EasySegments.RData"
chr.lengths <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv"
# output
plot.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.pdf"
txt.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.txt"
scanb.gene.cna <- "data/SCANB/4_CN/processed/CNA_genelevel.RData"
basis.gene.cna <- "data/BASIS/4_CN/processed/CNA_genelevel.RData"
#-------------------------------------------------------------------------------
# storing objects 
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
txt.out <- c() # object to store text output
# snippet to save plots (paste to end of script)
# pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)
# for(i in 1:length(plot.list)) { 
#   print(i)
#   print(plot.list[[i]])
# }
# dev.off()
# snippet to save text output (paste to end of script)
#writeLines(txt.out, txt.file)
#-------------------------------------------------------------------------------
# snipped to take exection time
#start.time <- Sys.time()
# CODE
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken

################################################################################
# required data
################################################################################
