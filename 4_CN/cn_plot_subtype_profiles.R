# Script: Copy number alterations in intrinsic subtypes (SCANB + BASIS)

# TODO:
# - 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run (for plot naming)
cohort <- "COMBINED" 

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(rtracklayer)
library(liftOver)
library(GenomicRanges)
library(rtracklayer)
#library(R.utils)
library(Repitools)

################################################################################
# loading all CN data
################################################################################

load("data/BASIS/4_CN/processed/LumA_CollectedFrequencyData.RData")
cn.luma <- my.frequency.list
load("data/BASIS/4_CN/processed/LumB_CollectedFrequencyData.RData")
cn.lumb <- my.frequency.list
load("data/BASIS/4_CN/processed/ERpHER2p_CollectedFrequencyData.RData")
cn.her2p <- my.frequency.list
load("data/BASIS/4_CN/processed/")
cn.her2e <- my.frequency.list

# l <- list(cn.luma,cn.lumb,cn.her2p)
# for (i in 1:3) {
#     str(l[i])
# }

################################################################################
# liftover of BASIS SNP-A positions to build 38 which SCAN-B is on.
################################################################################

# get the right chain file
#url <- "https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz"
destfile <- "data/BASIS/4_CN/raw/hg19xToHg38.over.chain.gz"
#download.file(url, destfile)
#gunzip(destfile)
#ch <- import.chain(gsub('.{3}$', '', destfile))
ch <- import.chain("data/BASIS/4_CN/raw/hg19ToHg38.over.chain")

# get df in right format to convert to granges object (only have to do that once for all files)
df <- cn.luma$fData %>% 
    mutate(end = centerPosition) %>% 
    dplyr::rename(start = centerPosition) %>% mutate(chromosome = paste("chr",as.character(chromosome),sep = "")) %>% 
    mutate(chromosome = if_else(chromosome == "chr23","chrX",chromosome)) # change name to same as in chain file

# make granges object
gr.df <- makeGRangesFromDataFrame(df,seqnames.field = "chromosome",start.field = "start",end.field = "end", keep.extra.columns = TRUE,ignore.strand = TRUE) #starts.in.df.are.0based=FALSE # think not

# liftover
hg38.pos <- liftOver(gr.df,ch) 
hg38.pos <- unlist(hg38.pos) # convert to granges object

# compare to see if the coordinates changed -> it worked :D
df <- annoGR2DF(gr.df)
#df %>% filter(reporterId == "SNP_A-8575395") # 564477
hg38.pos.df <- annoGR2DF(hg38.pos) %>% dplyr::select(-c("end","width")) 
hg38.pos.df$chr <- as.numeric(gsub("X",23,gsub('^.{3}','',hg38.pos.df$chr)))
#hg38.pos.df %>% filter(reporterId == "SNP_A-8575395") # 629097

################################################################################
# use liftover results to convert positions
################################################################################

cn.her2p.postlift <- postliftprocess(hg38.pos.df,cn.her2p$fData,cn.her2p$CN_Gain,cn.her2p$CN_Loss)
cn.luma.postlift <- postliftprocess(hg38.pos.df,cn.luma$fData,cn.luma$CN_Gain,cn.luma$CN_Loss)
cn.lumb.postlift <- postliftprocess(hg38.pos.df,cn.lumb$fData,cn.lumb$CN_Gain,cn.lumb$CN_Loss)

