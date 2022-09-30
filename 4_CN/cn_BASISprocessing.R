# Script: Process BASIS CN data to prepare for plotting (incl. Liftover)

# TODO:
# - convert count to % by dividing by # of samples

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run (for plot naming)
cohort <- "BASIS" 

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
# loading all CN data and convert GL freq. to %
################################################################################

load("data/BASIS/4_CN/processed/LumA_CollectedFrequencyData.RData")
cn.luma <- my.frequency.list %>% 
    mutate(CN_Gain = (CN_Gain/length(my.frequency.list$Samples))*100) %>% 
    mutate(CN_Loss = (CN_Loss/length(my.frequency.list$Samples))*100)
load("data/BASIS/4_CN/processed/LumB_CollectedFrequencyData.RData")
cn.lumb <- my.frequency.list %>% 
    mutate(CN_Gain = (CN_Gain/length(my.frequency.list$Samples))*100) %>% 
    mutate(CN_Loss = (CN_Loss/length(my.frequency.list$Samples))*100)
load("data/BASIS/4_CN/processed/ERpHER2p_CollectedFrequencyData.RData")
cn.her2p <- my.frequency.list %>% 
    mutate(CN_Gain = (CN_Gain/length(my.frequency.list$Samples))*100) %>% 
    mutate(CN_Loss = (CN_Loss/length(my.frequency.list$Samples))*100)

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

cn.her2p.postlift <- postliftprocess(hg38.pos.df,cn.her2p$fData,cn.her2p$CN_Gain,cn.her2p$CN_Loss) %>% dplyr::rename(freqloss.HER2p = Loss, freqgain.HER2p = Gain, Position = position, Chr = chr, ProbeID = reporterId) 
cn.luma.postlift <- postliftprocess(hg38.pos.df,cn.luma$fData,cn.luma$CN_Gain,cn.luma$CN_Loss) %>% dplyr::rename(freqloss.LUMA = Loss, freqgain.LUMA = Gain, Position = position, Chr = chr, ProbeID = reporterId) %>% dplyr::select(-c(Chr,Position))
cn.lumb.postlift <- postliftprocess(hg38.pos.df,cn.lumb$fData,cn.lumb$CN_Gain,cn.lumb$CN_Loss) %>% dplyr::rename(freqloss.LUMB = Loss, freqgain.LUMB = Gain, Position = position, Chr = chr, ProbeID = reporterId) %>% dplyr::select(-c(Chr,Position))

# merge
cn.basis.subtypes <- merge(cn.her2p.postlift,cn.luma.postlift, by = "ProbeID")
cn.basis.subtypes <- as.data.frame(merge(cn.basis.subtypes,cn.lumb.postlift, by = "ProbeID"))

# process to get genome positions
# 1. load the chr lengths
destfile <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv"
chr.lengths <- as.data.frame(read.table(file = destfile, sep = '\t', header = FALSE))[1:23,] %>% 
    dplyr::rename(Chr=V1,length=V2) %>% mutate(genome = cumsum(as.numeric(length))) %>% 
    mutate(genome = lag(genome,default = 0)) # ALWAYS HAS TO BE ADDED TO THE PREVIOUS CHR LENGTH so i move the colum 1 down
chr.lengths$Chr <- as.numeric(gsub("X",23,gsub('^.{3}','',chr.lengths$Chr)))

# 2. process
cn.basis.subtypes <- as.data.frame(processing(cn.basis.subtypes, chr.lengths)) %>% 
    relocate(c(Chr,Position,genome_pos), .after = ProbeID)

# save
save(cn.basis.subtypes, file = paste(data.path,"subtypes_GLmatrix_HG38",sep=""))

#head(cn.basis.subtypes)


