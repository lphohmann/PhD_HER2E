# Script: Process SCANB CN data to prepare for plotting
# Purpose: This script processes the SCANB HER2E CN file to match the format of the BASIS files for the other subtypes.

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

gainloss.cn.scanb <- as.data.frame((read.table(
    file = 'data/SCANB/4_CN/raw/gainloss_matrix.tsv', sep = '\t', header = TRUE)))

################################################################################
# calc. the gain/loss freq. for HER2E 
################################################################################

# get chr lengths to get genome positions of probes
#url <- "https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv"
destfile <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv"
#download.file(url, destfile)
chr.lengths <- as.data.frame(read.table(file = destfile, sep = '\t', header = FALSE))[1:23,] %>% 
    dplyr::rename(Chr=V1,length=V2) %>% mutate(genome = cumsum(as.numeric(length))) %>% 
    mutate(genome = lag(genome,default = 0)) # ALWAYS HAS TO BE ADDED TO THE PREVIOUS CHR LENGTH so i move the colum 1 down
chr.lengths$Chr <- as.numeric(gsub("X",23,gsub('^.{3}','',chr.lengths$Chr)))

# process
gainloss.cn.scanb <- as.data.frame(processing(gainloss.cn.scanb, chr.lengths))

# calc loss/gain freqs per group
gainloss.cn.scanb$freqloss.HER2E <- apply(
    gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)], 1, function(x) (
        length(which(x==-1))/ncol(gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)]))*-100) # i add a minus to make it easier for plotting

gainloss.cn.scanb$freqgain.HER2E <- apply(
    gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)], 1, function(x) (
        length(which(x==1))/ncol(gainloss.cn.scanb[,5:ncol(gainloss.cn.scanb)]))*100)

head(gainloss.cn.scanb)
# remove the indiv. sample data
gainloss.cn.scanb <- gainloss.cn.scanb %>% 
    relocate(freqgain.HER2E, .after=genome_pos) %>% 
    relocate(freqloss.HER2E, .after=genome_pos) 
gainloss.cn.scanb[,7:ncol(gainloss.cn.scanb)] <- NULL

#head(gainloss.cn.scanb)

# save
save(gainloss.cn.scanb, file = paste(data.path,"processed_gainloss_matrix",sep=""))

################################################################################
# plot the gain/loss freq. for HER2E 
################################################################################

# 
# pdf(file = paste(output.path,"HER2E_CNprofile.pdf", sep =""), height = 21.0, width = 72.0)
# 
# plot <- ggplot(gainloss.cn.scanb, aes(x=genome_pos)) + 
#     ggtitle("Frequency of gain/loss CN alterations in luminal HER2E BC") +
#     geom_line(aes(y = freqgain.HER2E, color = "red"),size=6) + 
#     geom_line(aes(y = freqloss.HER2E, color = "darkgreen"),size=6) + 
#     scale_colour_manual(name = NULL, values = c("#d72e2b","#07a109")) + # get corect colors +
#     #scale_colour_manual(name = "PAM50", values = c("#d72e2b", "#07a109", "#2176d5"), labels = c("HER2", "LUMA", "LUMB")) + # get corect colors
#     geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=2) +
#     scale_x_continuous(name = "Genome position (chromosome)",
#                        breaks = chr.lengths$genome,
#                        labels = c(as.character(1:22),"X"),
#                        limits = c(min(gainloss.cn.scanb$genome_pos),max(gainloss.cn.scanb$genome_pos)+50000000),
#                        expand = c(0, 0)) +
#     scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
#                        breaks = c(seq(-100,100,25)),
#                        labels = c(100,75,50,25,0,25,50,75,100),
#                        expand = c(0, 0),
#                        limits = c(-100,100)) +
#     theme(text=element_text(size=35),
#           legend.title = element_blank(),
#           axis.title.y = element_text(vjust = 0.5),
#           legend.position = "none") +
#     annotate(x=min(gainloss.cn.scanb$genome_pos)+30000000,y=c(-50,50), label=c("Loss","Gain"), 
#              geom="text", angle=90, hjust=0.5, size=15, colour=c("black","black")) 
# 
# 
# print(plot)
# 
# dev.off()
