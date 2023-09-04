# Script: Plot gain/loss frequencies over the whole genome for all subtypes

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
library(plyr)

###############################################################################
# loading data
###############################################################################

# SCANB (HER2E)
cn.scanb <- loadRData(file = "data/SCANB/4_CN/processed/CN_gainloss_frequencies.RData")

# BASIS (HER2p, LUMA, LUMB)
cn.basis <- loadRData(file = "data/BASIS/4_CN/processed/CN_gainloss_frequencies.RData")

# convert loss to negative values like in HER2E
cn.basis[c("freqloss.HER2p","freqloss.LUMA","freqloss.LUMB")] <- lapply(cn.basis[c("freqloss.HER2p","freqloss.LUMA","freqloss.LUMB")], FUN = function(x) {x*-1})

# merge
cn.data <- rbind.fill(cn.basis,cn.scanb)
cn.data <- cn.data[order(cn.data$Genome_pos),]

# chr lengths file
chr.lengths <- as.data.frame(read.table(file = "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv", sep = '\t', header = FALSE))[1:23,] %>% 
    dplyr::rename(Chr=V1,length=V2) %>% mutate(genome = cumsum(as.numeric(length))) # dont shift by 1 here
chr.lengths$Chr <- as.numeric(gsub("X",23,gsub('^.{3}','',chr.lengths$Chr)))

# vector that has always the middle value for each chromosome (for plotting)
chr.lengths <- chr.lengths %>% mutate(chrbreaks = genome - length/2)

#View(cn.data)
###############################################################################
# plot
###############################################################################

pdf(file = paste(output.path,cohort,"_GLsubtypeprofiles.pdf", sep =""), height = 21.0, width = 72.0)

plot <- ggplot() +  
    ggtitle("Genome-wide frequency of gain/loss CN alterations") +
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqgain.HER2E)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqgain.HER2E),]$freqgain.HER2E, 
        color = "HER2E"),size=4) + 
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqloss.HER2E)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqloss.HER2E),]$freqloss.HER2E, 
        color = "HER2E"),size=4) + 
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqgain.LUMA)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqgain.LUMA),]$freqgain.LUMA, 
        color = "LUMA"),size=4) + 
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqloss.LUMA)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqloss.LUMA),]$freqloss.LUMA, 
        color = "LUMA"),size=4) + 
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqgain.LUMB)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqgain.LUMB),]$freqgain.LUMB, 
        color = "LUMB"),size=4) + 
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqloss.LUMB)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqloss.LUMB),]$freqloss.LUMB, 
        color = "LUMB"),size=4) + 
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqgain.HER2p)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqgain.HER2p),]$freqgain.HER2p, 
        color = "HER2p"),size=4) + 
    geom_line(aes(
        x = cn.data[which(!is.na(cn.data$freqloss.HER2p)),]$Genome_pos, 
        y = cn.data[!is.na(cn.data$freqloss.HER2p),]$freqloss.HER2p, 
        color = "HER2p"),size=4) + 
    scale_colour_manual(name="Subtype", values = c("HER2E"="#d334eb", "LUMA"="#2176d5", "LUMB"="#34c6eb", "HER2p"="#d72e2b")) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=2) +
    scale_x_continuous(name="Genome position (chromosome)",
                       breaks=chr.lengths$chrbreaks, 
                       labels=c(as.character(1:22),"X"),
                       limits = c(0,max(chr.lengths$genome)), #+50000000
                       expand = c(0, 0)) +
    scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100),
                       expand = c(0, 0),
                       limits = c(-100,100)) +
    theme(text=element_text(size=30),
          legend.title = element_blank(),
          axis.title.y = element_text(vjust = 0.5),
          legend.position = c(0.97, 0.95)) + #legend.position = "none") +
    annotate(x=min(cn.data$Genome_pos)+30000000,
             y=c(-50,50), label=c("Loss","Gain"), 
             geom="text", angle=90, hjust=0.5, 
             size=9, colour=c("black","black")) 
print(plot)
dev.off()