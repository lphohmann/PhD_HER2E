# Script: Process SCANB CN data to prepare for plotting mutation waterfall
# Purpose: 
# 1. map the probes to genes for the driver mutation plots 

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
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

################################################################################
# loading data
################################################################################

gainloss.cn.scanb <- loadRData('data/SCANB/4_CN/raw/CN_gainloss.RData')

################################################################################
# get gene positions
################################################################################

# chr seq need to be in format chr1 etc.
gainloss.cn.scanb$Chr <- sub("^", "chr", gainloss.cn.scanb$Chr)

# create granges objects
probes <- GRanges(seqnames = gainloss.cn.scanb$Chr,
                           ranges = IRanges(gainloss.cn.scanb$Position),
                           ProbeID = gainloss.cn.scanb$ProbeID)
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) # meta = entrez ID

# convert to hgnc symbols
ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
gene.anno <- as.data.frame(genes)
probe.anno <- as.data.frame(probes)

###############################################################################

# which 
overlap.res <- findOverlaps(genes,probes)
# queryHits(): indexes of the gene coordinates that overlap the corresponding 
# subjectHits(): indexes of the probes
# line up the query column identifier (gene) that overlaps each probe
f1 <- factor(subjectHits(overlap.res), levels=seq_len(subjectLength(overlap.res)))
# use of factor() with exactly as many levels as there are subjects ensures that the splitAsList() command returns a 1:1 mapping between the subjects (probes) and the genes in the corresponding CharacterList
overlap.list <- splitAsList(mcols(genes)[["SYMBOL"]][queryHits(overlap.res)], f1) # split the column of gene symbols into lists corresponding to the regions of overlap
mcols(probes) <- overlap.list

key.df <- merge(as.data.frame(probes),probe.anno,
                by=c("seqnames","start","end","width","strand")) %>% 
  dplyr::rename(Gene_symbol=X,Position=start,Chr=seqnames) %>% 
  dplyr::select(Chr,ProbeID,Gene_symbol)
head(key.df) # why multiple genes per probe??

#test
t <- merge(as.data.frame(probes),probe.anno,
         by=c("seqnames","start","end","width","strand"))
t[1,]$X
probes.df <- as.data.frame(probes)
probes.df[which(probes.df$seqnames=="chr1" & probes.df$start == 100000723),]
genes.df <- as.data.frame(genes)
head(genes.df)
genes.df[which(genes.df$seqnames=="chr1" & genes.df$start >= 100000723 & genes.df$end <= 100000723),]
genes.df[which(genes.df$SYMBOL =="SLC35A3" | genes.df$SYMBOL == "MFSD14A"),]
# MFSD14A has the wrong starting coordinates?????
