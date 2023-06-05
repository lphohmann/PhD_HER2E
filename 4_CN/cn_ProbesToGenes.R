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
                           probeID = gainloss.cn.scanb$ProbeID)
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) # meta = entrez ID

# convert to hgnc symbols
ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
mcols(genes)

################################################################################

# which gene coordinates overlap which copy number variant coordinates
overlap.res <- findOverlaps(probes, genes)
# queryHits(): indexes of the probe coordinates that overlap the corresponding 
# subjectHits(): indexes of the genes
# line up the query column identifier (probeID) that overlaps each gene
f1 <- factor(subjectHits(overlap.res), levels=seq_len(subjectLength(overlap.res)))
# use of factor() with exactly as many levels as there are subjects ensures that the splitAsList() command returns a 1:1 mapping between the subjects (genes) and the probes in the corresponding CharacterList
overlap.list <- splitAsList(mcols(probes)[["probeID"]][queryHits(overlap.res)], f1) # split the column of probe IDs into lists corresponding to the regions of overlap
mcols(genes) <- overlap.list
head(genes)
