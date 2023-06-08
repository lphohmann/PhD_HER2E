# Script: create key file for mapping probes to genes
# Purpose: 
# 1. map the probes to genes

# used before cn_HER2Eprocessing.R
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

probe.positions <- loadRData('data/SCANB/4_CN/raw/CN_gainloss.RData')[c("ProbeID","Chr","Position")]

#head(probe.positions)

################################################################################
# Make key file for which probes map to which genes
################################################################################

# chr seq need to be in format chr1 etc.
probe.positions$Chr <- sub("^", "chr", probe.positions$Chr)

# create granges objects
probes <- GRanges(seqnames = probe.positions$Chr,
                           ranges = IRanges(probe.positions$Position),
                           ProbeID = probe.positions$ProbeID)
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
#head(key.df) # why multiple genes per probe??

# add to matrix
probe.positions <- merge(probe.positions,key.df,by=c("Chr","ProbeID")) %>% 
  relocate(Gene_symbol,.after = ProbeID)

#head(probe.positions)

################ test: why multiple genes for some probe - wrong gene loc annotation?
# t <- merge(as.data.frame(probes),probe.anno,
#          by=c("seqnames","start","end","width","strand"))
# t[1,]$X
# probes.df <- as.data.frame(probes)
# probes.df[which(probes.df$seqnames=="chr1" & probes.df$start == 100000723),]
# gene.anno %>% filter(SYMBOL %in% c("SLC35A3","MFSD14A"))

###############################################################################

# replace character(0) with NA for later filtering of probes without annotation
probe.positions$Gene_symbol <- lapply(probe.positions$Gene_symbol, 
                                      function(x) {
                                        if (identical(x, character(0))) {
                                          return(NA)
                                          } else {return(x)}
                                        })

# filter NA rows
probe.positions <- probe.positions[which(!is.na(probe.positions$Gene_symbol)),]

# which probes have multiple mapped genes
#multi.genes <- probe.positions[which(lengths(probe.positions$Gene_symbol)>1),]
#nrow(multi.genes) #107629

####################

save(probe.positions, 
     file = paste(data.path,"CN_mapped_probes.RData",sep=""))


