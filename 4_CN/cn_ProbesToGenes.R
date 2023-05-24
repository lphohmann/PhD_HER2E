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
library(readxl)
library(GenomicRanges)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

################################################################################
# loading data
################################################################################

gainloss.cn.scanb <- loadRData('data/SCANB/4_CN/raw/CN_gainloss.RData')

################################################################################
# get gene positions
# Question: to map genes do i need chr locations or genome locations
################################################################################

# chr seq need to be in format chr1 etc.
gainloss.cn.scanb$Chr <- sub("^", "chr", gainloss.cn.scanb$Chr)

# create granges object of probe positions
probes <- GRanges(gainloss.cn.scanb$Chr,IRanges(gainloss.cn.scanb$Position))
#target_range <- GRanges(target$chromosome, IRanges(start=target$start, end=target$end))

# genes of build hg38
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# intersect your regions of interest with the annotation
# 
overlaps <- findOverlaps(query = probes, 
                         subject = genes,
                         type="within", 
                         select="all")

res.df <- data.frame(query_hits = queryHits(overlaps),
                     subject_hits = subjectHits(overlaps)) 

# multiple probes map to one gene

head(res.df)
View(res.df)

f1 <- factor(subjectHits(overlaps), levels=seq_len(subjectLength(overlaps)))
splitAsList(mcols(x)[[column]][queryHits(overlaps)], f1)

#findOverlaps() returns a 'Hits' object that has two parallel vectors. 
# The vectors can be extracted with queryHits() and subjectHits(). 
# queryHits() are the indexes of the queries (the probe coordinates) that overlap the corresponding subjectHits(), i.e., the indexes of the subjects, the genes.


# alternative try
#Find which genes overlap which copy number regions
#splitByOverlap is that we can find which gene coordinates overlap which copy number variant coordinates, and then split the column of gene identifiers into lists corresponding to the regions of overlap
#query' is the gene coordinates, 'subject' the copy number coordinates. 'column' needs to match 'column' in the geneRanges() function.
geneRanges <- function(db, column="ENTREZID") {
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
}

splitColumnByOverlap <- function(query, subject, column="ENTREZID", ...) {
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }

symInCnv = splitColumnByOverlap(probes, genes, "gene_id",type="within", 
                                select="all")
symInCnv # CharacterList (list of character vectors) where each element contains the genes overlapping the corresponding CNV region.
genes[genes$gene_id %in% symInCnv[[1]]]

# or swapped
symInCnv = splitColumnByOverlap(genes,probes, "gene_id",type="within", 
                                select="all")
genes[genes$gene_id %in% symInCnv[[1]]]








#############################


#This returns an object telling which region overlaps with which annotation. The last step is to extract the gene names, collate them, and prepare your final output.

overlaps <- findOverlaps(query = probes, 
                         subject = genes,
                         type="within", 
                         select="all")

genes <- extractList(genes$gene_id, as(overlaps, "List"))
genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
res <- paste(as.character(probes), genes)
