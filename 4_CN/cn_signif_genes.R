# Script: testing which genes are significantly different in terms of CNA freq between subtypes

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" 

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for data
data.path <- "data/SCANB/4_CN/processed/"
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_signifprobes.pdf",sep = "")

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
#library(data.table)
library(purrr)
library(readxl)
library(IRanges)

################################################################################
# SCANB get CN data on gene level (1 row per gene)
################################################################################

# get gene positions & convert to hgnc symbols
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) # meta = entrez ID)
ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
# convert to df and filter out NA gene symbols
genes <- annoGR2DF(genes) %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(chr = gsub("^chr","",chr)) %>% 
  filter(chr %in% c(1:22)) %>% 
  dplyr::select(-c(gene_id,strand, width)) %>% 
  mutate_at(c('chr','start','end'), as.numeric)

# get segment data
cn.scanb.segments <- loadRData("data/SCANB/4_CN/processed/Segment_CN_states.RData")

# 
head(cn.scanb.segments[[1]])
head(genes)

# result storing objects
gl.names <- c("gene", "chr", "start", "end",
              c(names(cn.scanb.segments)))
amp.names <- c("gene", "chr", "start", "end",
               c(names(cn.scanb.segments)))
gl.df <- as.data.frame(matrix(nrow=length(genes$SYMBOL),ncol=length(gl.names)))
names(gl.df) <- gl.names
amp.df <- as.data.frame(matrix(nrow=length(genes$SYMBOL),ncol=length(amp.names)))
names(amp.df) <- amp.names

# loop over genes
pb = txtProgressBar(min = 0, max = length(genes$SYMBOL), initial = 0, style = 3)
for (i in 1:nrow(genes)) {
  setTxtProgressBar(pb,i)
  
  gene.dat <- genes[i,]
  #print(gene.dat)
  # gene.chr <- genes$chr[i]
  # gene.start <- genes$start[i]
  # gene.end <- genes$end[i]

  # get the gainloss and amp state of that gene for all samples
  gene.statuses <- sapply(cn.scanb.segments, function(segment.df) {
    segment.df <- segment.df %>% 
      mutate_at(c('chr','startpos','endpos','GainLoss','Amp'), as.numeric)
    # get the cn status of the segment the gene is on
    # only relevant segments
    chr.segments <- segment.df[segment.df$chr==gene.dat$chr,]
    
    segments.query <- with(chr.segments, IRanges(startpos, endpos))
    gene.subject <- with(gene.dat, IRanges(start, end))
    # check 
    chr.segments$overlap = countOverlaps(segments.query, gene.subject) != 0
    if (sum(chr.segments$overlap)==0) { # no segment covers the gene region
      gene.gl.state <- NA
      gene.amp.state <- NA
    } else {
      hits <- findOverlaps(segments.query, gene.subject)
      overlaps <- pintersect(segments.query[queryHits(hits)], gene.subject[subjectHits(hits)])
      percentOverlap <- width(overlaps) / width(gene.subject[subjectHits(hits)])
      chr.segments <- chr.segments[chr.segments$overlap == TRUE,]
      chr.segments$percentOverlap <- percentOverlap
      # only select the segment with the highest overlap proportion
      main.segment <- chr.segments[which.max(chr.segments$percentOverlap),]
      gene.gl.state <- main.segment$GainLoss
      gene.amp.state <- main.segment$Amp
      }
    
    return(list("GainLoss"=gene.gl.state,"Amp"=gene.amp.state))
    })

  # store results
  #print(gene.statuses) #527 #KCNJ18 #S003268 #numeric,0
  
  gl.df[i,] <- c(gene.dat$SYMBOL,
                 gene.dat$chr,
                 gene.dat$start,
                 gene.dat$end,
                 unname(unlist(gene.statuses["GainLoss",])))
  
  amp.df[i,] <- c(gene.dat$SYMBOL,
                  gene.dat$chr,
                  gene.dat$start,
                  gene.dat$end,
                  unname(unlist(gene.statuses["Amp",])))

  close(pb)
}

cn.list <- list("gainloss"=gl.df,"amp"=amp.df)

save(cn.list,
     file = "data/SCANB/4_CN/processed/CNA_genelevel.RData")



# add center position
# convert to genome position


################################################################################
# BASIS get CN data on gene level (1 row per gene)
################################################################################

# hg19 genes
# get gene positions & convert to hgnc symbols
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
# convert to df and filter out NA gene symbols
genes <- annoGR2DF(genes) %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(chr = gsub("^chr","",chr)) %>% 
  filter(chr %in% c(1:22)) %>% 
  dplyr::select(-c(gene_id,strand, width)) %>% 
  mutate_at(c('chr','start','end'), as.numeric)

# get segment data
cn.basis.segments <- loadRData("data/BASIS/4_CN/raw/ASCAT_CEL_Total/ASCAT_CEL_Total_EasySegments.RData")
#View(head(cn.basis.segments))
# organize in list of dfs (1 per sample)
cn.basis.segments <- with(cn.basis.segments, split(
  cn.basis.segments, list(SampleID)))
cn.basis.segments <- lapply(cn.basis.segments, function(x) {
  x <- x %>% dplyr::rename(GainLoss=CNA,startpos=start,endpos=stop)
  return(x)
})

# 
head(cn.basis.segments[[1]])
head(cn.scanb.segments[[1]])

head(genes)
# do same as above



# add center position
# convert to genome position

# liftover genome position to hg38



### 

# get gene positions & convert to hgnc symbols
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) # meta = entrez ID)
ENTREZID2SYMBOL <- select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
# convert to df and filter out NA gene symbols
genes <- annoGR2DF(genes) %>% 
  filter(!is.na(SYMBOL)) %>% 
  mutate(chr = gsub("^chr","",chr)) %>% 
  filter(chr %in% c(1:22)) %>% 
  dplyr::select(-c(gene_id,strand, width)) %>% 
  mutate_at(c('chr','start','end'), as.numeric)

# get segment data
cn.scanb.segments <- loadRData("data/SCANB/4_CN/processed/Segment_CN_states.RData")

# add genome positions (depends what position I have for the genes)
head(cn.scanb.segments[[1]])
head(genes)

# result storing objects
gl.names <- c("gene", "chr", "start", "end",
              c(names(cn.scanb.segments)))
amp.names <- c("gene", "chr", "start", "end",
               c(names(cn.scanb.segments)))
gl.df <- as.data.frame(matrix(nrow=length(genes$SYMBOL),ncol=length(gl.names)))
names(gl.df) <- gl.names
amp.df <- as.data.frame(matrix(nrow=length(genes$SYMBOL),ncol=length(amp.names)))
names(amp.df) <- amp.names

# loop over genes
pb = txtProgressBar(min = 0, max = length(genes$SYMBOL), initial = 0, style = 3)
for (i in 1:nrow(genes)) {
  setTxtProgressBar(pb,i)
  
  gene.dat <- genes[i,]
  
  # get the gainloss and amp state of that gene for all samples
  gene.statuses <- sapply(cn.scanb.segments, function(segment.df) {
    segment.df <- segment.df %>% 
      mutate_at(c('chr','startpos','endpos','GainLoss','Amp'), as.numeric)
    # get the cn status of the segment the gene is on
    # only relevant segments
    chr.segments <- segment.df[segment.df$chr==gene.dat$chr,]
    
    segments.query <- with(chr.segments, IRanges(startpos, endpos))
    gene.subject <- with(gene.dat, IRanges(start, end))
    # check 
    chr.segments$overlap = countOverlaps(segments.query, gene.subject) != 0
    if (sum(chr.segments$overlap)==0) { # no segment covers the gene region
      gene.gl.state <- NA
      gene.amp.state <- NA
    } else {
      hits <- findOverlaps(segments.query, gene.subject)
      overlaps <- pintersect(segments.query[queryHits(hits)], gene.subject[subjectHits(hits)])
      percentOverlap <- width(overlaps) / width(gene.subject[subjectHits(hits)])
      chr.segments <- chr.segments[chr.segments$overlap == TRUE,]
      chr.segments$percentOverlap <- percentOverlap
      # only select the segment with the highest overlap proportion
      main.segment <- chr.segments[which.max(chr.segments$percentOverlap),]
      gene.gl.state <- main.segment$GainLoss
      gene.amp.state <- main.segment$Amp 
    }
    
    return(list("GainLoss"=gene.gl.state,"Amp"=gene.amp.state))
  })
  
  # store results
  #print(gene.statuses) #527 #KCNJ18 #S003268 #numeric,0
  
  gl.df[i,] <- c(gene.dat$SYMBOL,
                 gene.dat$chr,
                 gene.dat$start,
                 gene.dat$end,
                 unname(unlist(gene.statuses["GainLoss",])))
  
  amp.df[i,] <- c(gene.dat$SYMBOL,
                  gene.dat$chr,
                  gene.dat$start,
                  gene.dat$end,
                  unname(unlist(gene.statuses["Amp",])))
  
  close(pb)
}


save(list(gl.df,amp.df), 
     file = "data/SCANB/4_CN/processed/CNA_genelevel.RData")

View(gl.df)


# add center position
# convert to genome position




################################################################################
# load scanb CNA data for all probes
cn.scanb <- loadRData("data/SCANB/4_CN/processed/CN_gainloss_genpos_genmap.RData")
# bring it down to gene level
cn.scanb <- unnest(cn.scanb,cols=c(Gene_symbol)) # duplicates rows where probe matched to two genes so its a 1:1 mapping
cn.scanb <- cn.scanb[!is.na(cn.scanb$Gene_symbol),] 
# genes with more than 1 probe mapping to them
n_occur <- data.frame(table(cn.scanb$Gene_symbol))
multiprobe.genes <- cn.scanb$Gene_symbol[
  cn.scanb$Gene_symbol %in% n_occur$Var1[n_occur$Freq > 1]]

# get cn data for genes with multiple probes
splitgenes <- c()
res.df <- as.data.frame(matrix(nrow=length(multiprobe.genes),ncol=length(colnames(cn.scanb))))
names(res.df) <- colnames(cn.scanb)
pb = txtProgressBar(min = 0, max = length(multiprobe.genes), initial = 0, style = 3)
for (i in 1:length(multiprobe.genes)) {
  setTxtProgressBar(pb,i)
  gene <- multiprobe.genes[i]
  gene.dat <- cn.scanb[which(cn.scanb$Gene_symbol==gene),]
  # Chr, ProbeID, Gene_symbol, Position, Genome_pos
  # genome positions are defined as the mean of the probe positions
  res.df$Chr[i] <- gene.dat$Chr[1]
  #res.df$ProbeID[i] <- NA
  res.df$Gene_symbol[i] <- gene
  res.df$Position[i] <- mean(gene.dat$Position)
  res.df$Genome_pos[i] <- mean(gene.dat$Genome_pos)
  
  for (j in 1:length(colnames(cn.scanb)[grepl("^S",colnames(cn.scanb))])) {
    sample <- colnames(cn.scanb)[grepl("^S",colnames(cn.scanb))][j]
    # see that is the majority CNA and take that
    sample.dat <- gene.dat[sample]
    tbl <- table(sample.dat[[sample]])
    most.common <- sort(names(tbl[tbl==max(tbl)]),decreasing = TRUE)
    if (length(most.common) > 1) {
      cn.state <- most.common[1] # takes the highest value then
      splitgenes <- append(splitgenes, gene) # save the genes that have 2 even cn states
    } else { cn.state <- most.common }
    
    res.df[[sample]][i] <- cn.state

  }
  close(pb)
}

save(res.df, 
     file = "data/SCANB/4_CN/processed/CNA_genelevel_prelim.RData")

# combine with the rest of the data
cn.scanb.final <- rbind(cn.scanb[which(!(cn.scanb$Gene_symbol %in% multiprobe.genes)),],res.df)
save(cn.scanb.final, 
     file = "data/SCANB/4_CN/processed/CNA_genelevel.RData")
a <- loadRData("data/SCANB/4_CN/processed/CNA_genelevel.RData")
View(cn.scanb.final)

# these are the probes mapped to the genes for basis
cn.basis <- loadRData(
  "data/BASIS/4_CN/processed/CN_gainloss_frequencies_genpos_genmap.RData")[1:5]
# ACHTUNG: these are post liftover, so yeah dont use them to get directly the segment cna data

# map the CN data to probes
# check that i use the pre-liftover positions

# 1. map pre-liftover probes to pre-liftover segments to annotate with CNA on probe level
pre.probes <- loadRData("data/BASIS/4_CN/processed/LumA_CollectedFrequencyData.RData")[["fData"]]
pre.segments <- loadRData(
  "./data/BASIS/4_CN/raw/ASCAT_CEL_Total/ASCAT_CEL_Total_EasySegments.RData") %>% 
  dplyr::select(c("SampleID","chr","start","stop","CNA"))
#View(head(pre.segments))
#View(head(pre.probes))

# need to do this with apply (runs over all rows or columns) or lapply (runs over all elements, returns list)

# df to store res
sample.cna.list <- list()

pb = txtProgressBar(min = 0, max = length(unique(pre.segments$SampleID)), initial = 0, style = 3)
for (i in 1:length(unique(pre.segments$SampleID))) {
  setTxtProgressBar(pb,i)
  
  sample <- unique(pre.segments$SampleID)[i]
  pre.segments.sample <- pre.segments[pre.segments$SampleID==sample,]
  
  # store res
  res <- c()
  
  for (probe in pre.probes$reporterId[1:50]) { #pre.probes$reporterId) {
    
    # probe data
    probe.pos <- pre.probes[pre.probes$reporterId==probe,]$centerPosition
    probe.chr <- pre.probes[pre.probes$reporterId==probe,]$chromosome
    # get segment cn state
    res <- append(res, 
                  pre.segments.sample[which(
                    pre.segments.sample$chr == probe.chr &
                    pre.segments.sample$start <= probe.pos &
                    pre.segments.sample$stop >= probe.pos),]$CNA)
  }
  
  # store result in list
  sample.cna.list <- append(sample.cna.list, as.data.frame(res))
  
  close(pb)
}

sample.cna.df <- as.data.frame(do.call(cbind, sample.cna.list))
names(sample.cna.df) <- unique(pre.segments$SampleID)

#View(sample.cna.df)

#pre.probes.cna <- cbind(pre.probes,sample.cna.df) # this should work when i run over all probes
pre.probes.cna <- cbind.fill(pre.probes,sample.cna.df)
pre.probes.cna <- as.data.frame(pre.probes.cna)
#View(pre.probes.cna)

save(pre.probes.cna, 
     file = "data/BASIS/4_CN/processed/probe_CNA.RData")
load("data/BASIS/4_CN/processed/probe_CNA.RData")
View(head(pre.probes.cna))

# 2. use liftover result to update the positions of the probes to hg38 like scanb
post.probes <- loadRData(file = "data/BASIS/4_CN/processed/CN_gainloss_frequencies.RData")[1:4]

# add position
pre.probes.cna$Position.pl <- post.probes$Position[match(
  pre.probes.cna$reporterId, post.probes$ProbeID)]
pre.probes.cna$Genome_pos.pl <- post.probes$Genome_pos[match(
  pre.probes.cna$reporterId, post.probes$ProbeID)]
pre.probes.cna$centerPosition <- NULL
post.probes.cna <- pre.probes.cna %>% 
  relocate(Position.pl, .after = chromosome) %>% 
  relocate(Genome_pos.pl, .after = Position.pl)

save(post.probes.cna, 
     file = "data/BASIS/4_CN/processed/probe_CNA.RData")
# 3. check how many genome positions match to do statistics
# BUT what to do with the rest that doesnt match up between scanb and basis?
# how many positional matches that allow statistics between scanb and basis?

View(head(cn.scanb))
View(head(post.probes.cna))

(length(intersect(cn.scanb$Genome_pos,post.probes.cna$Genome_pos.pl))/length(cn.scanb$Genome_pos))*100 # only 12% of scanb probes have an exact match in basis

# any already mapped rpobes?
s <- loadRData("data/BASIS/4_CN/processed/LumA_CollectedFrequencyData.RData")
str(s)               
View(s)               
