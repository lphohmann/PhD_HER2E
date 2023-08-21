# Script: testing which cn probes are significant

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" #  

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
# how to handle CN status for multiple probes matching to one gene: MAJORITY SIMPLY?
# do column wise for the samples and then just drop the ProbeIDs
for (gene in multiprobe.genes) {
  for (sample in colnames(cn.scanb)[grepl("^S",colnames(cn.scanb))]) {
    # see that is the majority CNA and take that
    
    # create extra df and then later rbind it with all the genes 
    # that are only have 1 mapped probe
  }
}
#View(head(cn.scanb))


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
