# Script: 
# test gene level CN matrices for SCANB and BASIS for singificant differences in CN gain/loss alterations

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

#-------------------------------------------------------------------------------
# packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape2)
#library(data.table)
#library(purrr)
#library(readxl)

#-------------------------------------------------------------------------------
# set/create output directories
# for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)
# for data
data.path <- "data/SCANB/4_CN/processed/"
dir.create(data.path)
#-------------------------------------------------------------------------------
# input & output file paths
# input
scanb.gene.cna <- "data/SCANB/4_CN/processed/CNA_genelevel.RData"
basis.gene.cna <- "data/BASIS/4_CN/processed/CNA_genelevel_all.RData"
basis.anno <- "data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
# output
plot.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.pdf"
txt.file <- "output/plots/4_CN/COMBINED_HER2n_signifgenes.txt"
signif.genes <- "data/COMBINED/4_CN/processed/CNA_genelevel.RData"
#-------------------------------------------------------------------------------
# storing objects 
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
txt.out <- c() # object to store text output
# snippet to save plots (paste to end of script)
# pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)
# for(i in 1:length(plot.list)) { 
#   print(i)
#   print(plot.list[[i]])
# }
# dev.off()
# snippet to save text output (paste to end of script)
#writeLines(txt.out, txt.file)
#-------------------------------------------------------------------------------
# snipped to take exection time
#start.time <- Sys.time()
# CODE
#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken

################################################################################
# required data
################################################################################
scanb.cna <- loadRData(scanb.gene.cna)[["gainloss"]]
basis.cna <- loadRData(basis.gene.cna) %>% 
  dplyr::select(-c(chr,centerPos,Genome_pos))
View(basis.cna)
basis.anno <- loadRData(basis.anno) %>% 
  filter(final.ER=="positive" & final.HER2=="negative") %>%
  dplyr::filter(PAM50_AIMS %in% c("LumA","LumB")) %>%
  dplyr::rename(PAM50=PAM50_AIMS) %>% 
  dplyr::select(sample_name,PAM50)

luma.ids <- basis.anno %>% filter(PAM50 == "LumA") %>% pull(sample_name)
lumb.ids <- basis.anno %>% filter(PAM50 == "LumB") %>% pull(sample_name)

#setdiff(luma.ids,names(basis.cna)) #1 luma sample doesnt have wgs data

# include only common genes
common.genes <- intersect(scanb.cna$gene,basis.cna$gene)
scanb.cna <- scanb.cna %>% filter(gene %in% common.genes)
basis.cna <- basis.cna %>% filter(gene %in% common.genes)

# add gene annotation to basis from scanb
basis.cna$centerPos <- scanb.cna$centerPos[match(basis.cna$gene, scanb.cna$gene)]
basis.cna$Genome_pos <- scanb.cna$Genome_pos[match(basis.cna$gene, scanb.cna$gene)]
basis.cna$chr <- scanb.cna$chr[match(basis.cna$gene, scanb.cna$gene)]
basis.cna <- basis.cna %>% relocate(c(chr,centerPos,Genome_pos), .after=gene)

# save obj
# gene LumA.Gain.pval LumB.Gain.pval HER2p.Gain.pval LumA.Loss.pval LumB.Loss.pval HER2p.Loss.pval
i=1

# which test: chi2?
for(i in 1:length(common.genes)) {
  
  # get data
  gene <- common.genes[i]
  basis.gene.dat <- basis.cna[basis.cna$gene==gene,][,!names(basis.cna) %in% c(
    "gene","chr","centerPos","Genome_pos")]
  # 
  her2e.dat <- as.numeric(scanb.cna[scanb.cna$gene==gene,][,!names(scanb.cna) %in% c("gene","chr","centerPos","Genome_pos")])
  luma.dat <- as.numeric(basis.gene.dat[,names(basis.gene.dat) %in% luma.ids]) 
  lumb.dat <- as.numeric(basis.gene.dat[,names(basis.gene.dat) %in% lumb.ids]) 
  
  comp.list <- list("LumA"=luma.dat,"LumB"=lumb.dat)
  # test for all pairs
  j = 1
  for(j in 1:length(comp.list)) {
    comp.dat <- comp.list[[j]]
    # get cross tables
    gain.tbl <- data.frame(her2e=c(sum(her2e.dat>0),sum(her2e.dat<1)), 
                           comp=c(sum(comp.dat>0),sum(comp.dat<1)), 
                           row.names = c("gain","no_gain"))
    loss.tbl <- data.frame(her2e=c(sum(her2e.dat<0),sum(her2e.dat>-1)), 
                           comp=c(sum(comp.dat<0),sum(comp.dat>-1)), 
                           row.names = c("loss","no_loss"))
    # test gain
    if(min(chisq.test(gain.tbl)$expected)<5) {
      gain.pval <- fisher.test(gain.tbl)$p.value
    } else {
      gain.pval <- chisq.test(gain.tbl)$p.value 
    }
    # test loss
    if(min(chisq.test(loss.tbl)$expected)<5) {
      loss.pval <- fisher.test(loss.tbl)$p.value
    } else {
      loss.pval <- chisq.test(loss.tbl)$p.value 
    }
    # save in vectors
    
    
  }
}

# adjust pval

# gene LumA.Gain.pval LumB.Gain.pval HER2p.Gain.pval LumA.Loss.pval LumB.Loss.pval HER2p.Loss.pval


#x <- loadRData("data/BASIS/4_CN/processed/ERpHER2p_CollectedFrequencyData.RData")

l