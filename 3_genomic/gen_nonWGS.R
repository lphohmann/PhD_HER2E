# Script: Mutational analyses in Metabric (not WGS-based)

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "METABRIC"

# set/create output directory for plots
output.path <- "output/plots/3_genomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_driverWF.pdf",sep = "")
txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_HER2n_driverBarplots.txt", sep="")

#packages
source("scripts/3_genomic/src/gen_functions.R")
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
#library(Hmisc)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)
library(ggstatsplot)
library(GenVisR)

#######################################################################
# functions
#######################################################################

# function to convert the freq table to long format
countsToCases <- function(x, countcol = "Freq") {
    # Get the row indices to pull from x
    idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
    # Drop count column
    x[[countcol]] <- NULL
    # Get the rows from x
    x[idx, ]
}

# 1 mut per sample
count.sample <- function(data,gene) {
  data <- data %>% 
    distinct(Sample,Gene, .keep_all = TRUE)
  data[data$Gene==gene,] %>% 
    dplyr::count(PAM50) %>% 
    mutate(Freq=case_when(PAM50=="Her2e" ~ (n/pam50.n["Her2e"])*100,
                          PAM50=="LumA" ~ (n/pam50.n["LumA"])*100,
                          PAM50=="LumB" ~ (n/pam50.n["LumB"])*100))
}

#######################################################################
# 2. load data
#######################################################################

# load data
mut.data <- as.data.frame(read.delim('data/METABRIC/3_genomic/raw/data_mutations_extended.txt', header = FALSE, sep = "\t", dec = "."))

mut.data <- mut.data[-1,] %>% 
  row_to_names(row_number = 1) %>% 
  dplyr::rename(gene=Hugo_Symbol,variant_class=Variant_Classification,
                sample=Tumor_Sample_Barcode) %>% 
  dplyr::select(sample,gene,variant_class)

# load annotation data
anno <- loadRData("data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>% 
  filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
  filter(grepl('ERpHER2n', ClinGroup)) %>% 
  dplyr::rename(sample=METABRIC_ID) %>% 
  dplyr::select(sample,PAM50)
    
mut.data <- mut.data %>% filter(sample %in% anno$sample)

#######################################################################
# 3. Waterfall plots
#######################################################################

# Create a vector to save mutation priority order for plotting
mutation.priority <- as.character(unique(mut.data$variant_class))

mutation.priority <- c("Nonsense_Mutation","Missense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins","Splice_Region","Intron","3'UTR","Nonstop_Mutation","Translation_Start_Site","3'Flank","Silent")

# her2 wf plot
Her2.samples <- anno %>% filter(PAM50=="Her2") %>% pull(sample)
mut.data.her2 <- mut.data %>% filter(sample %in% Her2.samples)


wf.plot <- waterfall(mut.data.her2, 
                    fileType = "Custom", 
                    variant_class_order = mutation.priority,
                    mainGrid = TRUE,
                    #mainPalette = custom.pallete, # adapt this to match SCANB plot
                    main_geneLabSize = 15,
                    mainRecurCutoff = 0,
                    maxGenes = 30,
                    mainDropMut = TRUE, # drop unused mutation types from legend
                    #rmvSilent = TRUE,
                    out= "grob",
                    mutBurdenLayer = layer,
                    plotMutBurden = FALSE) #

# append to list
plot.list <- append(plot.list,list(wf.plot))

#######################################################################
# plot genes of interest: mutation frequencies
#######################################################################

# sample pam50, gene
head(mut.data)
head(anno)

# annotate samples with pam50
mut.data$PAM50 <- anno$PAM50[match(mut.data$sample,anno$sample)]

# adapt to scanb colnames (to use same functions)
mut.data <- mut.data %>% 
  dplyr::rename(Sample=sample,Gene=gene) %>% 
  mutate(PAM50 = ifelse(PAM50=="Her2","Her2e",PAM50))

# sample numbers
pam50.n <- table(
  mut.data[!duplicated(mut.data[,c("Sample")]),]$PAM50) 

# selected genes
gene.vec <-c("ERBB2","ESR1","TP53","PIK3CA","FGFR4")

for (g in gene.vec) {
  
  # Plots freq samples mutated in PAM50 groups
  p2 <- ggplot(count.sample(mut.data, gene=g), aes(fill=as.factor(PAM50),x=PAM50,y=Freq))+
    geom_bar(position="stack", stat="identity") +
    ggtitle(g) +
    theme_bw() +
    theme(aspect.ratio=1/1,
          legend.position = "none",
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black",linewidth = 2),
          axis.ticks = element_line(colour = "black", linewidth = 1),
          axis.ticks.length=unit(0.2, "cm")) +
    scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                      c("Her2e","LumA","LumB"))) +
    scale_y_continuous(breaks = scales::breaks_pretty(10),
                       limits = c(0,80)) +
    ylab("Mutation frequency (%)") +
    xlab("PAM50 subtype")
  
  plot.list <- append(plot.list,list(p2)) 
  
}

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE, height = 10, width = 15)

for (i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
  
  #print(plot.list[[i]])
}

dev.off()

#######################################################################
# stat testing mut. freqs
#######################################################################

# need df in binary format with samples as columns and genes as rows
sample.anno <- mut.data[c("Sample","PAM50")] %>% distinct(Sample, .keep_all = TRUE)
all.dmut.binary <- mut.data %>% 
  dplyr::select(-c(variant_class, PAM50)) %>% 
  distinct(Sample,Gene, .keep_all = TRUE) %>% 
  mutate(Mutation = 1) %>%  
  pivot_wider(names_from = Sample, values_from = Mutation) %>% 
  replace(is.na(.), 0) %>% 
  filter(Gene %in% gene.vec)

res.df <- data.frame()

pb = txtProgressBar(min = 0, max = nrow(all.dmut.binary), initial = 0, style = 3)
for (i in 1:nrow(all.dmut.binary)) { #nrow(gex.data)
  setTxtProgressBar(pb,i)
  
  # gene to test
  gene <- all.dmut.binary$Gene[i]
  
  # mutation counts
  her2e.n.mut <- sum(all.dmut.binary[all.dmut.binary$Gene==gene,
                                     sample.anno[sample.anno$PAM50=="Her2e","Sample"]]==1)
  luma.n.mut <- sum(all.dmut.binary[all.dmut.binary$Gene==gene,
                                    sample.anno[sample.anno$PAM50=="LumA","Sample"]]==1)
  lumb.n.mut <- sum(all.dmut.binary[all.dmut.binary$Gene==gene,
                                    sample.anno[sample.anno$PAM50=="LumB","Sample"]]==1)
  
  # make tbl
  freq.tbl <- data.frame(her2e = c(her2e.n.mut,
                                   length(sample.anno[sample.anno$PAM50=="Her2e","Sample"])-
                                     her2e.n.mut), 
                         luma = c(luma.n.mut,
                                  length(sample.anno[sample.anno$PAM50=="LumA","Sample"])-
                                    luma.n.mut), 
                         lumb = c(lumb.n.mut,
                                  length(sample.anno[sample.anno$PAM50=="LumB","Sample"])-
                                    lumb.n.mut),
                         row.names = c("mutation", "no_mutation"))
  
  # Her2 vs LumA
  luma.freq.tbl <- freq.tbl[,c("her2e","luma")]
  luma.res <- fisher.test(luma.freq.tbl)
  
  # for Her2 vs LumB
  lumb.freq.tbl <- freq.tbl[,c("her2e","lumb")]
  lumb.res <- fisher.test(lumb.freq.tbl)
  
  # save results gene | luma_pval | lumb_pval
  res.df <- rbind(res.df,c(gene,luma.res$p.value,lumb.res$p.value))
  txt.out <- append(txt.out,c(gene,capture.output(luma.res)))
  txt.out <- append(txt.out,c(gene,capture.output(lumb.res)))
  close(pb)
}

# name output columns
names(res.df) <- c("gene","luma.pval","lumb.pval")

# save text output
writeLines(txt.out, txt.file)