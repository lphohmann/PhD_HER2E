# Script: Copy number alteration analyses in the Metabric cohort

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "METABRIC"

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/",sep="")
dir.create(data.path)

# store plots
plot.list <- list()

#packages
source("scripts/4_CN/src/cn_functions.R")
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(car)
library(ggsignif)
library(biomaRt)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# 2.1 load data
cn.data <- read.delim("./data/METABRIC/4_CN/raw/data_CNA.txt",header = TRUE, sep = "\t", dec = ".") %>% 
  dplyr::select(-c(Entrez_Gene_Id)) %>% 
  column_to_rownames(var = "Hugo_Symbol")

# annotation data
anno <- loadRData("./data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>% 
  as.data.frame(.) %>% 
  filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
  filter(grepl('ERpHER2n', ClinGroup)) %>% 
  dplyr::rename(sampleID=METABRIC_ID) %>% 
  dplyr::select(sampleID,PAM50) %>% 
  mutate(sampleID = gsub("-",".",sampleID)) # have to get identifiers in same format

cn.data <- cn.data  %>% 
  dplyr::select(any_of(anno$sampleID))

anno <- anno %>% 
  filter(sampleID %in% colnames(cn.data))

### quick FGFR4 check - put in manuscript
#prop.table(table(t(cn.data["FGFR4",c(anno[which(anno$PAM50=="Her2"),]$sampleID)])))*100
#prop.table(table(t(cn.data["FGFR4",c(anno[which(anno$PAM50=="LumA"),]$sampleID)])))*100
#prop.table(table(t(cn.data["FGFR4",c(anno[which(anno$PAM50=="LumB"),]$sampleID)])))*100
  
##########################################################################
# Part 1: general gain/loss, whole genome, freq. altered per group boxplot
##########################################################################

# alterations per sample 
res <- colSums(cn.data != 0) %>% as.data.frame() %>% rownames_to_column(var = "sampleID") %>% dplyr::rename(Nalt = 2)

#as percent
res$Nalt <- (res$Nalt/nrow(cn.data))*100

# add to anno
p1.anno <- merge(anno,res,by="sampleID") %>% drop_na(Nalt)

# calc. significance
# luma
p1.anno$PAM50 <- as.factor(p1.anno$PAM50)
adata <- p1.anno %>% filter(PAM50 %in% c("LumA", "Her2"))
bdata <- p1.anno %>% filter(PAM50 %in% c("LumB", "Her2"))
# check equal variance
leveneTest(Nalt~PAM50, data=adata) # signif.
t.test(Nalt~PAM50, data=adata, var.equal = FALSE)$p.value # signif. 

# lumb
leveneTest(Nalt~PAM50, data=bdata) # not signif.
t.test(Nalt~PAM50, data=bdata, var.equal = TRUE)$p.value # 0.17-> not signif. 

# plot
plot.list <- append(plot.list, list(
  three_boxplot(p1.anno,
              group.var = "PAM50",
              test.var = "Nalt",
              g1="Her2",g2="LumA",g3="LumB",
              g3.pos = 100, g3.sign = "ns",
              g2.pos = 110, g2.sign = "****",
              ylim = c(0,120),
              ylab = "Genome altered (%)", 
              title = "Genome alterations in PAM50 subtypes (ERpHER2n)",
              break.step = 10)))

#######################################################################
# Part 2: spatial gain/loss - get chromosomal position and cytoband from gene symbols 
#######################################################################

# apprach 1: - Biomart & Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://grch37.ensembl.org") # GRCh37

# Only use standard human chromosomes
normal.chroms <- c(1:22, "X", "Y")

# Filter on HGNC symbol and chromosome, retrieve genomic location and band
my.symbols <- row.names(cn.data)

my.regions <- getBM(
  c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
  filters = c("hgnc_symbol", "chromosome_name"),
  values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
  mart = ensembl)

# get data into the right format
cn.data.meta <- as.data.frame(merge(my.regions,rownames_to_column(cn.data,var = "hgnc_symbol"), by="hgnc_symbol")) %>% 
  dplyr::rename(chr=chromosome_name,start=start_position,end=end_position) %>%
  mutate(pos=end-((end-start)/2)) %>% 
  relocate(pos, .after=band) %>% 
  select(-c("start","end"))

# rename X chromosome to 23
cn.data.meta$chr[cn.data.meta$chr == "X"] <- "23" 
cn.data.meta$chr <- as.numeric(cn.data.meta$chr)

# sort
cn.data.meta <- cn.data.meta[with(cn.data.meta, order(chr, pos)),] 

# get chromosome lengths (i simply take the last pos on each)
chr.lengths <- cn.data.meta %>% 
  group_by(chr) %>% 
  summarise(length = max(pos)) %>% 
  as.data.frame()
chr.lengths <- chr.lengths[with(chr.lengths, order(chr)),] %>% 
  mutate(genome = cumsum(length)) %>% #dplyr::rename(chromosome=chr)
  mutate(chrbreaks = genome - length/2)



# add genome pos 
cn.data.meta <- cn.data.meta %>% 
  rowwise() %>% 
  mutate(genome = ifelse(chr==1,pos,pos + sum(chr.lengths$length[1:chr-1]))) %>%
  relocate(genome, .after=pos)

# for each group calc. the frequency of gain and loss for each gene
# sampleIDs
HER2.samples <- anno[which(anno$PAM50 == "Her2"),]$sampleID
LUMA.samples <- anno[which(anno$PAM50 == "LumA"),]$sampleID
LUMB.samples <- anno[which(anno$PAM50 == "LumB"),]$sampleID

# calc loss/gain freqs per group
cn.data.meta$freqloss.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x<0))/length(HER2.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqloss.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x<0))/length(LUMA.samples))*-100)
cn.data.meta$freqloss.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x<0))/length(LUMB.samples))*-100)
cn.data.meta$freqgain.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x>0))/length(HER2.samples))*100)
cn.data.meta$freqgain.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x>0))/length(LUMA.samples))*100)
cn.data.meta$freqgain.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x>0))/length(LUMB.samples))*100)

# plot
# her2
plot.list <- append(plot.list, list(
  ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of gain/loss CN alterations in luminal HER2-enriched BC") +
    geom_line(aes(y = freqgain.H2), color = "darkgreen",linewidth=3) + 
    geom_line(aes(y = freqloss.H2), color = "darkred",linewidth=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",linewidth=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$chrbreaks,
                       labels=c(1:22,"X"), #chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100)) +
    theme(text=element_text(size=27))))

# combined
plot.list <- append(plot.list, list(
  ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of gain/loss CN alterations (ERpHER2n)") +
    geom_line(aes(y = freqgain.H2, color = "red"),linewidth=4) + 
    geom_line(aes(y = freqloss.H2, color = "red"),linewidth=4) + 
    geom_line(aes(y = freqgain.LUMA, color = "green"),linewidth=4) + 
    geom_line(aes(y = freqloss.LUMA, color = "green"),linewidth=4) + 
    geom_line(aes(y = freqgain.LUMB, color = "blue"),linewidth=4) + 
    geom_line(aes(y = freqloss.LUMB, color = "blue"),linewidth=4) + 
    scale_colour_manual(name="PAM50", values = c("#d334eb", "#2176d5", "#34c6eb"), labels = c("HER2", "LUMA", "LUMB")) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=2) +
    scale_x_continuous(name="Genome position (chromosome)",
                       breaks=chr.lengths$chrbreaks,
                       labels=c(as.character(1:22),"X"),
                       limits = c(0,max(cn.data.meta$genome)+50000000),
                       expand = c(0, 0)) +
    scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100),
                       expand = c(0, 0),
                       limits = c(-80,80)) +
    theme(text=element_text(size=30),
          legend.title = element_blank(),
          axis.title.y = element_text(vjust = 0.5),
          legend.position = c(0.97, 0.90)) +
    annotate(x=45000000,y=c(-50,50), label=c("Loss","Gain"), 
    geom="text", angle=90, hjust=0.5, size=9, colour=c("black","black"))))

##################################################################
# Part 3: spatial driver gain/loss - plot CN driver variation of each group
##################################################################

# calc loss/gain freqs per group
cn.data.meta$freqloss.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x<(-1)))/length(HER2.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqloss.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x<(-1)))/length(LUMA.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqloss.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x<(-1)))/length(LUMB.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqgain.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x>1))/length(HER2.samples))*100)
cn.data.meta$freqgain.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x>1))/length(LUMA.samples))*100)
cn.data.meta$freqgain.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x>1))/length(LUMB.samples))*100)

# plot
#her2
plot.list <- append(plot.list, list(ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of driver gain/loss CN alterations in luminal HER2-enriched BC") +
    geom_line(aes(y = freqgain.H2), color = "darkgreen",size=3) + 
    geom_line(aes(y = freqloss.H2), color = "darkred",size=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$chrbreaks,
                       labels=chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,5)),
                       labels=c(seq(-100,100,5))) +
    theme(text=element_text(size=27)))) #+

# combined
plot.list <- append(plot.list, list(ggplot(cn.data.meta, aes(x=genome)) + #ggtitle("Frequency of driver gain/loss CN alterations in luminal BC") +
    geom_line(aes(y = freqgain.H2, color = "red"),size=4) + 
    geom_line(aes(y = freqloss.H2, color = "red"),size=4) + 
    geom_line(aes(y = freqgain.LUMA, color = "green"),size=4) + 
    geom_line(aes(y = freqloss.LUMA, color = "green"),size=4) + 
    geom_line(aes(y = freqgain.LUMB, color = "blue"),size=4) + 
    geom_line(aes(y = freqloss.LUMB, color = "blue"),size=4) + 
    scale_colour_manual(name="PAM50", values = c("#d334eb", "#2176d5", "#34c6eb"),
                        labels = c("HER2", "LUMA", "LUMB")) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=2) +
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="Genome position (chromosome)",
                       breaks=chr.lengths$chrbreaks,
                       labels=c(as.character(1:22),"X"),
                       limits = c(0,max(cn.data.meta$genome)+50000000),
                       expand = c(0, 0)) +
    scale_y_continuous(name="Alteration frequency (%)", 
                       breaks=c(seq(-20,50,10)),
                       labels=c(20,10,0,10,20,30,40,50),
                       expand = c(0, 0),
                       limits = c(-15,45)) +
    theme(text=element_text(size=30),
          legend.title = element_blank(),
          axis.title.y = element_text(vjust = 0.5),
          legend.position = c(0.97, 0.90)) +
    annotate(x=45000000,y=c(-10,10), label=c("Loss","Gain"), 
             geom="text", angle=90, hjust=0.5, size=9, colour=c("black","black"))))

#plot
pdf(file = paste(output.path,cohort,"_HER2n_CNprofiles.pdf", sep=""), 
    onefile = TRUE, width = 29.7, height = 21) 

for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}

dev.off()

##################################################################
# calc. singificance of CN alterations between groups
##################################################################

# genes
genes <- cn.data.meta$hgnc_symbol

# objects to store results
res.matrix <- data.frame(matrix(
  ncol = 9, 
  nrow = length(cn.data.meta$hgnc_symbol))) %>% 
  dplyr::rename(Gene=1, Loss_LUMA_pval=2, Loss_LUMB_pval=3, 
                Gain_LUMA_pval=4, Gain_LUMB_pval=5,
                DLoss_LUMA_pval=6, DLoss_LUMB_pval=7, 
                DGain_LUMA_pval=8, DGain_LUMB_pval=9)

# main loop testing all genes for gain/loss normal and driver
pb = txtProgressBar(min = 0, max = length(cn.data.meta$hgnc_symbol), initial = 0, style = 3) 
for (i in 1:length(cn.data.meta$hgnc_symbol)) {
  setTxtProgressBar(pb,i)
  res.matrix$Gene[i] <- cn.data.meta$hgnc_symbol[i]
  # transpose data # samples geneval pam50 columns
  cn.data <- cn.data.meta[i,] %>% column_to_rownames(var="hgnc_symbol") %>% select(-c("chr","band","pos","genome"))
  cn.data <- as.data.frame(t(cn.data)) %>% rownames_to_column(var="sampleID") 
  cn.data <- merge(cn.data,anno,by="sampleID") #add pam50anno as column
  
  # group data
  LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
  LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
  
  # normal
  # 1. gain
  LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
  LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
  LUMA.data[[2]][LUMA.data[[2]] <= 0] <- 0
  LUMA.data[[2]][LUMA.data[[2]] >= 1] <- 1
  LUMB.data[[2]][LUMB.data[[2]] <= 0] <- 0
  LUMB.data[[2]][LUMB.data[[2]] >= 1] <- 1
  
  # make count tables
  # lumA
  LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
  if(nrow(LUMA.cont.table)==2) {
    res.matrix$Gain_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
  } else {res.matrix$Gain_LUMA_pval[i] <- 1}
  # lumB
  LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
  if(nrow(LUMB.cont.table)==2) {
    res.matrix$Gain_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
  } else {res.matrix$Gain_LUMB_pval[i] <- 1}
  # 2. loss
  LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
  LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
  LUMA.data[[2]][LUMA.data[[2]] <= -1] <- -1
  LUMA.data[[2]][LUMA.data[[2]] >= 0] <- 0
  LUMB.data[[2]][LUMB.data[[2]] <= -1] <- -1
  LUMB.data[[2]][LUMB.data[[2]] >= 0] <- 0
  # lumA
  LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
  if(nrow(LUMA.cont.table)==2) {
    res.matrix$Loss_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
  } else {res.matrix$Loss_LUMA_pval[i] <- 1}
  # lumB
  LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
  if(nrow(LUMB.cont.table)==2) {
    res.matrix$Loss_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
  } else {res.matrix$Loss_LUMB_pval[i] <- 1}
  
  # 3. driver gain
  LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
  LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
  LUMA.data[[2]][LUMA.data[[2]] <= 1] <- 0
  LUMA.data[[2]][LUMA.data[[2]] >= 2] <- 2
  LUMB.data[[2]][LUMB.data[[2]] <= 1] <- 0
  LUMB.data[[2]][LUMB.data[[2]] >= 2] <- 2
  
  # make count tables
  # lumA
  LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
  if(nrow(LUMA.cont.table)==2) {
    res.matrix$DGain_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
  } else {res.matrix$DGain_LUMA_pval[i] <- 1}
  # lumB
  LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
  if(nrow(LUMB.cont.table)==2) {
    res.matrix$DGain_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
  } else {res.matrix$DGain_LUMB_pval[i] <- 1}
  # 4. driver loss
  LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
  LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
  LUMA.data[[2]][LUMA.data[[2]] <= -2] <- -2
  LUMA.data[[2]][LUMA.data[[2]] >= -1] <- 0
  LUMB.data[[2]][LUMB.data[[2]] <= -2] <- -2
  LUMB.data[[2]][LUMB.data[[2]] >= -1] <- 0
  
  # make count tables
  # lumA
  LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
  if(nrow(LUMA.cont.table)==2) {
    res.matrix$DLoss_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
  } else {res.matrix$DLoss_LUMA_pval[i] <- 1}
  # lumB
  LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
  if(nrow(LUMB.cont.table)==2) {
    res.matrix$DLoss_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
  } else {res.matrix$DLoss_LUMB_pval[i] <- 1}
}

# significant genes
# adjust
res.matrix$Loss_LUMA_padj <- p.adjust(
  res.matrix$Loss_LUMA_pval,method = "bonferroni")
res.matrix$Gain_LUMA_padj <- p.adjust(
  res.matrix$Gain_LUMA_pval,method = "bonferroni")
res.matrix$Loss_LUMB_padj <- p.adjust(
  res.matrix$Loss_LUMB_pval,method = "bonferroni")
res.matrix$Gain_LUMB_padj <- p.adjust(
  res.matrix$Gain_LUMB_pval,method = "bonferroni")
res.matrix$DLoss_LUMA_padj <- p.adjust(
  res.matrix$DLoss_LUMA_pval,method = "bonferroni")
res.matrix$DGain_LUMA_padj <- p.adjust(
  res.matrix$DGain_LUMA_pval,method = "bonferroni")
res.matrix$DLoss_LUMB_padj <- p.adjust(
  res.matrix$DLoss_LUMB_pval,method = "bonferroni")
res.matrix$DGain_LUMB_padj <- p.adjust(
  res.matrix$DGain_LUMB_pval,method = "bonferroni")

# modify my.regions to include chromosome name in the cytoband column
my.regions <- my.regions %>% mutate(band = paste(chromosome_name,band,sep=""))
my.regions <- dplyr::rename(my.regions,Gene=hgnc_symbol)
res.matrix <- merge(res.matrix,my.regions,by="Gene") %>% as.data.frame()

save(res.matrix,file = "./data/METABRIC/4_CN/processed/CN_matrix.RData")
#View(res.matrix)

# gain
Gain_LUMA <- res.matrix %>% filter(Gain_LUMA_padj  <= 0.05) %>% dplyr::select(c(Gene,band))
Gain_LUMA_CB <- unique(Gain_LUMA$band)
Gain_LUMB <- res.matrix %>% filter(Gain_LUMB_padj <= 0.05) %>% select(c(Gene,band))
Gain_LUMB_CB <- unique(Gain_LUMB$band)
# loss
Loss_LUMA <- res.matrix %>% filter(Loss_LUMA_padj  <= 0.05) %>% select(c(Gene,band))
Loss_LUMA_CB <- unique(Loss_LUMA$band)
Loss_LUMB <- res.matrix %>% filter(Loss_LUMB_padj  <= 0.05) %>% select(c(Gene,band))
Loss_LUMB_CB <- unique(Loss_LUMB$band)
# gain
DGain_LUMA <- res.matrix %>% filter(DGain_LUMA_padj <= 0.05) %>% select(c(Gene,band))
DGain_LUMA_CB <- unique(DGain_LUMA$band)
DGain_LUMB <- res.matrix %>% filter(DGain_LUMB_padj  <= 0.05) %>% select(c(Gene,band))
DGain_LUMB_CB <- unique(DGain_LUMB$band)
# loss
DLoss_LUMA <- res.matrix %>% filter(DLoss_LUMA_padj <= 0.05) %>% select(c(Gene,band))
DLoss_LUMA_CB <- unique(DLoss_LUMA$band)
DLoss_LUMB <- res.matrix %>% filter(DLoss_LUMB_padj <= 0.05) %>% select(c(Gene,band))
DLoss_LUMB_CB <- unique(DLoss_LUMB$band)