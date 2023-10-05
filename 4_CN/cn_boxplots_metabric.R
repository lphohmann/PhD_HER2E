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
txt.file <- paste(output.path,cohort,"_HER2n_cnplots_tests.txt", sep="")
txt.out <- c() # object to store text output

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
head(cn.freq)

cn.freq <- apply(cn.data, 2, function(x) {(length(x[which(x != 0 & !is.na(x))])/length(x[!is.na(x)]))*100}) %>% 
  as.data.frame() %>% 
  dplyr::rename("freq.altered"=".") %>% 
  rownames_to_column(var="sampleID") 

# add to anno
cn.freq$PAM50 <- anno$PAM50[match(cn.freq$sampleID,anno$sampleID)]

# plot
p <- ggplot(cn.freq, aes(x=PAM50,y=as.numeric(freq.altered),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Genome altered (%)") +
  xlab("PAM50 subtype") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) +
  scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("Her2","LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10))


print(p)
plot.list <- append(plot.list,list(p))

################################################################################
freqs <- cn.freq
# Statistics
# vars
her2e.data <- unname(unlist(subset(
  freqs, PAM50 == "Her2", select=c(freq.altered))))
for (subtype in c("LumA","LumB")) {
  # signature data for subtype
  comp.data <- unname(unlist(subset(freqs, PAM50 == subtype, select=c(freq.altered))))
  
  test.type <- if(diff(range(unlist(her2e.data))) == 0 |
                  diff(range(unlist(comp.data))) == 0 ) {
    # all values identical for one or both -> cant be normally distrib.
    "non-parametric"
  } else if(shapiro.test(unlist(her2e.data))$p.value > 0.05 & 
            shapiro.test(unlist(comp.data))$p.value > 0.05)  {
    # parametric: both normally distributed
    "parametric"
  } else {     
    # at least one group not normally distributed
    "non-parametric"
  }
  
  if(test.type=="parametric") { 
    # equal variance check
    levenes.res <- var.test(unlist(her2e.data), unlist(comp.data), alternative = "two.sided")
    # t.test
    if (!is.nan(levenes.res$p.value) & levenes.res$p.value <= 0.05) {
      res <- t.test(her2e.data, comp.data, var.equal = FALSE)
    } else if (is.nan(levenes.res$p.value)) {
      res <- t.test(her2e.data, comp.data, var.equal = TRUE)
    } else { res <- t.test(her2e.data, comp.data, var.equal = TRUE) }
  } else { 
    # non-parametric
    # mann whitney u test
    res <- wilcox.test(unlist(her2e.data), unlist(comp.data), exact=TRUE)
  }
  # save result#
  print(res)
  txt.out <- append(txt.out,c("Genome altered comparison ",subtype,capture.output(res)))
}

################################################################################
# CN boxplot freq differences of significant genes
################################################################################
cn.data.test <- cn.data %>% na.omit()
#cn.data.test[1:5,1:5]
#row.names(cn.data.test)

her2e.data <- cn.data.test[,colnames(cn.data.test) %in% 
                            anno[anno$PAM50=="Her2",]$sampleID]
luma.data <- cn.data.test[,colnames(cn.data.test) %in% 
                            anno[anno$PAM50=="LumA",]$sampleID]
lumb.data <- cn.data.test[,colnames(cn.data.test) %in% 
                            anno[anno$PAM50=="LumB",]$sampleID]
#View(her2e.dat)

# save obj
# gene LumA.Gain.pval LumB.Gain.pval LumA.Loss.pval LumB.Loss.pval
LumA.Gain.pval <- c()
LumB.Gain.pval <- c()
LumA.Loss.pval <- c()
LumB.Loss.pval <- c()

pb = txtProgressBar(min = 0, max = nrow(cn.data.test), initial = 0, style = 3)
for(i in 1:length(row.names(cn.data.test))) {
  setTxtProgressBar(pb,i)
  # get data
  gene <- row.names(cn.data.test)[i]
  #print(i)
  #print(gene)
  her2e.dat <- as.numeric(her2e.data[gene,])
  luma.dat <- as.numeric(luma.data[gene,])
  lumb.dat <- as.numeric(lumb.data[gene,])
  
  comp.list <- list("LumA"=luma.dat,"LumB"=lumb.dat)
  # test for all pairs
  for(j in 1:length(comp.list)) {
    comp.dat <- comp.list[[j]] 
    comp.group <- names(comp.list)[j]
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
    if (comp.group=="LumA") {
      LumA.Gain.pval[i] <- gain.pval
      LumA.Loss.pval[i] <- loss.pval
    } else if (comp.group=="LumB") {
      LumB.Gain.pval[i] <- gain.pval
      LumB.Loss.pval[i] <- loss.pval
    }
    
  }
  close(pb)
}


gene.test.df <- data.frame("gene"=row.names(cn.data.test), 
                           "LumA.Gain.pval"=LumA.Gain.pval,
                           "LumB.Gain.pval"=LumB.Gain.pval,
                           "LumA.Loss.pval"=LumA.Loss.pval,
                           "LumB.Loss.pval"=LumB.Loss.pval)

# adjust pval
gene.test.df$LumA.Gain.padj <- p.adjust(gene.test.df$LumA.Gain.pval, method = "fdr") # "bonferroni", "fdr"
gene.test.df$LumB.Gain.padj <- p.adjust(gene.test.df$LumB.Gain.pval, method = "fdr")
gene.test.df$LumA.Loss.padj <- p.adjust(gene.test.df$LumA.Loss.pval, method = "fdr")
gene.test.df$LumB.Loss.padj <- p.adjust(gene.test.df$LumB.Loss.pval, method = "fdr")
# View(gene.test.df)
# length(gene.test.df %>% filter(LumA.Gain.padj<=0.05) %>% pull(gene))
# length(gene.test.df %>% filter(LumB.Gain.padj<=0.05) %>% pull(gene))
# length(gene.test.df %>% filter(LumA.Loss.padj<=0.05) %>% pull(gene))
# length(gene.test.df %>% filter(LumB.Loss.padj<=0.05) %>% pull(gene))

#save(gene.test.df, file= "data/METABRIC/4_CN/processed/Metabric_gene_cna_signif.RData")

# for the signfiicant genes, calc the difference in CNfreq and then visualize
genes <- loadRData("data/METABRIC/4_CN/processed/Metabric_gene_cna_signif.RData")
#View(genes)
#colnames(genes)
luma.genes.loss <- genes %>% 
  filter(LumA.Loss.padj <= 0.05) %>% 
  pull(gene)
luma.genes.gain <- genes %>% 
  filter(LumA.Gain.padj <= 0.05) %>% 
  pull(gene)
lumb.genes.loss <- genes %>% 
  filter(LumB.Loss.padj <= 0.05) %>% 
  pull(gene)
lumb.genes.gain <- genes %>% 
  filter(LumB.Gain.padj <= 0.05) %>% 
  pull(gene)

#cn.data.test[1:5,1:5]
gainloss.freqs <- cn.data.test %>% 
  rownames_to_column(var="gene") %>% 
  dplyr::select(gene)

# calc. loss/gain freqs per group
# her2e
gainloss.freqs$freq.loss.her2e <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(her2e.data)], 1, 
  function(x) {
    (length(which(x<=-1))/length(x))*-100
    }
  ) 

gainloss.freqs$freq.gain.her2e  <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(her2e.data)], 1, 
  function(x) {
    (length(which(x>=1))/length(x))*100
  }
)

# luma
gainloss.freqs$freq.loss.luma <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(luma.data)], 1, 
  function(x) {
    (length(which(x<=-1))/length(x))*-100
  }
) 

gainloss.freqs$freq.gain.luma <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(luma.data)], 1, 
  function(x) {
    (length(which(x>=1))/length(x))*100
  }
)

# lumb
gainloss.freqs$freq.loss.lumb <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(lumb.data)], 1, 
  function(x) {
    (length(which(x<=-1))/length(x))*-100
  }
) 

gainloss.freqs$freq.gain.lumb <- apply(
  cn.data.test[,colnames(cn.data.test) %in% colnames(lumb.data)], 1, 
  function(x) {
    (length(which(x>=1))/length(x))*100
  }
)

# save
#save(gainloss.freqs,
#     file = "data/METABRIC/4_CN/processed/CN_gainloss_frequencies.RData")

gene.freqs <- gainloss.freqs

lumb.diff.loss <- gene.freqs %>% 
  filter(gene %in% lumb.genes.loss) %>% 
  mutate(diff = abs(freq.loss.her2e - freq.loss.lumb)) %>% 
  pull(diff)
lumb.diff.gain <- gene.freqs %>% 
  filter(gene %in% lumb.genes.gain) %>% 
  mutate(diff = abs(freq.gain.her2e - freq.gain.lumb)) %>% 
  pull(diff)
lumb.diff <- as.data.frame(c(lumb.diff.loss,lumb.diff.gain)) %>% 
  dplyr::rename("Freq.diff" = "c(lumb.diff.loss, lumb.diff.gain)") %>% 
  mutate(PAM50 = "LumB")

luma.diff.loss <- gene.freqs %>% 
  filter(gene %in% luma.genes.loss) %>% 
  mutate(diff = abs(freq.loss.her2e - freq.loss.luma)) %>% 
  pull(diff)
luma.diff.gain <- gene.freqs %>% 
  filter(gene %in% luma.genes.gain) %>% 
  mutate(diff = abs(freq.gain.her2e - freq.gain.luma)) %>% 
  pull(diff)
luma.diff <- as.data.frame(c(luma.diff.loss,luma.diff.gain)) %>% 
  dplyr::rename("Freq.diff" = "c(luma.diff.loss, luma.diff.gain)") %>% 
  mutate(PAM50 = "LumA")
freq.diffs <- rbind(lumb.diff,luma.diff)

# plot
p <- ggplot(freq.diffs, aes(x=PAM50,y=as.numeric(Freq.diff),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Alteration frequency diff (%)") +
  xlab("PAM50 subtype") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length = unit(0.5, "cm")) +
  scale_fill_manual(values = setNames(c("#2176d5","#34c6eb"),
                                      c("LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10))

#print(p)
#boxplot(Freq.diff~PAM50,data=freq.diffs,col=c("#2176d5","#34c6eb"))
#p <- recordPlot()
#plot.new()
plot.list <- append(plot.list,list(p))









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

chr.lengths <- as.data.frame(read.table(file = "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv", sep = '\t', header = FALSE))[1:23,] %>% 
  dplyr::rename(Chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% # dont shift by 1 here
  filter(Chr != "chrX")

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


# save plots
pdf(file = plot.file, onefile = TRUE)#, height = 5, width = 5)

for(i in 1:length(plot.list)) { 
  print(i)
  print(plot.list[[i]])
}

dev.off()

# save text output
writeLines(txt.out, txt.file)