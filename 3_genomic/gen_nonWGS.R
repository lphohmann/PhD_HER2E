# Script: Mutational enrichment analyses that are not based on WGS data in SCANB (mutational calling in RNA Seq) and Metabric (panel)

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB or METABRIC

# set/create output directory for plots
output.path <- "output/plots/3_genomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# store plots
pdf(file = paste(output.path,cohort,"_mutfreqs_nonWGS.pdf", sep =""),onefile = TRUE)

#packages
#source("scripts/3_genomic/src/mut_functions.R")
source("scripts/2_transcriptomic/src/tscr_functions.R")
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

#######################################################################
# SCAN-B
#######################################################################

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for scanb cohort
if (cohort=="SCANB") {

# load annotation data
anno <- as.data.frame(
  read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx")) %>%
  filter(Follow.up.cohort==TRUE) %>% 
  filter(NCN.PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
  filter(ER=="Positive" & HER2=="Negative") %>% 
  dplyr::rename(sampleID = Sample, PAM50 = NCN.PAM50)

# load data
mut.data <- read.csv('./data/SCANB/3_genomic/raw/SCANBrel4_ExprSomMutations.csv') %>% 
  as.data.frame(.) %>% 
  column_to_rownames(var = "Gene") %>% 
  dplyr::select(any_of(anno$sampleID))

anno <- anno %>% filter(sampleID %in% names(mut.data))

#######################################################################
# 3. Check which genes significantly differ in their mut pattern between groups
#######################################################################

# groups
HER2.samples <- anno %>% filter(PAM50 == "Her2") %>% pull(sampleID)
LUMA.samples <- anno %>% filter(PAM50 == "LumA") %>% pull(sampleID)
LUMB.samples <- anno %>% filter(PAM50 == "LumB") %>% pull(sampleID)

# pre-filter because it is too slow
# mut.data$total <- rowSums(mut.data == 1)
# mut.data <- mut.data %>% filter(total > 5) %>% dplyr::select(-c("total"))

# genes
genes <- rownames(mut.data)

# objects to store results
res.matrix <- data.frame(matrix(ncol = 3, nrow = length(genes))) %>% dplyr::rename(Gene=1, LUMA_pval=2,LUMB_pval=3)

# merge anno and mut data to create cont tables
mut.anno.data <- as.data.frame(t(mut.data)) %>% rownames_to_column(var="sampleID")
mut.anno.data <- as.data.frame(merge(mut.anno.data,anno,by="sampleID")) %>% dplyr::select(-c(sampleID))
LUMA.mut.data <- subset(mut.anno.data,PAM50 %in% c("Her2","LumA")) %>% droplevels()
LUMB.mut.data <- subset(mut.anno.data,PAM50 %in% c("Her2","LumB")) %>% droplevels()

# main loop testing all genes

pb = txtProgressBar(min = 0, max = length(genes), initial = 0, style = 3) 
for (i in 1:length(genes)) {
    setTxtProgressBar(pb,i)
    # make count tables
    # name results row
    res.matrix$Gene[i] <- genes[i]
    # lumA
    LUMA.cont.table <- table(LUMA.mut.data[[genes[i]]],LUMA.mut.data$PAM50)
    if(nrow(LUMA.cont.table)==2) {
        res.matrix$LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
    } else {res.matrix$LUMA_pval[i] <- 1}
    # lumB
    LUMB.cont.table <- table(LUMB.mut.data[[genes[i]]],LUMB.mut.data$PAM50)
    if(nrow(LUMB.cont.table)==2) {
        res.matrix$LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
    } else {res.matrix$LUMB_pval[i] <- 1}
}

res.matrix$LUMA_padj <- p.adjust(res.matrix$LUMA_pval,method = "bonferroni")
res.matrix$LUMB_padj <- p.adjust(res.matrix$LUMB_pval,method = "bonferroni")

# save
save(res.matrix, file = paste(data.path,"mut_enrichment_result.RData",sep = ""))

LUMA_signif_genes <- res.matrix %>% filter(LUMA_padj <= 0.05) %>% pull(Gene)
LUMB_signif_genes <- res.matrix %>% filter(LUMB_padj <= 0.05) %>% pull(Gene)

#table(LUMB.mut.data[["TP53"]],LUMB.mut.data$PAM50)
#table(LUMA.mut.data[["TP53"]],LUMA.mut.data$PAM50)
#table(mut.anno.data$PAM50)

#######################################################################
# 3. Quick look into ERBB2, AKT1, PTEN, PIK3CA, and TP53
#######################################################################

plot_gene <- function(gene, mut.data) {
    cont.table <- table(mut.data[[gene]],mut.data$PAM50)
    test <- fisher.test(cont.table)
    df <- as.data.frame(countsToCases(as.data.frame(cont.table))) %>% dplyr::rename(PAM50_subtype=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))
    ggbarstats(
        df, Mutated, PAM50_subtype,
        title=paste(gene,"mutation status in SCAN-B",sep = " "),
        results.subtitle = FALSE,
        subtitle = paste0(
            "Fisher's exact test", ", p-value = ",
            ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
        )
    )
}

plot_gene("ERBB2",LUMA.mut.data) 
plot_gene("ERBB2",LUMB.mut.data) 
plot_gene("AKT1",LUMA.mut.data) 
plot_gene("AKT1",LUMB.mut.data) 
plot_gene("PIK3CA",LUMA.mut.data) 
plot_gene("PIK3CA",LUMB.mut.data) 
plot_gene("TP53",LUMA.mut.data) 
plot_gene("TP53",LUMB.mut.data) 
plot_gene("PTEN",LUMA.mut.data) 
plot_gene("PTEN",LUMB.mut.data) 
plot_gene("MTOR",LUMA.mut.data)
plot_gene("MTOR",LUMB.mut.data)
plot_gene("ESR1",LUMA.mut.data) 
plot_gene("ESR1",LUMB.mut.data) 
plot_gene("PTEN",LUMA.mut.data) 
plot_gene("PTEN",LUMB.mut.data) 
plot_gene("PIK3R1",LUMA.mut.data) 
plot_gene("PIK3R1",LUMB.mut.data) 

# make plot with erbb2, akt1, tp53, pik3ca
# columns: gene, pam50, mutcount, total, mutperc

genes <- c("ERBB2","AKT1","PIK3CA","TP53","PIK3R1","ESR1","PTEN")

getformat <- function(gene,mut.anno.data) {
    #gene="ERBB2"
    data <- as.data.frame(matrix(ncol=5,nrow=3,dimnames=list(NULL,c("gene", "pam50", "mutcount", "total", "mutperc"))))
    ct <- as.data.frame(table(mut.anno.data[[gene]],mut.anno.data$PAM50))
    #ct
    data$gene[1:3] <- gene
    data$pam50[1:3] <- c("Her2","LumA","LumB")
    
    for (i in 1:length(data$pam50)) {
        #subtype="LumA"
        #print("Hi")
        subtype <- data$pam50[i]
        data[which(data$pam50==subtype),]$mutcount <- ct[which(ct$Var2==subtype & ct$Var1==1),]$Freq
        data[which(data$pam50==subtype),]$total <- ct[which(ct$Var2==subtype & ct$Var1==1),]$Freq + ct[which(ct$Var2==subtype & ct$Var1==0),]$Freq
        data[which(data$pam50==subtype),]$mutperc <- (data[which(data$pam50==subtype),]$mutcount/data[which(data$pam50==subtype),]$total)*100
        #print(subtype)
    }
    return(data)
}

erbb2.data <- getformat("ERBB2",mut.anno.data)
akt1.data <- getformat("AKT1",mut.anno.data)
tp53.data <- getformat("TP53",mut.anno.data)
pik3ca.data <- getformat("PIK3CA",mut.anno.data)
esr1.data <- getformat("ESR1",mut.anno.data)
pik3r1.data <- getformat("PIK3R1",mut.anno.data)
pten.data <- getformat("PTEN",mut.anno.data)

comb.data <- rbind(erbb2.data,akt1.data,tp53.data,pik3ca.data,esr1.data,pik3r1.data,pten.data)

# erbb2
ggplot(comb.data[which(comb.data$gene=="ERBB2"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,10)) +
    scale_fill_manual(values=c("#d334eb", "#2176d5", "#34c6eb")) +
    ggtitle("ERBB2 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="*", tip_length = 0.02, vjust=0.01, y_position = 9.5, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="*", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

# akt1
ggplot(comb.data[which(comb.data$gene=="AKT1"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,14),
                       breaks=seq(0,12.5,2.5)) +
    scale_fill_manual(values=c("#d334eb", "#2176d5", "#34c6eb")) +
    ggtitle("AKT1 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="*", tip_length = 0.02, vjust=0.01, y_position = 13, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

# pik3cs
ggplot(comb.data[which(comb.data$gene=="PIK3CA"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,60),
                       breaks=seq(0,60,10)) +
    scale_fill_manual(values=c("#d334eb", "#2176d5", "#34c6eb")) +
    ggtitle("PIK3CA mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 58, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="***", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

# tp53
ggplot(comb.data[which(comb.data$gene=="TP53"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,60),
                       breaks=seq(0,60,10)) +
    scale_fill_manual(values=c("#d334eb", "#2176d5", "#34c6eb")) +
    ggtitle("TP53 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="***", tip_length = 0.02, vjust=0.01, y_position = 56, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="***", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

# esr1
ggplot(comb.data[which(comb.data$gene=="ESR1"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,10),
                       breaks=seq(0,10,2.5)) +
    scale_fill_manual(values=c("#d334eb", "#2176d5", "#34c6eb")) +
    ggtitle("ESR1 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 7, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

# pten
ggplot(comb.data[which(comb.data$gene=="PTEN"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,10),
                       breaks=seq(0,10,2.5)) +
    scale_fill_manual(values=c("#d334eb", "#2176d5", "#34c6eb")) +
    ggtitle("PTEN mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 9, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

# pik3r1
ggplot(comb.data[which(comb.data$gene=="PIK3R1"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,5),
                       breaks=seq(0,5,2.5)) +
    scale_fill_manual(values=c("#d334eb", "#2176d5", "#34c6eb")) +
    ggtitle("PIK3R1 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 3, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

#######################################################################
# 3. Quick look into erbb2, pik3ca, akt1, pten - pi3k/akt pathway
#######################################################################
#If you would put erbb2, pik3ca, akt1, pten together, how many cases (%) is affected in your groups?

# define the mutation groups
tp53.pos.samples <- as.data.frame(t(mut.data)) %>% filter(TP53==1)
tp53.pos.samples <- row.names(tp53.pos.samples)
pik3ca.pos.samples <- as.data.frame(t(mut.data)) %>% filter(PIK3CA==1)
pik3ca.pos.samples <- row.names(pik3ca.pos.samples)
erbb2.pos.samples <- as.data.frame(t(mut.data)) %>% filter(ERBB2==1)
erbb2.pos.samples <- row.names(erbb2.pos.samples)
akt1.pos.samples <- as.data.frame(t(mut.data)) %>% filter(AKT1==1)
akt1.pos.samples <- row.names(akt1.pos.samples)
pten.pos.samples <- as.data.frame(t(mut.data)) %>% filter(PTEN==1)
pten.pos.samples <- row.names(pten.pos.samples)

# add to anno
anno <- anno %>% mutate(TP53 = ifelse(sampleID %in% tp53.pos.samples,1,0)) %>% 
    mutate(PIK3CA = ifelse(sampleID %in% pik3ca.pos.samples,1,0)) %>%
    mutate(ERBB2 = ifelse(sampleID %in% erbb2.pos.samples,1,0)) %>%
    mutate(AKT1 = ifelse(sampleID %in% akt1.pos.samples,1,0)) %>%
    mutate(PTEN = ifelse(sampleID %in% pten.pos.samples,1,0))

# make a table with 0,1 if they are mutated in any of these genes
# luma
LUMA.mut.data <- anno %>% mutate(combMUT = ifelse(
  sampleID %in% anno[which(anno$PIK3CA == 1),]$sampleID | 
    sampleID %in% anno[which(anno$ERBB2 == 1),]$sampleID | 
    sampleID %in% anno[which(anno$AKT1 == 1),]$sampleID | 
    sampleID %in% anno[which(anno$PTEN == 1),]$sampleID, 1, 0))
# lumb
LUMB.mut.data <- anno %>% mutate(combMUT = ifelse(
  sampleID %in% anno[which(anno$PIK3CA == 1),]$sampleID | 
    sampleID %in% anno[which(anno$ERBB2 == 1),]$sampleID | 
    sampleID %in% anno[which(anno$AKT1 == 1),]$sampleID | 
    sampleID %in% anno[which(anno$PTEN == 1),]$sampleID, 1, 0))

# lumA
LUMA.cont.table <- table(LUMA.mut.data[["combMUT"]],LUMA.mut.data$PAM50)
atest <- fisher.test(LUMA.cont.table)
adf <- as.data.frame(countsToCases(as.data.frame(LUMA.cont.table))) %>% dplyr::rename(PAM50=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))

# lumB
LUMB.cont.table <- table(LUMB.mut.data[["combMUT"]],LUMB.mut.data$PAM50)
btest <- fisher.test(LUMB.cont.table)
bdf <- as.data.frame(countsToCases(as.data.frame(LUMB.cont.table))) %>% dplyr::rename(PAM50=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))

ggbarstats(
    bdf, Mutated, PAM50,
    title="pi3k/akt pathway mutation status in SCAN-B",
    results.subtitle = FALSE,
    caption = "genes: erbb2,pik3ca,akt1,pten",
    subtitle = paste0(
        "Fisher's exact test", ": LumB pval = ",
        ifelse(btest$p.value < 0.001, "< 0.001", round(btest$p.value, 3)),
        ", LumA pval = ",
        ifelse(atest$p.value < 0.001, "< 0.001", round(atest$p.value, 3)))
    )

#######################################################################
# 3. Quick look into p53 pathway - tp53,mdm2,mdm4,cdkn2a
#######################################################################

# define the mutation groups
tp53.pos.samples <- as.data.frame(t(mut.data)) %>% filter(TP53==1)
tp53.pos.samples <- row.names(tp53.pos.samples)
cdkn2a.pos.samples <- as.data.frame(t(mut.data)) %>% filter(CDKN2A==1)
cdkn2a.pos.samples <- row.names(cdkn2a.pos.samples)
mdm2.pos.samples <- as.data.frame(t(mut.data)) %>% filter(MDM2==1)
mdm2.pos.samples <- row.names(mdm2.pos.samples)
mdm4.pos.samples <- as.data.frame(t(mut.data)) %>% filter(MDM4==1)
mdm4.pos.samples <- row.names(mdm4.pos.samples)

# add to anno
anno <- anno %>% mutate(TP53 = ifelse(sampleID %in% tp53.pos.samples,1,0)) %>% 
    mutate(MDM2 = ifelse(sampleID %in% mdm2.pos.samples,1,0)) %>%
    mutate(MDM4 = ifelse(sampleID %in% mdm4.pos.samples,1,0)) %>%
    mutate(CDKN2A = ifelse(sampleID %in% cdkn2a.pos.samples,1,0))
    
# make a table with 0,1 if they are mutated in any of these genes
# luma
LUMA.mut.data <- anno %>% mutate(combP53 = ifelse(
    sampleID %in% anno[which(anno$TP53 == 1),]$sampleID | 
        sampleID %in% anno[which(anno$MDM2 == 1),]$sampleID | 
        sampleID %in% anno[which(anno$MDM4 == 1),]$sampleID | 
        sampleID %in% anno[which(anno$CDKN2A == 1),]$sampleID, 1, 0))
# lumb
LUMB.mut.data <- anno %>% mutate(combP53 = ifelse(
    sampleID %in% anno[which(anno$TP53 == 1),]$sampleID | 
        sampleID %in% anno[which(anno$MDM2 == 1),]$sampleID | 
        sampleID %in% anno[which(anno$MDM4 == 1),]$sampleID | 
        sampleID %in% anno[which(anno$CDKN2A == 1),]$sampleID, 1, 0))

# plot
# lumA
LUMA.cont.table <- table(LUMA.mut.data[["combP53"]],LUMA.mut.data$PAM50)
atest <- fisher.test(LUMA.cont.table)
adf <- as.data.frame(countsToCases(as.data.frame(LUMA.cont.table))) %>% dplyr::rename(PAM50=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))

# lumB
LUMB.cont.table <- table(LUMB.mut.data[["combP53"]],LUMB.mut.data$PAM50)
btest <- fisher.test(LUMB.cont.table)
bdf <- as.data.frame(countsToCases(as.data.frame(LUMB.cont.table))) %>% dplyr::rename(PAM50=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))

ggbarstats(
    bdf, Mutated, PAM50,
    title="p53 pathway mutation status in SCAN-B",
    results.subtitle = FALSE,
    caption = "genes: tp53,mdm2,mdm4,cdkn2a",
    subtitle = paste0(
      "Fisher's exact test", ": LumB pval = ",
      ifelse(btest$p.value < 0.001, "< 0.001", round(btest$p.value, 3)),
      ", LumA pval = ",
      ifelse(atest$p.value < 0.001, "< 0.001", round(atest$p.value, 3))))

#######################################################################
# 3. Waterfall plots
#######################################################################

# get right input format
# need file with the column sample, gene, variant_class
index <- as.data.frame(which(mut.data==1, arr.ind=TRUE))
matrix.nf <- data.frame(matrix(ncol = 2, nrow = nrow(index))) %>% dplyr::rename(gene=1, sample=2)
pb = txtProgressBar(min = 0, max = nrow(index), initial = 0, style = 3) 
for(i in 1:nrow(index)) {
    setTxtProgressBar(pb,i)
    matrix.nf$gene[i] <- rownames(mut.data)[index$row[i]]
    matrix.nf$sample[i] <- colnames(mut.data)[index$col[i]]
}

# add variant_class column and PAM50 annotation
matrix.nf$variant_class <- "mutation"

# Create a vector to save mutation priority order for plotting
mutation_priority <- as.character(unique(matrix.nf$variant_class))

# plotting parameters
# 1. mainRecurCutoff accepts a numeric value between 0 and 1, and will only plot genes with mutations in x proportion of samples.
# 2. if there are specific genes of interest those can be specified directly via the plotGenes parameter. Input to plotGenes should be a character vector of a list of genes that are desireable to be shown and is case sensitive. 
# 3. plot only specific samples. This can be achieved via the parameter plotSamples
# 4. the maxGenes parameter will only plot the top x genes and takes an integer value. This is usefull for example if when using the mainRecurCutoff parameter a vector of genes have values at x cutoff and all of them are not desired. 
# 5. the rmvSilent parameter will remove all silent mutations from the data.
library(GenVisR)

library(reshape2)
# Melt the clinical data into 'long' format.
clinical.data <- melt(anno, id.vars = c("sample")) # HERERERERE
new_samp_order <- as.character(unique(clinical.data[order(clinical.data$variable, clinical.data$value),]$sample))

waterfall(matrix.nf, fileType = "Custom", variant_class_order=c("mutation"), clinDat = clinical.data, mainRecurCutoff = 0.05, maxGenes = 25, plotMutBurden = FALSE, mainGrid = FALSE,sampOrder = new_samp_order)

# plto for each pam50 subtype
# her2
Her2_samples <- anno %>% filter(PAM50=="Her2") %>% pull(sample)
matrix.nf.her2 <- matrix.nf %>% filter(sample %in% Her2_samples)
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/waterfall_plot_her2.pdf", sep =""),width=20, height=10)
waterfall(matrix.nf.her2, fileType = "Custom", variant_class_order=c("mutation"), clinDat = clinical.data, mainRecurCutoff = 0.05, maxGenes = 25, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15)
dev.off()
# luma
Luma_samples <- anno %>% filter(PAM50=="LumA") %>% pull(sample)
matrix.nf.luma <- matrix.nf %>% filter(sample %in% Luma_samples)
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/waterfall_plot_luma.pdf", sep =""),width=20, height=10)
waterfall(matrix.nf.luma, fileType = "Custom", variant_class_order=c("mutation"), clinDat = clinical.data, mainRecurCutoff = 0.05, maxGenes = 25, plotMutBurden = FALSE, mainGrid = FALSE,main_geneLabSize=15)
dev.off()
#lumb
Lumb_samples <- anno %>% filter(PAM50=="LumB") %>% pull(sample)
matrix.nf.lumb <- matrix.nf %>% filter(sample %in% Lumb_samples)
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/waterfall_plot_lumb.pdf", sep =""),width=20, height=10)
waterfall(matrix.nf.lumb, fileType = "Custom", variant_class_order=c("mutation"), clinDat = clinical.data, mainRecurCutoff = 0.05, maxGenes = 25, plotMutBurden = FALSE, mainGrid = FALSE,main_geneLabSize=15)
dev.off()


#######################################################################
# 4. Investigating P53 and PI3K mutation groups metagene scores
#######################################################################
# define the mutation groups
tp53.pos.samples <- as.data.frame(t(mut.data)) %>% filter(TP53==1)
tp53.pos.samples <- row.names(tp53.pos.samples)
pik3ca.pos.samples <- as.data.frame(t(mut.data)) %>% filter(PIK3CA==1)
pik3ca.pos.samples <- row.names(pik3ca.pos.samples)

anno <- anno %>% mutate(TP53 = ifelse(sampleID %in% tp53.pos.samples,1,0)) %>% mutate(PIK3CA = ifelse(sampleID %in% pik3ca.pos.samples,1,0))


# load metagene scores
load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))

mut.mg_scores <-  scaled_mg_scores %>% rownames_to_column(var="sampleID") %>% 
    mutate(mut_status = case_when(
        sampleID %in% anno[which(anno$TP53 == 1),]$sampleID & sampleID %in% anno[which(anno$PIK3CA == 1),]$sampleID ~ "TP53 & PIK3CA",
        sampleID %in% anno[which(anno$TP53 == 1),]$sampleID & sampleID %in% anno[which(anno$PIK3CA == 0),]$sampleID ~ "TP53",
        sampleID %in% anno[which(anno$TP53 == 0),]$sampleID & sampleID %in% anno[which(anno$PIK3CA == 1),]$sampleID ~ "PIK3CA",
        sampleID %in% anno[which(anno$TP53 == 0),]$sampleID & sampleID %in% anno[which(anno$PIK3CA == 0),]$sampleID ~ "none")) %>% drop_na(mut_status) %>% column_to_rownames(var="sampleID")

# create pdf with all the boxplots
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/tp53_pik3ca_metagene_plots.pdf", sep =""), onefile= TRUE)

for(i in c(1:8)) {
    plot <- ggplot(mut.mg_scores, aes(x=as.factor(mut_status),y=mut.mg_scores[,i])) +
        geom_boxplot(fill="slateblue",alpha=0.2) +
        xlab("Mutation status") +
        ylab(paste(colnames(mut.mg_scores)[i],"metagene score", sep = " ")) +
        ggtitle(colnames(mut.mg_scores)[i]) +
        geom_text(data=as.data.frame(dplyr::count(x=mut.mg_scores, mut_status)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -2,nudge_x = 0.3,size=5) +
        theme(axis.text.x = element_text(size = 15),
              axis.title.x = element_text(size = 20),
              axis.text.y = element_text(size = 15),
              axis.title.y = element_text(size = 20))
    print(plot)
}

dev.off()

#######################################################################
# METABRIC cohort
#######################################################################

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################
#cohort<-"Metabric"
library(janitor)
} else if (cohort=="Metabric") {
cohort="Metabric"
# load data - which one do i need?
# extended
mut.data2 <- as.data.frame(read.delim('Data/Metabric/data_mutations_extended.txt', header = FALSE, sep = "\t", dec = "."))
mut.data2 <- mut.data2[-1,]
mut.data2 <- mut.data2 %>% row_to_names(row_number = 1)
# mskcc
mut.data3 <- as.data.frame(read.delim('Data/Metabric/data_mutations_mskcc.txt', header = FALSE, sep = "\t", dec = "."))
mut.data3 <- mut.data3[-1,]
mut.data3 <- mut.data3 %>% row_to_names(row_number = 1)
#all.equal(mut.data2,mut.data3)
    
# i'll go with extended for now
mut.data <- mut.data2 %>% dplyr::rename(gene=Hugo_Symbol,variant_class=Variant_Classification,sample=Tumor_Sample_Barcode) %>% dplyr::select(sample,gene,variant_class)

# load annotation data
load("Data/Metabric/Annotations/Merged_annotations.RData")
anno <- as.data.frame(anno) %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% filter(grepl('ERpHER2n', ClinGroup)) %>% dplyr::rename(sample=METABRIC_ID) %>% dplyr::select(sample,PAM50)
    
mut.data <- mut.data %>% filter(sample %in% anno$sample)
# samples are not in the mut.data so i exclude them from the anno data as well
# anno <- anno %>% filter(sample %in% mut.data$sample)# this doesnt make sense

#######################################################################
# 3. Waterfall plots
#######################################################################
    
# Create a vector to save mutation priority order for plotting: CHANGE THIS
mutation_priority <- as.character(unique(mut.data$variant_class))

mutation_priority <- c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Splice_Region","Intron","3'UTR","Nonstop_Mutation","Translation_Start_Site","3'Flank","Silent")
# plotting parameters
# 1. mainRecurCutoff accepts a numeric value between 0 and 1, and will only plot genes with mutations in x proportion of samples.
# 2. if there are specific genes of interest those can be specified directly via the plotGenes parameter. Input to plotGenes should be a character vector of a list of genes that are desireable to be shown and is case sensitive. 
# 3. plot only specific samples. This can be achieved via the parameter plotSamples
# 4. the maxGenes parameter will only plot the top x genes and takes an integer value. This is usefull for example if when using the mainRecurCutoff parameter a vector of genes have values at x cutoff and all of them are not desired. 
# 5. the rmvSilent parameter will remove all silent mutations from the data.
    
# clinical data
# Melt the clinical data into 'long' format.
clinical.data <- melt(anno, id.vars = c("sample"))
new_samp_order <- as.character(unique(clinical.data[order(clinical.data$variable, clinical.data$value),]$sample))
    
# plot for each pam50 subtype
# her2
Her2_samples <- anno %>% filter(PAM50=="Her2") %>% pull(sample)
mut.data.her2 <- mut.data %>% filter(sample %in% Her2_samples)
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/waterfall_plot_her2.pdf", sep =""),width=20, height=10)
waterfall(mut.data.her2, fileType = "Custom", variant_class_order=mutation_priority, clinDat = clinical.data, mainRecurCutoff = 0.05, maxGenes = 25, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15)
dev.off()
    
# luma
Luma_samples <- anno %>% filter(PAM50=="LumA") %>% pull(sample)
mut.data.luma <- mut.data %>% filter(sample %in% Luma_samples)
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/waterfall_plot_luma.pdf", sep =""),width=20, height=10)
waterfall(mut.data.luma, fileType = "Custom", variant_class_order=mutation_priority, clinDat = clinical.data, mainRecurCutoff = 0.05, maxGenes = 25, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15)
dev.off()

#lumb
Lumb_samples <- anno %>% filter(PAM50=="LumB") %>% pull(sample)
mut.data.lumb <- mut.data %>% filter(sample %in% Lumb_samples)
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/waterfall_plot_lumb.pdf", sep =""),width=20, height=10)
waterfall(mut.data.lumb,  fileType = "Custom", variant_class_order=mutation_priority, clinDat = clinical.data, mainRecurCutoff = 0.05, maxGenes = 25, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15)
dev.off()

#######################################################################
# 3. Quick look into ERBB2
#######################################################################

#View(head(mut.data))
str(anno)
ERBB2.data <- mut.data %>% filter(gene=="ERBB2")
table(ERBB2.data$variant_class)
    
# create a binary mutation matrix
binary.mut.data <- mut.data %>% mutate(mutation = 1) %>% dplyr::select(-c(variant_class))
binary.mut.data <- distinct(binary.mut.data) # remove duplicates
binary.mut.data$mutation <- NULL
binary.mut.data <- binary.mut.data %>% group_by(gene, sample) %>% ungroup() %>%
    spread(sample, 1, fill = 0) %>% as.data.frame() %>% 
    column_to_rownames(var="gene")
binary.mut.data[binary.mut.data != 0] <- 1
#View(binary.mut.data)
# add the other samples without mutations to it
length(anno$sample)
ncol(binary.mut.data)
    
# merge anno and muta data to create cont tables
mut.anno.data <- as.data.frame(t(binary.mut.data)) %>% rownames_to_column(var="sample")
mut.anno.data <- as.data.frame(merge(mut.anno.data,anno,by="sample")) #%>% dplyr::select(-c(sampleID))
LUMA.mut.data <- subset(mut.anno.data,PAM50 %in% c("Her2","LumA")) %>% droplevels()
LUMB.mut.data <- subset(mut.anno.data,PAM50 %in% c("Her2","LumB")) %>% droplevels()

pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/ERBB2_mutation_status.pdf", sep =""),onefile = TRUE)
    
# lumA
LUMA.cont.table <- table(LUMA.mut.data[["ERBB2"]],LUMA.mut.data$PAM50)
atest <- fisher.test(LUMA.cont.table)
adf <- as.data.frame(countsToCases(as.data.frame(LUMA.cont.table))) %>% dplyr::rename(PAM50_subtype=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))
ggbarstats(
    adf, Mutated, PAM50_subtype,
    title="ERBB2 mutation status in METABRIC",
    results.subtitle = FALSE,
    subtitle = paste0(
        "Fisher's exact test", ", p-value = ",
        ifelse(atest$p.value < 0.001, "< 0.001", round(atest$p.value, 3))
    )
)
    
# lumB
LUMB.cont.table <- table(LUMB.mut.data[["ERBB2"]],LUMB.mut.data$PAM50)
btest <- fisher.test(LUMB.cont.table)
bdf <- as.data.frame(countsToCases(as.data.frame(LUMB.cont.table))) %>% dplyr::rename(PAM50_subtype=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))
ggbarstats(
    bdf, Mutated, PAM50_subtype,
    title="ERBB2 mutation status in METABRIC",
    results.subtitle = FALSE,
    subtitle = paste0(
        "Fisher's exact test", ", p-value = ",
        ifelse(btest$p.value < 0.001, "< 0.001", round(btest$p.value, 3))
    )
)
dev.off()


#######################################################################
# a look into erbb2 etc.
#######################################################################

plot_gene <- function(gene, mut.data) {
    cont.table <- table(mut.data[[gene]],mut.data$PAM50)
    test <- fisher.test(cont.table)
    df <- as.data.frame(countsToCases(as.data.frame(cont.table))) %>% dplyr::rename(PAM50_subtype=Var2,Mutated=Var1) %>% mutate(Mutated = if_else(Mutated==1, "Yes", "No"))
    ggbarstats(
        df, Mutated, PAM50_subtype,
        title=paste(gene,"mutation status in SCAN-B",sep = " "),
        results.subtitle = FALSE,
        subtitle = paste0(
            "Fisher's exact test", ", p-value = ",
            ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
        )
    )
}

#pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/selected_genes_mutation_status.pdf", sep =""),onefile = TRUE)
plot_gene("ERBB2",LUMA.mut.data) # include 0.044 *
plot_gene("ERBB2",LUMB.mut.data) # include 0.041 *
plot_gene("AKT1",LUMA.mut.data) # include 0.59 ns
plot_gene("AKT1",LUMB.mut.data) # include 0.39 ns    0.05 0.01 0.001 0.0001
plot_gene("PIK3CA",LUMA.mut.data) # include 0.16 ns
plot_gene("PIK3CA",LUMB.mut.data) # include 0.08 ns
plot_gene("TP53",LUMA.mut.data) # include 0.001 ***
plot_gene("TP53",LUMB.mut.data) # include 0.001 ***
plot_gene("PTEN",LUMA.mut.data) 
plot_gene("PTEN",LUMB.mut.data) 
plot_gene("MTOR",LUMA.mut.data)
plot_gene("MTOR",LUMB.mut.data)
plot_gene("ESR1",LUMA.mut.data)
plot_gene("ESR1",LUMB.mut.data)


#plot_gene("ESR1",LUMA.mut.data) #ns
#plot_gene("ESR1",LUMB.mut.data) #ns
plot_gene("PTEN",LUMA.mut.data) #ns
plot_gene("PTEN",LUMB.mut.data) #ns
plot_gene("PIK3R1",LUMA.mut.data) #ns
plot_gene("PIK3R1",LUMB.mut.data) #ns
#dev.off()

# make plot with erbb2, akt1, tp53, pik3ca
# columns: gene, pam50, mutcount, total, mutperc
# its ugly but im too tired to think
#genes <- c("ERBB2","AKT1","PIK3CA","TP53")

getformat <- function(gene,mut.anno.data) {
    #gene="ERBB2"
    data <- as.data.frame(matrix(ncol=5,nrow=3,dimnames=list(NULL,c("gene", "pam50", "mutcount", "total", "mutperc"))))
    ct <- as.data.frame(table(mut.anno.data[[gene]],mut.anno.data$PAM50))
    #ct
    data$gene[1:3] <- gene
    data$pam50[1:3] <- c("Her2","LumA","LumB")
    
    for (i in 1:length(data$pam50)) {
        #subtype="LumA"
        #print("Hi")
        subtype <- data$pam50[i]
        data[which(data$pam50==subtype),]$mutcount <- ct[which(ct$Var2==subtype & ct$Var1==1),]$Freq
        data[which(data$pam50==subtype),]$total <- ct[which(ct$Var2==subtype & ct$Var1==1),]$Freq + ct[which(ct$Var2==subtype & ct$Var1==0),]$Freq
        data[which(data$pam50==subtype),]$mutperc <- (data[which(data$pam50==subtype),]$mutcount/data[which(data$pam50==subtype),]$total)*100
        #print(subtype)
    }
    return(data)
}

erbb2.data <- getformat("ERBB2",mut.anno.data)
akt1.data <- getformat("AKT1",mut.anno.data)
tp53.data <- getformat("TP53",mut.anno.data)
pik3ca.data <- getformat("PIK3CA",mut.anno.data)
#esr1.data <- getformat("ESR1",mut.anno.data) # doesnt exist in mb
pik3r1.data <- getformat("PIK3R1",mut.anno.data)
pten.data <- getformat("PTEN",mut.anno.data)
comb.data <- rbind(erbb2.data,akt1.data,tp53.data,pik3ca.data,pik3r1.data,pten.data)
comb.data

#here #############################

# erbb2
ggplot(comb.data[which(comb.data$gene=="ERBB2"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,10)) +
    ggtitle("ERBB2 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="*", tip_length = 0.02, vjust=0.01, y_position = 9.5, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="*", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/ERBB2_mut.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm",)

# akt1
ggplot(comb.data[which(comb.data$gene=="AKT1"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,10),
                       breaks=seq(0,12.5,2.5)) +
    ggtitle("AKT1 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 9.5, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/AKT1_mut.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm",)

# pik3cs
ggplot(comb.data[which(comb.data$gene=="PIK3CA"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,70),
                       breaks = seq(0,70,10)) +
    ggtitle("PIK3CA mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 66, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/PIK3CA_mut.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm",)

# tp53
ggplot(comb.data[which(comb.data$gene=="TP53"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,60),
                       breaks=seq(0,60,10)) +
    ggtitle("TP53 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="***", tip_length = 0.02, vjust=0.01, y_position = 57, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="***", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/TP53_mut.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm",)



# pten
ggplot(comb.data[which(comb.data$gene=="PTEN"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,10),
                       breaks=seq(0,10,2.5)) +
    ggtitle("PTEN mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 7, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/PTEN_mut.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm")


# pik3r1
ggplot(comb.data[which(comb.data$gene=="PIK3R1"),], aes(x=as.factor(pam50),y=mutperc,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Mutation frequency (%)",
                       limits = c(0,5),
                       breaks=seq(0,5,2.5)) +
    ggtitle("PIK3R1 mutation frequency") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 3, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/PIK3R1_mut.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm")








#######################################################################
# 4. Investigating P53 and PI3K mutation groups metagene scores - ADAPT to MB
#######################################################################

# define the mutation groups
tp53.pos.samples <- as.data.frame(t(binary.mut.data)) %>% filter(TP53==1)
tp53.pos.samples <- row.names(tp53.pos.samples)
pik3ca.pos.samples <- as.data.frame(t(binary.mut.data)) %>% filter(PIK3CA==1)
pik3ca.pos.samples <- row.names(pik3ca.pos.samples)

anno <- anno %>% mutate(TP53 = ifelse(sample %in% tp53.pos.samples,1,0)) %>% mutate(PIK3CA = ifelse(sample %in% pik3ca.pos.samples,1,0))

# load metagene scores
load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))

mut.mg_scores <-  scaled_mg_scores %>% rownames_to_column(var="sampleID") %>% 
    mutate(mut_status = case_when(
        sampleID %in% anno[which(anno$TP53 == 1),]$sample & sampleID %in% anno[which(anno$PIK3CA == 1),]$sample ~ "TP53 & PIK3CA",
        sampleID %in% anno[which(anno$TP53 == 1),]$sample & sampleID %in% anno[which(anno$PIK3CA == 0),]$sample ~ "TP53",
        sampleID %in% anno[which(anno$TP53 == 0),]$sample & sampleID %in% anno[which(anno$PIK3CA == 1),]$sample ~ "PIK3CA",
        sampleID %in% anno[which(anno$TP53 == 0),]$sample & sampleID %in% anno[which(anno$PIK3CA == 0),]$sample ~ "none")) %>% drop_na(mut_status) %>% column_to_rownames(var="sampleID")

# create pdf with all the boxplots
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Mutations/",cohort,"/tp53_pik3ca_metagene_plots.pdf", sep =""), onefile= TRUE)

for(i in c(1:8)) {
    plot <- ggplot(mut.mg_scores, aes(x=as.factor(mut_status),y=mut.mg_scores[,i])) +
        geom_boxplot(fill="slateblue",alpha=0.2) +
        xlab("Mutation status") +
        ylab(paste(colnames(mut.mg_scores)[i],"metagene score", sep = " ")) +
        ggtitle(colnames(mut.mg_scores)[i]) +
        geom_text(data=as.data.frame(dplyr::count(x=mut.mg_scores, mut_status)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -2,nudge_x = 0.3,size=5) +
        theme(axis.text.x = element_text(size = 15),
              axis.title.x = element_text(size = 20),
              axis.text.y = element_text(size = 15),
              axis.title.y = element_text(size = 20))
    print(plot)
}

dev.off()

############################### ADAPT STILL ########################################
# 3. Check which genes significantly differ in their mut pattern between groups
#######################################################################

mut.data <- binary.mut.data
View(head(anno))

# groups
HER2.samples <- anno %>% filter(PAM50 == "Her2") %>% pull(sample)
LUMA.samples <- anno %>% filter(PAM50 == "LumA") %>% pull(sample)
LUMB.samples <- anno %>% filter(PAM50 == "LumB") %>% pull(sample)

# pre-filter because it is too slow
mut.data$total <- rowSums(mut.data == 1)
mut.data <- mut.data %>% filter(total > 5) %>% dplyr::select(-c("total"))

# genes
genes <- rownames(mut.data)

# objects to store results
res.matrix <- data.frame(matrix(ncol = 3, nrow = length(genes))) %>% dplyr::rename(Gene=1, LUMA_pval=2,LUMB_pval=3)

# merge anno and muta data to create cont tables
mut.anno.data <- as.data.frame(t(mut.data)) %>% rownames_to_column(var="sample")
mut.anno.data <- as.data.frame(merge(mut.anno.data,anno,by="sample"))
#mut.anno.data <- as.data.frame(t(mut.anno.data)) %>% row_to_names(row_number =1)
LUMA.mut.data <- subset(mut.anno.data,PAM50 %in% c("Her2","LumA")) %>% droplevels()
LUMB.mut.data <- subset(mut.anno.data,PAM50 %in% c("Her2","LumB")) %>% droplevels()

# main loop testing all genes

pb = txtProgressBar(min = 0, max = length(genes), initial = 0, style = 3) 
for (i in 1:length(genes)) {
    setTxtProgressBar(pb,i)
    # make count tables
    # name results row
    res.matrix$Gene[i] <- genes[i]
    # lumA
    LUMA.cont.table <- table(LUMA.mut.data[[genes[i]]],LUMA.mut.data$PAM50)
    if(nrow(LUMA.cont.table)==2) {
        res.matrix$LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
    } else {res.matrix$LUMA_pval[i] <- 1}
    # lumB
    LUMB.cont.table <- table(LUMB.mut.data[[genes[i]]],LUMB.mut.data$PAM50)
    if(nrow(LUMB.cont.table)==2) {
        res.matrix$LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
    } else {res.matrix$LUMB_pval[i] <- 1}
}

# save
save(res.matrix, file = paste("~/Desktop/MTP_project/Output/Mutations/",cohort,"/mut_enrichment_result.RData",sep = ""))


#load(paste("~/Desktop/MTP_project/Output/Mutations/",cohort,"/mut_enrichment_result.RData",sep = ""))
res.matrix$LUMA_padj <- p.adjust(res.matrix$LUMA_pval,method = "fdr")
res.matrix$LUMB_padj <- p.adjust(res.matrix$LUMB_pval,method = "fdr")

LUMA_signif_genes <- res.matrix %>% filter(LUMA_pval <= 0.01)
LUMB_signif_genes <- res.matrix %>% filter(LUMB_pval <= 0.01)

#View(LUMA_signif_genes)
#View(LUMB_signif_genes)
#table(LUMB.mut.data[["TP53"]],LUMB.mut.data$PAM50)
#table(LUMA.mut.data[["TP53"]],LUMA.mut.data$PAM50)
#table(mut.anno.data$PAM50)
# 
} 

dev.off()
