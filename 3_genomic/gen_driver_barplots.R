# Script: Driver mutations barplots for SCANB and BASIS

#TODO:  add ttest for mut freqs

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" 

# set/create output directory for plots
output.path <- "output/plots/3_genomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_mutBarplots.pdf",sep = "")
txt.out <- c() # object to store text output
txt.file <- paste(output.path,cohort,"_HER2n_mutBarplots.txt", sep="")

#packages
source("scripts/3_genomic/src/gen_functions.R")
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(GenVisR)
library(reshape2)
library(ggstatsplot)
#library(data.table)
library(grid)
library(gridExtra)
library(R.utils)
library(vcfR)

#######################################################################
# all mut data
#######################################################################

# scanb
# ID key file
id.key <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "Samples")) %>% 
  dplyr::select(c("Sample","Tumour")) %>% dplyr::rename(sample=Sample)

# get full list from VCF files
scanb.muts <- read_table("data/SCANB/3_genomic/raw/vcf/WGS_epb3d_003_vs_epb3db_003/caveman/epb3d_003_vs_epb3db_003.annot.muts.vcf.gz")
# get normal IDs
# id.key <- read_excel("data/SCANB/3_genomic/raw/JVCimpression2022-11-09_ERpos_WGS_Batch1-3.xlsx") %>% 
#   dplyr::select(c(SENT.TUMOR,SENT.TUMOR.aliquot,TUMOR.alias)) %>% 
#   mutate(SENT.TUMOR.aliquot = gsub("\\.","_",SENT.TUMOR.aliquot))

# compile all files into one R object
caveman.files <- list.files(
  path="./data/SCANB/3_genomic/raw/vcf/",
  pattern="*.annot.muts.vcf.gz$", full.names = TRUE, recursive = TRUE)
pindel.files <- list.files(
  path="./data/SCANB/3_genomic/raw/vcf/",
  pattern="*.annot.vcf.gz$", full.names = TRUE, recursive = TRUE)

# load the files
#caveman.files <- lapply(caveman.files,read.table)  
caveman.files <- lapply(caveman.files, function(x) {
  # load
  file <- read.vcfR(x)
  sample <- gsub(".*caveman/(.+)_vs.*", "\\1", x)
  file.dat <- as.data.frame(file@fix) %>% 
    # make CLPM and ASMD into own column
    mutate(ASMD = as.numeric(gsub(".*ASMD=(.+);VT.*", "\\1", INFO))) %>% 
    mutate(CLPM = as.numeric(gsub(".*CLPM=(.+);ASMD.*", "\\1", INFO))) %>% 
    filter(CLPM==0 & ASMD >=140) #Caveman counts (CLPM=0,ASMD>=140)_final
  
  return(c(sample, nrow(file.dat)))
})  
cv <- as.data.frame(do.call(rbind,caveman.files))
names(cv) <- c("sample","caveman_count")

pindel.files <- lapply(pindel.files, function(x) {
  # load
  file <- read.vcfR(x)
  sample <- gsub(".*pindel/(.+)_vs.*", "\\1", x)
  file.dat <- as.data.frame(file@fix) %>% 
    # make repeats column
    mutate(Repeats = as.numeric(gsub(".*REP=(.+);VT.*", "\\1", INFO))) %>% 
    filter(QUAL>=250 & Repeats>10) #Pindel counts (QUAL>=250,Repeats>10)_final
  
  return(c(sample, nrow(file.dat)))
})
pd <- as.data.frame(do.call(rbind,pindel.files))
names(pd) <- c("sample","pindel_count")

# merge and correct IDs
res <- merge(cv,pd,by="sample")
save(res, file="./data/SCANB/3_genomic/processed/mutation_counts.RData")


# correct ids now






# get down to one file per sample
c(caveman.files, pindel.files)

# convert the file names to S identifiers
# sample ids used in the ascat files
ascat.ids <- as.data.frame(unlist(lapply(temp, function(x) {
  gsub("_vs.*","",
       gsub("./data/SCANB/4_CN/raw/to_lennart//ascat.","",x))}))) %>% dplyr::rename(Ascat.id=1)
# 1. convert epb IDs to normal sample IDs
ascat.ids$Sample <- id.key$SENT.TUMOR[match(
  ascat.ids$Ascat.id, id.key$SENT.TUMOR.aliquot)] 
# 2. convert the other IDs to normal sample IDs
ascat.ids$Sample2 <- id.key$SENT.TUMOR[match(
  ascat.ids$Ascat.id, id.key$TUMOR.alias)] # 2. convert the other IDs to normal sample IDs
ascat.ids$Sample <- ifelse(is.na(ascat.ids$Sample), ascat.ids$Sample2, ascat.ids$Sample)
ascat.ids$Sample2 <- NULL

# name the segments dataframes
names(segment.files) <- ascat.ids$Sample



# basis mut data # HERE GET TEH RIGHT FILE
basis.muts <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  mutate(N_mut = nbrIndels + nbrSubs) %>% 
  dplyr::select(c("sample_name","PAM50_AIMS","N_mut")) %>% 
  dplyr::rename(sample=sample_name,PAM50=PAM50_AIMS)

# combined
mutation.sample.counts <- rbind(basis.muts,scanb.muts)
View(mutation.sample.counts)
#View(basis.muts)

#######################################################################
# plot TMB
#######################################################################

# add back the pam50 annot
#tmb.dat$PAM50 <- all.dmut$PAM50[match(tmb.dat$Sample,all.dmut$Sample)]

p <- ggplot(mutation.sample.counts, aes(x=PAM50,y=N_mut,fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Mutational burden") +
  xlab("PAM50 subtype") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 60,margin = margin(t=10)),
        axis.title.x = element_text(size = 60),
        axis.text.y = element_text(size = 55,margin = margin(r=10)),
        axis.title.y = element_text(size = 60),
        plot.title = element_text(size=50),
        legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) +
  scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("Her2e","LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10)) + 
  ggtitle("TMB based on SBS and indels") +
  


print(p)
plot.list <- append(plot.list,list(p))

#######################################################################
# driver data
#######################################################################

# scanb mut data
scanb.dmut <- loadRData("data/SCANB/3_genomic/processed/driver_mutations_all.RData") %>% 
  mutate(PAM50 = "Her2e") %>% 
  dplyr::rename(Sample=sample, Gene=gene, Effect=variant_class) %>% 
  # split the mutation_type column in two columns (add effect column)
  mutate(Mutation_Type = sub("_$","",str_extract(Effect,"^.*?_"))) %>% 
  mutate(Effect = gsub("^.*?_","",Effect)) %>% 
  mutate(Mutation_Type = if_else(Effect=="amplified","CN",Mutation_Type)) 

# basis mut data
basis.dmut <- as.data.frame(read_excel("data/BASIS/3_genomic/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx", sheet = "COMBINED_EVENTS")) %>% 
  filter(Mutation_Type %in% c("Sub","Ins","Del","CopyNumber")) %>% 
  # split into two rows where #Chr8:(ZNF703/FGFR1)
  mutate(Gene = gsub("^Chr8:","",Gene)) %>% #Chr8:(ZNF703/FGFR1)
  mutate(Gene = gsub("[()]", "", Gene)) %>%
  mutate(Gene = strsplit(Gene, "/")) %>% 
  unnest(Gene) %>% 
  mutate(Mutation_Type = case_when(Mutation_Type=="Sub" ~ "sub",
                                   Mutation_Type=="Ins" ~ "indel",
                                   Mutation_Type=="Del" ~ "indel",
                                   Mutation_Type=="CopyNumber" ~ "CN")) %>% 
  mutate(Effect = if_else(Effect=="amp","amplified",Effect)) %>% 
  dplyr::select(-c(Tumour_name))

# get pam50 annotation
basis.anno <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  dplyr::select(c("sample_name","PAM50_AIMS")) %>% 
  dplyr::rename(Sample=sample_name,PAM50=PAM50_AIMS)


basis.dmut <- merge(basis.anno,basis.dmut,by="Sample")

# combined
all.dmut <- rbind(basis.dmut,scanb.dmut)

#######################################################################
# plot genes of interest: mutation frequencies
#######################################################################

pam50.n <- table(
  all.dmut[!duplicated(all.dmut[,c("Sample")]),]$PAM50) 

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

gene.vec <-c("ERBB2","ESR1","TP53","PIK3CA","FGFR4") 

for (g in gene.vec) {
  
  # Plots freq samples mutated in PAM50 groups
  p2 <- ggplot(count.sample(all.dmut, gene=g), aes(fill=as.factor(PAM50),x=PAM50,y=Freq))+
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

# save plots
pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)

for(i in 1:length(plot.list)) { 
  print(plot.list[[i]])
}

dev.off()


#######################################################################
# stat testing mut. freqs
#######################################################################

# need df in binary format with samples as columns and genes as rows
sample.anno <- all.dmut[1:2] %>% distinct(Sample, .keep_all = TRUE)
all.dmut.binary <- all.dmut %>% 
  dplyr::select(-c(Mutation_Type, Effect, PAM50)) %>% 
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