# Script: Driver mutations barplots for SCANB and BASIS

#TODO:  

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
library(car)

#######################################################################
# load TMB data (extract from VCF files)
#######################################################################

# scanb
# ID key file
id.key <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "Samples")) %>% 
  dplyr::select(c("Sample","Tumour")) %>% dplyr::rename(sample=Sample)

load("./data/SCANB/3_genomic/processed/mutation_counts.RData")
scanb.muts$caveman_count <- NULL
scanb.muts$pindel_count <- NULL
# basis mut data
basis.muts <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  mutate(N_mut = nbrIndels + nbrSubs) %>% 
  dplyr::select(c("sample_name","PAM50_AIMS","N_mut")) %>% 
  dplyr::rename(sampleID=sample_name,PAM50=PAM50_AIMS)

# combined
mutation.sample.counts <- rbind(basis.muts,scanb.muts)
#mutation.sample.counts <- rbind(basis.muts[c("PAM50","N_mut")],scanb.muts)

#View(mutation.sample.counts)
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
  theme(legend.position = "none",
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) +
  scale_fill_manual(values=setNames(c("#d334eb","#2176d5","#34c6eb"),
                                    c("HER2E","LumA","LumB"))) +  
  scale_y_continuous(breaks = scales::breaks_pretty(10)) + 
  coord_cartesian(ylim=c(0, 12000)) + 
  ggtitle("TMB based on SBS and indels") 
  

pdf(file = plot.file, onefile = TRUE)
plot.list <- append(plot.list,list(p))


#######################################################################
# stat testing TMB 
#######################################################################

# Her2 vs LumA
a.dat <- mutation.sample.counts[which(
  mutation.sample.counts$PAM50=="LumA"),]$N_mut
h.dat <- mutation.sample.counts[which(
  mutation.sample.counts$PAM50=="HER2E"),]$N_mut
# assumptions
#hist(a.dat) # not normal
#hist(h.dat) # not normal
lev.res <- leveneTest(N_mut~as.factor(PAM50),
                      data=mutation.sample.counts[which(
  mutation.sample.counts$PAM50 %in% c("LumA","HER2E")),]) # unequal variances
#luma.res <- t.test(h.dat, a.dat) # assumptions not fulfilled
# Mann-Whitney (Wilcoxon) test
luma.res <- wilcox.test(h.dat, a.dat, exact=TRUE)

# for Her2 vs LumB
b.dat <- mutation.sample.counts[which(
  mutation.sample.counts$PAM50=="LumB"),]$N_mut
h.dat <- mutation.sample.counts[which(
  mutation.sample.counts$PAM50=="HER2E"),]$N_mut
# assumptions
#hist(b.dat) # not normal
#hist(h.dat) # not normal
lev.res <- leveneTest(N_mut~as.factor(PAM50),
                      data=
                        mutation.sample.counts[which(mutation.sample.counts$PAM50 %in% c("LumB","HER2E")),]) # unequal variances
#lumb.res <- t.test(h.dat, b.dat) # assumptions not fulfilled
# Mann-Whitney (Wilcoxon) test
lumb.res <- wilcox.test(h.dat, b.dat, exact=TRUE)
 
# save result
txt.out <- append(txt.out,c("TMB-luma",capture.output(luma.res)))
txt.out <- append(txt.out,c("TMB-lumb",capture.output(lumb.res)))

# save text output
#writeLines(txt.out, txt.file)

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

gene.vec <-c("ERBB2","ESR1","TP53","PIK3CA","FGFR4","MYC") 

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
  txt.out <- append(txt.out, c(g))
  txt.out <- append(txt.out, c(capture.output(
    count.sample(all.dmut, gene=g))))
}

# save plots
#pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)

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



#######################################################################
# check HRD frequencies
#######################################################################

# scanb
HRD.df <- as.data.frame(read_excel("data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "HRDetect")) %>% 
  mutate(HRDetect = ifelse(Probability >= 0.7,"HRD-high","HRD-low")) %>% 
  mutate(PAM50="HER2E") %>% 
  dplyr::select(PAM50,HRDetect)

# basis
HRD.basis <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  dplyr::select(c("PAM50_AIMS","HRDetect")) %>% 
  dplyr::rename(PAM50=PAM50_AIMS)

HRD.all <- rbind(HRD.df,HRD.basis)

# freqs
her2e.dat <- HRD.all[HRD.all$PAM50 == "HER2E",]$HRDetect
her2e.freq <- (length(her2e.dat[her2e.dat == "HRD-high"])/length(her2e.dat))*100
luma.dat <- HRD.all[HRD.all$PAM50 == "LumA",]$HRDetect
luma.freq <- (length(luma.dat[luma.dat == "HRD-high"])/length(luma.dat))*100
lumb.dat <- HRD.all[HRD.all$PAM50 == "LumB",]$HRDetect
lumb.freq <- (length(lumb.dat[lumb.dat == "HRD-high"])/length(lumb.dat))*100

# chi2
txt.out <- append(txt.out,c("Testing HRD Frequency"))
txt.out <- append(txt.out,c("Kataegis frequencies:",
                            "  HER2E =",her2e.freq,
                            "  LumA =",luma.freq,
                            "  LumB =",lumb.freq))

res <- fisher.test(table(
  HRD.all[which(
    HRD.all$PAM50 %in% c("HER2E","LumA")),c("PAM50", "HRDetect")]))

txt.out <- append(txt.out,c(capture.output(res)))
res <- fisher.test(table(
  HRD.all[which(
    HRD.all$PAM50 %in% c("HER2E","LumB")),c("PAM50", "HRDetect")]))
txt.out <- append(txt.out,c(capture.output(res)))












# save text output
writeLines(txt.out, txt.file)

