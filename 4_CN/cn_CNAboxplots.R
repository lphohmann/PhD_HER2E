# Script: plotting CN boxplot % of probes gain/loss; ploidy vs. PAM50; freq diff gain loss vs. pam50

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
################################################################################

# 1. ploidy

# scanb
scanb.ploidy <- loadRData("./data/SCANB/4_CN/processed/HER2Esample_ploidy.RData") %>% 
  mutate(PAM50 = "Her2e") %>% dplyr::rename(Ploidy=ploidy)

# basis
basis.ploidy <- as.data.frame(read_excel("data/BASIS/4_CN/raw/Supplementary Table 5.Ploidy.AberrantCellFraction.090402015.v1.xlsx"))[1:2]
basis.anno <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  dplyr::select(c("sample_name","PAM50_AIMS")) %>% 
  dplyr::rename(Sample=sample_name,PAM50=PAM50_AIMS)
basis.ploidy <- merge(basis.anno,basis.ploidy,by="Sample")

all.ploidy <- rbind(scanb.ploidy, basis.ploidy)

p <- ggplot(all.ploidy, aes(x=PAM50,y=as.numeric(Ploidy),fill=PAM50)) +
  geom_boxplot(size=2.5, outlier.size = 7) +
  ylab("Tumor ploidy") +
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
  scale_y_continuous(breaks = scales::breaks_pretty(10))


print(p)
plot.list <- append(plot.list,list(p))#

################################################################################
# CN boxplot % of Genes gain/loss 
################################################################################

# WHERE IS THE FILES FOR THE GENES CHECK HOW I GOT THE AMP STATUS
scanb.probes <- loadRData("data/SCANB/4_CN/processed/CN_gainloss_genpos_genmap.RData") %>% 
  dplyr::select(-c(Gene_symbol)) %>% mutate(Chr = gsub("chr","",Chr))
View(head(scanb.probes))

basis.probes <- loadRData("data/BASIS/4_CN/processed/probe_CNA.RData")
basis.probes <- loadRData("data/BASIS/4_CN/summarized")
View(head(basis.probes))



