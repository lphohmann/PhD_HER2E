# Script: Checking PAM50 correlation

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/1_clinical/processed/",sep="")
dir.create(data.path)

#packages
#source("scripts/1_clinical/src/clin_functions.R")
library(ggplot2)
library(tidyverse)
library(reshape)

#######################################################################
# Load data
#######################################################################

# load data
load("data/SCANB/1_clinical/raw/pam50ann_correlations_summarizedByMean_REL4.RData")
# load anno
clin.rel4 <- as.data.frame(read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))

#######################################################################
# PAM50 centroid correlation in ERpHER2nHER2E samples
#######################################################################

# data processing already done: calc. mean (or median) for each sample for each pam50 group based on the 100 runs

# select subgroup data
anno <- clin.rel4 %>% 
    filter(Follow.up.cohort == TRUE) %>% 
    filter(ER=="Positive" & HER2=="Negative") %>% 
    filter(NCN.PAM50 == "Her2") 

data <- pam50ann_correlations %>% 
    filter(Assay %in% anno$GEX.assay) %>% 
    distinct() #S003571.l.r.m.c.lib.g.k2.a.t occurs 2 times for some reason, the data is identical (duplicated) so I exclude one

# creating modified dataframe
data_mod <- melt(data, measure.vars=colnames(data)[5:9])

# creating boxplot
ggplot(data_mod) + 
    geom_boxplot(aes(x=variable, y=value,fill=as.factor(variable)),alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 Centroid") +
    ylab("Correlation") +
    ggtitle("PAM50 centroid correlation in ERpHER2nHER2E samples") +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 25),
          legend.position = "none") +
    scale_fill_manual(values=c(meanLumA = "#2176d5", meanLumB = "#34c6eb", meanHer2 ="#d334eb", meanBasal ="#524645", meanNormal="#c41b0e"))

ggsave(filename=paste(output.path,cohort,"HER2E_PAM50correlations.pdf",sep=""), 
       width = 260,
       height = 210,
       units = "mm")

#######################################################################
# Crosstable with 2nd best PAM50 match for all samples
#######################################################################

#######################################################################
# Distinctiveness: Boxplot the difference between the best and second 
# best correlation matches
#######################################################################

#######################################################################
# HER2E specific distinctiveness: Boxplot for the her2e samples the difference between the best (her2e) and second best correlation matches
# so three groups in boxplot: basal, luma, lumb, normal
#######################################################################
