# Script: Survival analyses in the SCANB cohort for HER2pHER2E vs. HER2p non-HER2E

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
source("scripts/1_clinical/src/clin_functions.R")
library(ggplot2)
library(ggfortify)
library(survival)
library(tidyverse)
library(survminer)
library(grid)
library(readxl)

#######################################################################
# data preprocessing including selection of  
# the clinical HER2p subtyped samples
#######################################################################

# load data
clin.rel4 <- as.data.frame(read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))

clin.rel4 <- clin.rel4 %>% 
    filter(Follow.up.cohort==TRUE) %>% 
    filter(HER2 == "Positive") %>% 
    filter(grepl("Immu",TreatGroup)) %>% 
    dplyr::rename(OS_rel4 = OS_days,
                  OSbin_rel4 = OS_event,
                  RFI_rel4 = RFi_days,
                  RFIbin_rel4 = RFi_event,
                  DRFI_rel4 = DRFi_days,
                  DRFIbin_rel4 = DRFi_event,
                  PAM50 = NCN.PAM50) %>% 
    mutate(PAM50= ifelse(PAM50=="Her2","HER2E","nonHER2E"))
 
#table(clin.rel4$PAM50) 
#######################################################################
# Survival Analyses
#######################################################################

pdf(file = paste(output.path,cohort,"_","HER2p_KMplots.pdf", sep=""), onefile = TRUE, width = 21, height = 14.8) 

# OS
# ET OS rel4
HER2p_KMplot(group.cohort.version = "(HER2p; cohort: rel4)",
       OMstring = "Overall survival",
       OM = clin.rel4$OS_rel4,
       OMbin = clin.rel4$OSbin_rel4,
       sdata = clin.rel4)

# RFI
# ET RFI rel4
HER2p_KMplot(group.cohort.version = "(HER2p; cohort: rel4)",
       OMstring = "Recurrence-free interval",
       OM = clin.rel4$RFI_rel4,
       OMbin = clin.rel4$RFIbin_rel4,
       sdata = clin.rel4)

# DRFI
# ET DRFI rel4
HER2p_KMplot(group.cohort.version = "(HER2p; cohort: rel4)",
       OMstring = "Distant recurrence-free interval",
       OM = clin.rel4$DRFI_rel4,
       OMbin = clin.rel4$DRFIbin_rel4,
       sdata = clin.rel4)

dev.off()
