# Script: SA comparison in different SCANB data releases

# TODO: 
# - no endpoints fro release 2 rfi??

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/SCANB_release_comparison/"
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
# load data and combine into one file
#######################################################################

# SCANB review release 1
load("data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_with_ExternalReview_Bosch_data.RData")
clin.rev1 <- pam50.frame

# SCANB review release 2
clin.rev2 <- read.table("data/SCANB/1_clinical/raw/Project51_ERp_AB_transformed_SPSSselection.txt",sep="\t",fill=TRUE, header = TRUE) %>% filter(str_detect(Case_selected, "^C")) # weird import bug which imported cells with invis characters

# I HAVE NO ENDPOINTS FOR REL2 RFI WITHOUT EVENTS SO I SET THEM MYSELF TO THE OS ENDPOINT
#sum(is.na(clin.rev2$Relapse))
clin.rev2 <- clin.rev2 %>% mutate(Days_until_relapse = if_else(Relapse,Days_until_relapse,Days_to_OS))
#clin.rev2$Days_until_relapse[is.na(clin.rev2$Days_until_relapse)] <- max(clin.rev2$Days_until_relapse[is.finite(clin.rev2$Days_until_relapse)]) 

# SCANB NPJ release (full)
clin.rel4 <- as.data.frame(read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))
#s <- clin.rel4 %>% filter(Case=="C006808") # checking a sample that was in review file but not in  follow up NPJ cohort

# add outcome columns to npj
clin.rev1 <- clin.rev1 %>%  
    filter(fuV8 == 1) %>% # to filter out duplicates
    filter(!is.na(treatment_Bosch)) %>% # this step is to only include the review data (just took one of the Bosch columns to check)
    dplyr::rename(
        Case = case_rel4, 
        OS_rev1 = OS, 
        OSbin_rev1 = OSbin, 
        RFIbin_rev1 = relapse_Bosch, 
        RFI_rev1 = relapseTime_Bosch) %>% 
    dplyr::select(Case,OS_rev1,OSbin_rev1,RFI_rev1,RFIbin_rev1) %>% 
    mutate(RFIbin_rev1 = ifelse(RFIbin_rev1,1,0))

clin.rev2 <- clin.rev2 %>% 
    dplyr::rename(
        Case = Case_selected, 
        OS_rev2 = Days_to_OS, 
        OSbin_rev2 = INCA_OS_outcome, 
        RFIbin_rev2 = Relapse, 
        RFI_rev2 = Days_until_relapse) %>% 
    dplyr::select(Case,OS_rev2,OSbin_rev2,RFI_rev2,RFIbin_rev2) %>%  
    mutate(RFIbin_rev2 = ifelse(RFIbin_rev2,1,0)) %>% 
    mutate(OSbin_rev2 = ifelse(OSbin_rev2 == 0,0,1)) # this is weird, why does the event have 3 levels

clin.rel4 <- clin.rel4 %>% 
    dplyr::rename(OS_rel4 = OS_days,
                  OSbin_rel4 = OS_event,
                  RFI_rel4 = RFi_days,
                  RFIbin_rel4 = RFi_event,
                  DRFI_rel4 = DRFi_days,
                  DRFIbin_rel4 = DRFi_event,
                  PAM50 = NCN.PAM50) %>% 
    dplyr::select(Case,ER,HER2, 
                  TreatGroup, PAM50,
                  Follow.up.cohort,
                  OS_rel4, OSbin_rel4, 
                  RFI_rel4, RFIbin_rel4, 
                  DRFI_rel4, DRFIbin_rel4) %>% 
    filter(Follow.up.cohort==TRUE) %>% dplyr::select(-c(Follow.up.cohort))

# checking
# res <- merge(clin.rel4,clin.rev1,by="Case",all=TRUE)
# length(unique(res$Case)) #6660
# res <- merge(clin.rel4,clin.rev2,by="Case",all=TRUE)
# length(unique(res$Case)) #6684 -> cases in review that are not in NPJ follow up cohort
# length(setdiff(clin.rev2$Case,clin.rel4$Case))
# setdiff(clin.rev2$Case,clin.rel4$Case)

# combine
comb.surv.data <- as.data.frame(merge(merge(clin.rel4,clin.rev1,by="Case",all.x=TRUE),clin.rev2,by="Case",all.x=TRUE))
save(comb.surv.data,file = paste(data.path,"surv_data_allreleases.RData",sep=""))

#######################################################################
# Survival Analyses
#######################################################################

load("data/SCANB/1_clinical/processed/surv_data_allreleases.RData")

# select data
comb.surv.data <- comb.surv.data %>% 
    filter(ER=="Positive" & HER2=="Negative") %>% 
    filter(PAM50 %in% c("LumA","LumB","Her2")) 

# create treatment dataframes
E.data <- comb.surv.data %>% filter(TreatGroup =="Endo")
CE.data <- comb.surv.data %>% filter(TreatGroup == "ChemoEndo")

# correct data types
comb.surv.data$PAM50 <- as.factor(comb.surv.data$PAM50)
comb.surv.data$PAM50 <- relevel(comb.surv.data$PAM50, ref = "Her2")
for(i in 6:ncol(comb.surv.data)) {
    comb.surv.data[,i] <- as.numeric(comb.surv.data[,i])
}

#######################################################################
# KM plots for complete SCANB rel4 
#######################################################################
pdf(file = paste(output.path,cohort,"_releases_KMplots.pdf", sep=""), onefile = TRUE, width = 21, height = 14.8) 

# OS
# ET OS rel4
KMplot(group.cohort.version = "ET (cohort: rel4)",
       OMstring = "Overall survival",
       OM = E.data$OS_rel4,
       OMbin = E.data$OSbin_rel4,
       sdata = E.data)

# CT+ET OS rel4
KMplot(group.cohort.version = "CT+ET (cohort: rel4)",
       OMstring = "Overall survival",
       OM = CE.data$OS_rel4,
       OMbin = CE.data$OSbin_rel4,
       sdata = CE.data)

# RFI
# ET RFI rel4
KMplot(group.cohort.version = "ET (cohort: rel4)",
       OMstring = "Recurrence-free interval",
       OM = E.data$RFI_rel4,
       OMbin = E.data$RFIbin_rel4,
       sdata = E.data)

# CT+ET RFI rel4
KMplot(group.cohort.version = "CT+ET (cohort: rel4)",
       OMstring = "Recurrence-free interval",
       OM = CE.data$RFI_rel4,
       OMbin = CE.data$RFIbin_rel4,
       sdata = CE.data)

# DRFI
# ET DRFI rel4
KMplot(group.cohort.version = "ET (cohort: rel4)",
       OMstring = "Distant recurrence-free interval",
       OM = E.data$DRFI_rel4,
       OMbin = E.data$DRFIbin_rel4,
       sdata = E.data)

# CT+ET DRFI rel4
KMplot(group.cohort.version = "CT+ET (cohort: rel4)",
       OMstring = "Distant recurrence-free interval",
       OM = CE.data$DRFI_rel4,
       OMbin = CE.data$DRFIbin_rel4,
       sdata = CE.data)

#######################################################################
# KM plots for complete SCANB clinical review 1 
#######################################################################

# OS
# ET OS rev1
KMplot(group.cohort.version = "ET (cohort: rev1)",
       OMstring = "Overall survival",
       OM = E.data$OS_rev1,
       OMbin = E.data$OSbin_rev1,
       sdata = E.data)

# CT+ET OS rev1
KMplot(group.cohort.version = "CT+ET (cohort: rev1)",
       OMstring = "Overall survival",
       OM = CE.data$OS_rev1,
       OMbin = CE.data$OSbin_rev1,
       sdata = CE.data)

# RFI
# ET RFI rev1
KMplot(group.cohort.version = "ET (cohort: rev1)",
       OMstring = "Recurrence-free interval",
       OM = E.data$RFI_rev1,
       OMbin = E.data$RFIbin_rev1,
       sdata = E.data)

# CT+ET RFI rev1
KMplot(group.cohort.version = "CT+ET (cohort: rev1)",
       OMstring = "Recurrence-free interval",
       OM = CE.data$RFI_rev1,
       OMbin = CE.data$RFIbin_rev1,
       sdata = CE.data)

#######################################################################
# KM plots for complete SCANB clinical review 2
#######################################################################

# OS
# ET OS rev2
KMplot(group.cohort.version = "ET (cohort: rev2)",
       OMstring = "Overall survival",
       OM = E.data$OS_rev2,
       OMbin = E.data$OSbin_rev2,
       sdata = E.data)

# CT+ET OS rev2
KMplot(group.cohort.version = "CT+ET (cohort: rev2)",
       OMstring = "Overall survival",
       OM = CE.data$OS_rev2,
       OMbin = CE.data$OSbin_rev2,
       sdata = CE.data)

# RFI
# ET RFI rev2
KMplot(group.cohort.version = "ET (cohort: rev2)",
       OMstring = "Recurrence-free interval",
       OM = E.data$RFI_rev2,
       OMbin = E.data$RFIbin_rev2,
       sdata = E.data)

# CT+ET RFI rev2
KMplot(group.cohort.version = "CT+ET (cohort: rev2)",
       OMstring = "Recurrence-free interval",
       OM = CE.data$RFI_rev2,
       OMbin = CE.data$RFIbin_rev2,
       sdata = CE.data)


# save
dev.off()