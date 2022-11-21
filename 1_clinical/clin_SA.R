# Script: Survival analyses in the Metabric and SCANB review 1 cohort based on RFI 

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # Metabric or SCANB 

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

# list to store plots
plot.list <- list()

#######################################################################
# Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for metabric cohort
if (cohort=="Metabric") {
    # load data
    load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
    anno$Chemotherapy[is.na(anno$Chemotherapy)] <- 0
    anno$Endocrine[is.na(anno$Endocrine)] <- 0
    
    # extract relevant variables
    survdata <- anno %>% 
        dplyr::select(
            METABRIC_ID, Age, Grade, TumSize,
            lymph_nodes_positive, ClinGroup, 
            PAM50, OS, OSbin, DSS, DSSbin, 
            DRFI, DRFIbin, RFI, RFIbin, 
            IDFS, IDFSbin, Chemotherapy, Endocrine) %>% 
        mutate(Treatment = case_when(Chemotherapy==1 & Endocrine==1 ~ "CE",
                                     Chemotherapy==0 & Endocrine==1 ~ "E")) %>% 
        mutate(LN = ifelse(lymph_nodes_positive > 0, "N+", "N0")) %>% 
        dplyr::select(-c(lymph_nodes_positive)) %>% 
        mutate(across(c(OS,DSS,DRFI,RFI,IDFS), (function(years) return(years*365)))) %>% 
        filter(grepl('ERpHER2n', ClinGroup))
    
    # getting correct structure
    survdata$DSS <- as.numeric(survdata$DSS)
    survdata$DSSbin <- as.numeric(survdata$DSSbin)
    survdata$DRFI <- as.numeric(survdata$DRFI)
    survdata$DRFIbin <- as.numeric(survdata$DRFIbin)
    survdata$IDFS <- as.numeric(survdata$IDFS)
    survdata$IDFSbin <- as.numeric(survdata$IDFSbin)
    
# for SCANB cohort
} else if (cohort=="SCANB") {
    # load data
    load("data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_with_ExternalReview_Bosch_data.RData")
    # extract relevant variables
    survdata <- pam50.frame %>% 
        filter(fuV8 == 1) %>% 
        filter(!is.na(treatment_Bosch)) %>% # only include the review data
        mutate(relapse_Bosch = ifelse(relapse_Bosch,1,0)) %>% 
        filter(ERpHER2n_Bosch==1) %>%
        mutate(Treatment = case_when(Chemo_Bosch & ET_Bosch ~ "CE",
                                     !Chemo_Bosch & ET_Bosch ~ "E")) %>% 
        dplyr::rename(
            RFI = relapseTime_Bosch,
            RFIbin = relapse_Bosch, 
            TumSize = TumSize_Bosch,
            PAM50 = PAM50_NCN_rel4,
            Grade = NHG_Bosch,
            LN = LNstatus_Bosch) %>%
        dplyr::select(
            rba_rel4, PAM50, OS, OSbin, TumSize, 
            Age, Grade, LN, RFI, RFIbin, Treatment)
}

# getting correct structure for common variables
survdata$PAM50 <- as.factor(survdata$PAM50)
survdata$Age <- as.numeric(survdata$Age)
survdata$TumSize <- as.numeric(survdata$TumSize)
survdata$Grade <- as.factor(survdata$Grade) 
survdata$LN <- as.factor(survdata$LN) 
survdata$LN <- relevel(survdata$LN, ref = "N0")

# outcome measures
survdata$OS <- as.numeric(survdata$OS)
survdata$OSbin <- as.numeric(survdata$OSbin)
survdata$RFI <- as.numeric(survdata$RFI)
survdata$RFIbin <- as.numeric(survdata$RFIbin)

#######################################################################
# Defining the PAM50 subtypes of interest
#######################################################################

# filter to only include subjects that are PAM50 == Her2 | LumA | LumB
survdata <- survdata %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) # basal
survdata$PAM50 <- droplevels(survdata$PAM50) # drop empty levels
#barplot(table(survdata$PAM50))

#######################################################################
# 3. Relevel to use HER2E as base for comparison
#######################################################################

# relevel and check
levels(survdata$PAM50)
survdata$PAM50 <- relevel(survdata$PAM50, ref = "Her2")
levels(survdata$PAM50)

#######################################################################
# 4. Investigate the EC treatment group
#######################################################################

# define the group
EC_group <- survdata %>% filter(Treatment == "CE")
EC_group.surv <- Surv(EC_group[["RFI"]], EC_group[["RFIbin"]])
#table(EC_group$PAM50) 

##########################

# univariate cox regression + forest plot
plot.list <- append(plot.list,list(unicox(EC_group,EC_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, RFI, cohort=",cohort,")"))))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(group.cohort.version = paste("CT+ET (cohort: ",cohort,")",sep=""),
       OMstring = "Recurrence-free interval",
       OM = EC_group$RFI,
       OMbin = EC_group$RFIbin,
       sdata = EC_group)))

##########################

# Multivariate Cox proportional hazards model
plot.list <- append(plot.list,list(mvcox(EC_group,EC_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, RFI, cohort=",cohort,")"))))

#######################################################################
# 7. Investigate the Endo treatment group
#######################################################################

# define the group
E_group <- survdata %>% filter(Treatment == "E")
E_group.surv <- Surv(E_group[["RFI"]], E_group[["RFIbin"]])
#table(E_group$PAM50)

##########################

# univariate cox regression + forest plot
plot.list <- append(plot.list,list(unicox(E_group,E_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=ET, RFI, cohort=",cohort,")"))))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(group.cohort.version = paste("ET (cohort: ",cohort,")",sep=""),
       OMstring = "Recurrence-free interval",
       OM = E_group$RFI,
       OMbin = E_group$RFIbin,
       sdata = E_group)))

##########################

# Multivariate Cox proportional hazards model
plot.list <- append(plot.list,list(mvcox(E_group,E_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=ET, RFI, cohort=",cohort,")"))))

#######################################################################
#######################################################################

# save plots
pdf(file = paste(output.path,cohort,"_SA.pdf", sep=""), 
    onefile = TRUE, width = 21, height = 14.8) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()
