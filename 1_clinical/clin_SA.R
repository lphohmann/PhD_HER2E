# Script: Survival analyses in the Metabric and SCANB review 1 cohort based on RFI 

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "Metabric" # Metabric or SCANB 

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
    OM="IDFS"# put this IN THE COHORT SPECIFC ONES I.E. IDFS FOR SCANB AND RFI FOR METABRIC
    OMbin="IDFSbin"
    
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
        dplyr::select(-c(lymph_nodes_positive)) %>% #mutate(across(c(OS,DSS,DRFI,RFI,IDFS), (function(years) return(years*365)))) %>% 
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
  OM="IDFS"# put this IN THE COHORT SPECIFC ONES I.E. IDFS FOR SCANB AND RFI FOR METABRIC
  OMbin="IDFSbin"
  # OS from whole release, RFI and RDFS from Bosch review
  
    # # load data
    # load("data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_with_ExternalReview_Bosch_data.RData")
    # # extract relevant variables
    # survdata <- pam50.frame %>% 
    #     filter(fuV8 == 1) %>% 
    #     filter(!is.na(treatment_Bosch)) %>% # only include the review data
    #     mutate(relapse_Bosch = ifelse(relapse_Bosch,1,0)) %>% 
    #     filter(ERpHER2n_Bosch==1) %>%
    #     mutate(Treatment = case_when(Chemo_Bosch & ET_Bosch ~ "CE",
    #                                  !Chemo_Bosch & ET_Bosch ~ "E")) %>% 
    #     dplyr::rename(
    #         RFI = relapseTime_Bosch,
    #         RFIbin = relapse_Bosch, 
    #         TumSize = TumSize_Bosch,
    #         PAM50 = PAM50_NCN_rel4,
    #         Grade = NHG_Bosch,
    #         LN = LNstatus_Bosch) %>%
    #     dplyr::select(
    #         rba_rel4, PAM50, OS, OSbin, TumSize, 
    #         Age, Grade, LN, RFI, RFIbin, Treatment)

  # load data and get right format
  survdata <- loadRData("data/SCANB/1_clinical/processed/SCANB_clinData.RData") %>% 
    filter(Follow.up.cohort == TRUE) %>% 
    filter(!is.na(Bosch_RS1)) %>% # only include the review data
    filter(ERpHER2n_Bosch==1) %>%
    mutate(Treatment = case_when(Chemo_Bosch & ET_Bosch ~ "CE",
                                 !Chemo_Bosch & ET_Bosch ~ "E")) %>% 
    dplyr::select(c("GEX.assay","Sample","Treatment","Age",
                    "ER","PR","HER2","NCN.PAM50","OS","OSbin",
                    "Bosch_RS1","DRFI_bin_Bosch","DRFI_Bosch","RFI_bin_Bosch","RFI_Bosch",
                    "IDFS_bin_Bosch","IDFS_Bosch","ERpHER2n_Bosch",
                    "Chemo_Bosch","ET_Bosch","AI_Bosch","Tamoxifen_Bosch","TumSize_Bosch","NHG_Bosch","LNstatus_Bosch","Ki67_Bosch","Ki67_Bosch_RS2","HER2_Low")) %>% 
    dplyr::rename(
      RFI = RFI_Bosch,
      RFIbin = RFI_bin_Bosch, 
      DRFI = DRFI_Bosch,
      DRFIbin = DRFI_bin_Bosch, 
      IDFS = IDFS_Bosch,
      IDFSbin = IDFS_bin_Bosch, 
      TumSize = TumSize_Bosch,
      PAM50 = NCN.PAM50,
      Grade = NHG_Bosch,
      LN = LNstatus_Bosch) 
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
survdata$IDFS <- as.numeric(survdata$IDFS)
survdata$IDFSbin <- as.numeric(survdata$IDFSbin)

#######################################################################
# Defining the PAM50 subtypes of interest
#######################################################################

# filter to only include subjects that are PAM50 == Her2 | LumA | LumB
survdata <- survdata %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) # basal
survdata$PAM50 <- droplevels(survdata$PAM50) # drop empty levels
#barplot(table(survdata$PAM50))

# define the other comparison groups 
survdata <- survdata %>% 
  mutate(Comp1_groups = if_else(PAM50=="LumA","LUMA","LUMB+HER2E")) %>% 
  mutate(Comp2_groups = if_else(PAM50=="Her2","HER2E","LUMA+LUMB"))

#######################################################################
# 3. Relevel to use HER2E as base for comparison
#######################################################################

# relevel and check
#levels(survdata$Comp2_groups)
survdata$PAM50 <- relevel(survdata$PAM50, ref = "LumA")
survdata$Comp1_groups <- as.factor(survdata$Comp1_groups)
survdata$Comp1_groups <- relevel(survdata$Comp1_groups, ref = "LUMA")
survdata$Comp2_groups <- as.factor(survdata$Comp2_groups)
survdata$Comp2_groups <- relevel(survdata$Comp2_groups, ref = "LUMA+LUMB")

#######################################################################
# 4. Investigate the EC treatment group
#######################################################################

# 4.1 HER2E vs. LUMA vs. LUMB

# define the group
EC_group <- survdata %>% filter(Treatment == "CE")
EC_group.surv <- Surv(EC_group[[OM]], EC_group[[OMbin]])
#table(EC_group$PAM50) 

##########################

# univariate cox regression + forest plot
plot.list <- append(plot.list,list(unicox(EC_group,EC_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, ",OM,", cohort=",cohort,")"))))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(group.cohort.version = paste("Comparison 1: CT+ET (cohort: ",cohort,")",sep=""),
       OMstring = OM,
       OM = EC_group[[OM]],
       OMbin = EC_group[[OMbin]],
       sdata = EC_group)))

##########################

# Multivariate Cox proportional hazards model
plot.list <- append(plot.list,list(mvcox(data=EC_group,
                                         surv=EC_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, ",OM,", cohort=",cohort,")"))))

# 4.2 LUMA vs. HER2E + LUMB 
source("scripts/1_clinical/src/clin_functions.R")

# need to do 2 gorup km function and 2 group mvcox function
plot.list <- append(plot.list,list(
  TwoGroup.KMplot(group.cohort.version = paste("Comparison 2: CT+ET (cohort: ",cohort,")",sep=""),
       OMstring = OM,
       OM = EC_group[[OM]],
       OMbin = EC_group[[OMbin]],
       sdata = EC_group,
       comp.var = "Comp1_groups",
       palette = c("#2176d5", "#bf80ff"),
       legend.labs = 
         c(paste(
            names(table(EC_group[!is.na(OM),]$Comp1_groups)[1]),
            " (",table(EC_group[!is.na(OM),]$Comp1_groups)[1],")",sep=""),
           paste(
             names(table(EC_group[!is.na(OM),]$Comp1_groups)[2]),
             " (",table(EC_group[!is.na(OM),]$Comp1_groups)[2],")",sep="")))))

# 4.2 HER2E vs. LUMA + LUMB 
plot.list <- append(plot.list,list(
  TwoGroup.KMplot(group.cohort.version = paste("Comparison 3: CT+ET (cohort: ",cohort,")",sep=""),
                  OMstring = OM,
                  OM = EC_group[[OM]],
                  OMbin = EC_group[[OMbin]],
                  sdata = EC_group,
                  comp.var = "Comp2_groups",
                  palette = c("#0066ff", "#d334eb"),
                  legend.labs = c(
                    paste(names(table(EC_group[!is.na(OM),]$Comp2_groups)[1])," (",table(EC_group[!is.na(OM),]$Comp2_groups)[1],")",sep=""),
                    paste(names(table(EC_group[!is.na(OM),]$Comp2_groups)[2])," (",table(EC_group[!is.na(OM),]$Comp2_groups)[2],")",sep="")))))

#######################################################################
# 7. Investigate the Endo treatment group
#######################################################################

# define the group
E_group <- survdata %>% filter(Treatment == "E")
E_group.surv <- Surv(E_group[[OM]], E_group[[OMbin]])
#table(E_group$PAM50)

##########################

# univariate cox regression + forest plot
plot.list <- append(plot.list,list(unicox(E_group,E_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=ET, ",OM,", cohort=",cohort,")"))))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(group.cohort.version = paste("Comparison 1: ET (cohort: ",cohort,")",sep=""),
       OMstring = OM,
       OM = E_group[[OM]],
       OMbin = E_group[[OMbin]],
       sdata = E_group)))

##########################

# Multivariate Cox proportional hazards model
plot.list <- append(plot.list,list(mvcox(E_group,E_group.surv,title=paste("Hazard ratios (ERpHER2n, treatment=ET, ",OM,", cohort=",cohort,")"))))

#######################################################################
#######################################################################

# save plots
pdf(file = paste(output.path,cohort,"_SA_new.pdf", sep=""), 
    onefile = TRUE, width = 21, height = 14.8) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()
