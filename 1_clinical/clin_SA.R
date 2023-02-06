# Script: Survival analyses in the Metabric and SCANB review 1 cohort based on RFI 

# TODO:

# per cohort output: 
# - Fig A: IDFS or RFI DRFI; 2x (KM + uniCox) + 2x FP of mvCox
# - Fig B: OS; 2x (KM + uniCox) + 2x FP of mvCox
# notes: 
# - add median IDFS for censored pts (bin=0)
# - for SCANB: Fig A based on RS1; Fig B based on rel4
# - for SCANB: mvCox based on RS1 Age, NHG, Size, LN

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
    
    # Fig A OM (Fig B is OS for both)
    OM <- "RFI"
    OMbin <- "RFIbin"
    
    # load data
    svdata <- loadRData("data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData")
    svdata$Chemotherapy[is.na(svdata$Chemotherapy)] <- 0
    svdata$Endocrine[is.na(svdata$Endocrine)] <- 0
    
    # extract relevant variables
    svdata <- svdata %>% 
        mutate(Treatment = case_when(Chemotherapy==1 & Endocrine==1 ~ "CE",
                                     Chemotherapy==0 & Endocrine==1 ~ "E")) %>% 
        mutate(LN = ifelse(lymph_nodes_positive > 0, "N+", "N0")) %>% 
      mutate(across(c(OS,OSbin,RFI,RFIbin,Age,TumSize), as.numeric)) %>% 
      mutate(Grade = factor(Grade, levels = c("1","2","3"))) %>% 
      mutate(PAM50 = factor(PAM50, levels = c("LumA","LumB","Her2"))) %>% 
      mutate(LN = factor(LN, levels = c("N0","N+"))) %>% 
      dplyr::select(
        METABRIC_ID, Age, Grade, TumSize,
        LN, PAM50, OS, OSbin, RFI, RFIbin, 
        Treatment)
    
# for SCANB cohort
} else if (cohort=="SCANB") {
  
  # Fig A OM
  OM <- "IDFS"
  OMbin <- "IDFSbin"
  
  # review 1 data
  svdata.rs1 <- loadRData(file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>% filter(!is.na(Bosch_RS1)) %>%
    mutate(Treatment = case_when(Chemo_Bosch & ET_Bosch ~ "CE",
                                 !Chemo_Bosch & ET_Bosch ~ "E")) %>% 
    dplyr::select(c("GEX.assay","Treatment","Age","NCN.PAM50",
                    "IDFS_bin_Bosch","IDFS_Bosch",
                    "TumSize_Bosch","NHG_Bosch","LNstatus_Bosch")) %>% 
    dplyr::rename(IDFS = IDFS_Bosch, IDFSbin = IDFS_bin_Bosch, 
      TumSize = TumSize_Bosch, PAM50 = NCN.PAM50,
      Grade = NHG_Bosch, LN = LNstatus_Bosch) %>% 
    mutate(across(c(IDFS,IDFSbin,Age,TumSize), as.numeric)) %>% 
    mutate(Grade = factor(Grade, levels = c("1","2","3"))) %>% 
    mutate(PAM50 = factor(PAM50, levels = c("LumA","LumB","Her2"))) %>% 
    mutate(LN = factor(LN, levels = c("N0","N+")))
  
  # rel4 data
  svdata.rel4 <- loadRData(file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>% 
    mutate(Treatment = case_when(TreatGroup=="ChemoEndo" ~ "CE",
                                 TreatGroup=="Endo" ~ "E")) %>% 
    mutate(LN = ifelse(LN > 0, "N+", "N0")) %>%
    dplyr::select(c("GEX.assay","Treatment","Age","NCN.PAM50",
                    "OS","OSbin","Size.mm","NHG","LN")) %>% 
    dplyr::rename(
      TumSize = Size.mm,
      PAM50 = NCN.PAM50,
      Grade = NHG) %>% 
    mutate(across(c(OS,OSbin,Age,TumSize), as.numeric)) %>% 
    mutate(Grade = factor(Grade, levels = c("1","2","3"))) %>% 
    mutate(PAM50 = factor(PAM50, levels = c("LumA","LumB","Her2"))) %>% 
    mutate(LN = factor(LN, levels = c("N0","N+")))
  
}

# getting correct structure for common variables
# survdata$PAM50 <- as.factor(survdata$PAM50)
# survdata$Age <- as.numeric(survdata$Age)
# survdata$TumSize <- as.numeric(survdata$TumSize)
# survdata$Grade <- as.factor(survdata$Grade) 
# survdata$LN <- as.factor(survdata$LN) 
# survdata$LN <- relevel(survdata$LN, ref = "N0")

#######################################################################
# 3. Relevel to use LumA as base for comparison
#######################################################################

# relevel and check
survdata$PAM50 <- relevel(survdata$PAM50, ref = "LumA")


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
