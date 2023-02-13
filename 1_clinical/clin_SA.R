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
library(gridExtra)
library(ggplotify)

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

#######################################################################
# Part 1: Figure A
#######################################################################

# list to store plots
plot.list <- list()

# Fig A: IDFS or RFI DRFI; 2x (KM + uniCox) + 2x FP of mvCox
# NOTE FOR SELF: DONT FORGET TO CALC THE MEDIAN FOR CENSORED

######################################
# Investigate the EC treatment group #
######################################

# define the group
EC.dat <- (if(cohort=="SCANB") svdata.rs1 else svdata) %>% 
  filter(Treatment == "CE")
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
#table(EC.dat$PAM50) 

##########################

# uv cox
out <- unicox(EC.dat,EC.surv,title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(group.cohort.version = paste("CT+ET (cohort: ",cohort,")",sep=""),
       OMstring = OM,
       OM = EC.dat[[OM]],
       OMbin = EC.dat[[OMbin]],
       sdata = EC.dat)))

# add the median stuff
median.EC <- median(EC.dat[which(EC.dat[[OMbin]]==0),][[OM]])
median.EC.plot <- ggplot() + annotate("text", x = 2, y = 2, label = paste(
             "Median ",OM, " for CT+ET censored patients = ",median.EC,sep="")) +
  theme_void()
plot.list <- append(plot.list,list(median.EC.plot))

##########################

# mv cox
out <- mvcox(EC.dat, surv=EC.surv,
             title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

#source("scripts/1_clinical/src/clin_functions.R")

########################################
# Investigate the Endo treatment group # CONTINUE HERE TO ADAPT LIKE ABOVE
########################################

# define the group
E.dat <- (if(cohort=="SCANB") svdata.rs1 else svdata) %>% 
  filter(Treatment == "E")
E.surv <- Surv(E.dat[[OM]], E.dat[[OMbin]])
##########################

# uv cox
out <- unicox(E.dat,E.surv,title=paste("Hazard ratios (ERpHER2n, treatment=ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(
  group.cohort.version = paste("ET (cohort: ",cohort,")",sep=""),
  OMstring = OM,
  OM = E.dat[[OM]],
  OMbin = E.dat[[OMbin]],
  sdata = E.dat)))

# add the median stuff
median.E <- median(E.dat[which(E.dat[[OMbin]]==0),][[OM]])
median.E.plot <- ggplot() + annotate("text", x = 2, y = 2, label = paste(
  "Median ",OM, " for ET censored patients = ",median.E,sep="")) +
  theme_void()
plot.list <- append(plot.list,list(median.E.plot))

##########################

# mv cox
out <- mvcox(E.dat, surv=E.surv,
             title=paste("Hazard ratios (ERpHER2n, treatment=ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

#######################################################################
#######################################################################

# save plots
pdf(file = paste(output.path,cohort,"_SA_outcome1.pdf", sep=""), 
    onefile = TRUE, width = 10.5, height = 7.4) 

for (i in 1:length(plot.list)) {
    print(plot.list[[i]])
}

dev.off()

#######################################################################
# Part 2: Figure B
#######################################################################

# - Fig B: OS; 2x (KM + uniCox) + 2x FP of mvCox
# notes: 
# - for SCANB: Fig B based on rel4

# list to store plots
plot.list <- list()

OM <- "OS"
OMbin <- "OSbin"

######################################
# Investigate the EC treatment group #
######################################

# define the group
EC.dat <- (if(cohort=="SCANB") svdata.rel4 else svdata) %>% 
  filter(Treatment == "CE")
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
#table(EC.dat$PAM50) 

##########################

# uv cox
out <- unicox(EC.dat,EC.surv,title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(group.cohort.version = paste("CT+ET (cohort: ",cohort,")",sep=""),
                            OMstring = OM,
                            OM = EC.dat[[OM]],
                            OMbin = EC.dat[[OMbin]],
                            sdata = EC.dat)))

# add the median stuff
median.EC <- median(EC.dat[which(EC.dat[[OMbin]]==0),][[OM]])
median.EC.plot <- ggplot() + annotate("text", x = 2, y = 2, label = paste(
  "Median ",OM, " for CT+ET censored patients = ",median.EC,sep="")) +
  theme_void()
plot.list <- append(plot.list,list(median.EC.plot))

##########################

# mv cox
out <- mvcox(EC.dat, surv=EC.surv,
             title=paste("Hazard ratios (ERpHER2n, treatment=CT+ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

#source("scripts/1_clinical/src/clin_functions.R")

########################################
# Investigate the Endo treatment group # CONTINUE HERE TO ADAPT LIKE ABOVE
########################################

# define the group
E.dat <- (if(cohort=="SCANB") svdata.rel4 else svdata) %>% 
  filter(Treatment == "E")
E.surv <- Surv(E.dat[[OM]], E.dat[[OMbin]])
##########################

# uv cox
out <- unicox(E.dat,E.surv,title=paste("Hazard ratios (ERpHER2n, treatment=ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

##########################

# KM plot
plot.list <- append(plot.list,list(KMplot(
  group.cohort.version = paste("ET (cohort: ",cohort,")",sep=""),
  OMstring = OM,
  OM = E.dat[[OM]],
  OMbin = E.dat[[OMbin]],
  sdata = E.dat)))

# add the median stuff
median.E <- median(E.dat[which(E.dat[[OMbin]]==0),][[OM]])
median.E.plot <- ggplot() + annotate("text", x = 2, y = 2, label = paste(
  "Median ",OM, " for ET censored patients = ",median.E,sep="")) +
  theme_void()
plot.list <- append(plot.list,list(median.E.plot))

##########################

# mv cox
out <- mvcox(E.dat, surv=E.surv,
             title=paste("Hazard ratios (ERpHER2n, treatment=ET, ",OM,", cohort=",cohort,")"))
plot.list <- append(plot.list,list(out$plot))
plot.list <- append(plot.list,list(out$result))

#######################################################################
#######################################################################

# save plots
pdf(file = paste(output.path,cohort,"_SA_outcome2.pdf", sep=""), 
    onefile = TRUE, width = 10.5, height = 7.4) #width = 21, height = 14.8 

for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}

dev.off()

