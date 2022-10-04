# Script: Survival analyses in the Metabric and clinical SCANB cohort based on RFI 

# TODO: - compare different data releases fro SCANB

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
#source("scripts/1_clinical/src/")
library(ggplot2)
library(ggfortify)
library(survival)
library(tidyverse)
library(survminer)
library(grid)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for metabric cohort
if (cohort=="Metabric") {
    # 2.1 load data
    load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
    # 2.2 extract relevant variables
    survdata <- anno %>% dplyr::select(METABRIC_ID, Age, Grade, TumSize,
                                lymph_nodes_positive, ClinGroup, 
                                PAM50, 
                                OS, OSbin, DSS, DSSbin, 
                                DRFI, DRFIbin, RFI, RFIbin, 
                                IDFS, IDFSbin, 
                                Chemotherapy, Endocrine)
    # 2.3 add treatment column
    survdata <- survdata %>% mutate(Treatment = case_when(Chemotherapy==1 & Endocrine==1 ~ "EC",
                                                  Endocrine==1 ~ "E"))
    # 2.4 change lymph node variable to status
    survdata <- survdata %>% mutate(LN = ifelse(lymph_nodes_positive > 0, "N+", "N0"))
    survdata$lymph_nodes_positive <- NULL
    # 2.5 change outcome data to days instead of years
    survdata <- survdata %>% mutate(across(c(OS,DSS,DRFI,RFI,IDFS), (function(years) return(years*365))))
    # 2.6 filter to only include cERpHER2n subjects
    survdata <- survdata %>% filter(grepl('ERpHER2n', ClinGroup))
    # 2.7 getting correct strucutre
    survdata$DSS <- as.numeric(survdata$DSS)
    survdata$DSSbin <- as.numeric(survdata$DSSbin)
    survdata$DRFI <- as.numeric(survdata$DRFI)
    survdata$DRFIbin <- as.numeric(survdata$DRFIbin)
    survdata$IDFS <- as.numeric(survdata$IDFS)
    survdata$IDFSbin <- as.numeric(survdata$IDFSbin)
    
# for SCANB cohort
} else if (cohort=="SCANB") {
    # 2.1 load data
    load("data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_with_ExternalReview_Bosch_data.RData")
    # 2.2 extract relevant variables
    survdata <- pam50.frame %>% dplyr::select(rba_rel4, fuV8,
                                       PAM50_NCN_ProSigna_rel4,
                                       treatment_Bosch, ERpHER2n_Bosch,
                                       Chemo_Bosch, OS, OSbin,
                                       TumSize_Bosch, Age, NHG_Bosch,
                                       LNstatus_Bosch, relapse_Bosch, 
                                       relapseTime_Bosch) 
    # 2.3 rename columns to match metabric annotation
    survdata <- survdata %>% dplyr::rename(RFI = relapseTime_Bosch,
                                    RFIbin = relapse_Bosch, 
                                    TumSize = TumSize_Bosch,
                                    PAM50 = PAM50_NCN_ProSigna_rel4,
                                    Grade = NHG_Bosch,
                                    LN = LNstatus_Bosch)
    # 2.4 filter to only include ERpHER2n subjects suited for outcome analysis
    survdata <- survdata %>% filter(ERpHER2n_Bosch==1 & fuV8==1)
    # 2.5 add treatment column (ASK: IS IT EC OR ONLY C??)
    survdata <- survdata %>% mutate(Treatment = case_when(Chemo_Bosch==1 ~ "EC",
                                                          Chemo_Bosch==0 & treatment_Bosch==1 ~ "E"))
} 

# 2.6 getting correct structure for common variables
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
# 2. Defining the PAM50 subtypes of interest
#######################################################################

# 2.1 filter to only include subjects that are PAM50 == Her2 | LumA | LumB
survdata <- survdata %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) # remove basal
survdata$PAM50 <- droplevels(survdata$PAM50) # drop empty levels
barplot(table(survdata$PAM50))

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
EC_group <- survdata %>% filter(Treatment == "EC")
EC_group.surv <- Surv(EC_group[["RFI"]], EC_group[["RFIbin"]])
table(EC_group$PAM50) 

##########################

# 5.1 Univariate Cox proportional hazards model

# 5.1.1 Model construction
EC_main.pam50 <- coxph(EC_group.surv~PAM50, data=EC_group)

# 5.1.2 Checking assumptions
# mention that prop. Hazard ratio is not fulfilled
cox.zph(EC_main.pam50, transform="km", global=TRUE)
plot(cox.zph(EC_main.pam50, transform="km", global=TRUE))

# 5.1.3 Result
cres <- summary(EC_main.pam50)
round(cres$coefficients[9],5) #lumA
round(cres$coefficients[10],5) #lumB
#ggforest(EC_main.pam50,fontsize = 1)

# 5.1.4 Plot

# add count
# add column n with counts for each group
EC_group <- EC_group %>% add_count(PAM50) #%>% mutate(PAM50_count = paste0(PAM50, ' (', n, ')'))
fit <- survfit(EC_group.surv~PAM50, data=EC_group, conf.type="log-log")
survdiff(EC_group.surv ~ PAM50, data = EC_group) # p=0.05 scnab, p=0.004 metabric

plot <- ggsurvplot(
    fit,
    censor.size = 6,
    censor.shape = "|",
    size = 3,
    risk.table = FALSE,       
    pval = TRUE,
    pval.size = 6,
    pval.coord = c(0,0.45),
    conf.int = FALSE,         
    xlim = c(0,max(EC_group$RFI[is.finite(EC_group$RFI)])),         
    xlab = "Relapse-free interval (days)",
    ylab = "Relapse-free interval probability",
    ylim = c(0.4,1),
    palette = c("#d334eb", "#2176d5", "#34c6eb"), 
    legend = c(0.85,0.90),
    break.time.by = 500,     # break X axis in time intervals by 500.
    ggtheme = theme(legend.title = element_text(size=20), #20
                    legend.key.size = unit(0.5,"cm"), 
                    legend.text = element_text(size = 20), #20
                    axis.text.x = element_text(size = 20), #20
                    axis.title.x = element_text(size = 25), #25
                    axis.text.y = element_text(size = 20), #20
                    axis.title.y = element_text(size = 25),
                    plot.title = element_text(size=22)),
    title="Survival analysis in the chemo + endocrine therapy treatment group",
    legend.title = "Subtypes",
    legend.labs = c(paste("HER2E"," (",(table(EC_group$PAM50)[1]),")",sep = ""),
                    paste("LUMA"," (",(table(EC_group$PAM50)[2]),")",sep = ""),
                    paste("LUMB"," (",(table(EC_group$PAM50)[3]),")",sep = "")),
    break.x.by = 500,
    break.y.by = 0.1)
plot
# save
ggsave(filename=paste(output.path,cohort,"_RFI_KM_EC.pdf",sep=""), #_basal
       width = 325,
       height = 210,
       units = "mm")

# forest 
ggforest(EC_main.pam50,fontsize = 1)

##########################

# 5.2 Multivariate Cox proportional hazards model

# 5.2.1 Model construction 
# parameters to incl: PAM50, Age, Grade, TumSize
EC_main.all <- coxph(EC_group.surv~PAM50+Age+LN+TumSize+Grade, data=EC_group) 

# 5.2.2 Checking assumptions
# Check for violation of proportional hazard (constant HR over time)
cox.zph(EC_main.all, transform="km", global=TRUE) 
plot(cox.zph(EC_main.all, transform="km", global=TRUE))

# 5.2.3 Result
summary(EC_main.all) 

# 5.2.4 Plot forest 
ggforest(EC_main.all,fontsize = 3,cpositions = c(0.01,0.13,0.35))

ggsave(filename=paste(output.path,cohort,"_RFI_forest_EC.pdf",sep=""),
       width = 560,
       height = 480,
       units = "mm")

#######################################################################
# 7. Investigate the Endo treatment group
#######################################################################

# define the group
E_group <- survdata %>% filter(Treatment == "E")
E_group.surv <- Surv(E_group[["RFI"]], E_group[["RFIbin"]])
table(E_group$PAM50)

##########################

# 7.1 Univariate Cox proportional hazards model

# 7.1.1 Model construction
E_main.pam50 <- coxph(E_group.surv~PAM50, data=E_group)

# 7.1.2 Checking assumptions
cox.zph(E_main.pam50, transform="km", global=TRUE)

# 7.1.3 Result
cres <- summary(E_main.pam50)
#cres
round(cres$coefficients[9],5) #luma
round(cres$coefficients[10],5) #lumb

# 7.1.4 Plots

# add column n with counts for each group
E_group <- E_group %>% add_count(PAM50) #%>% mutate(PAM50_count = paste0(PAM50, ' (', n, ')'))

fit <- survfit(E_group.surv~PAM50, data=E_group, conf.type="log-log")
survdiff(E_group.surv ~ PAM50, data = E_group) # p<0.00001 scanb, p= 0.0004 metabric

####################
plot <- ggsurvplot(
    fit,
    censor.size = 6,
    censor.shape = "|",
    size = 3,
    risk.table = FALSE,       
    pval = TRUE,
    pval.size = 6,
    pval.coord = c(0,0.45),
    conf.int = FALSE,         
    xlim = c(0,max(E_group$RFI[is.finite(E_group$RFI)])),     #3500    
    xlab = "Relapse-free interval (days)",
    ylab = "Relapse-free interval probability",
    ylim = c(0.4,1),
    palette = c("#d334eb", "#2176d5", "#34c6eb"), 
    legend = c(0.85,0.90),
    break.time.by = 500,     # break X axis in time intervals by 500.
    ggtheme = theme(legend.title = element_text(size=20), #20
                    legend.key.size = unit(0.5,"cm"), 
                    legend.text = element_text(size = 20), #20
                    axis.text.x = element_text(size = 20), #20
                    axis.title.x = element_text(size = 25), #25
                    axis.text.y = element_text(size = 20), #20
                    axis.title.y = element_text(size = 25),
                    plot.title = element_text(size=22)),
    title="Survival analysis in the endocrine therapy treatment group",
    legend.title = "Subtypes",
    legend.labs = c(paste("HER2E"," (",(table(E_group$PAM50)[1]),")",sep = ""),
                    paste("LUMA"," (",(table(E_group$PAM50)[2]),")",sep = ""),
                    paste("LUMB"," (",(table(E_group$PAM50)[3]),")",sep = "")),
    break.x.by = 500,
    break.y.by = 0.1)

plot

#save
ggsave(filename=paste(output.path,cohort,"_RFI_KM_E.pdf",sep=""), #_basal
       width = 325,
       height = 210,
       units = "mm")

# forest 
ggforest(E_main.pam50,fontsize = 1)

##########################

# 7.2 Multivariate Cox proportional hazards model
# 7.2.1 Model construction
# parameters to incl: PAM50, Age, NHG, LN, Size
E_main.all <- coxph(E_group.surv~PAM50+Age+LN+TumSize+Grade, data=E_group) 

# 7.2.2 Checking assumptions
# Check for violation of proportional hazard (constant HR over time)
cox.zph(E_main.all, transform="km", global=TRUE) 
plot(cox.zph(E_main.all, transform="km", global=TRUE))

# 7.2.3 Result
summary(E_main.all) 

# 7.2.4 Plot forest 
ggforest(E_main.all,fontsize = 3,cpositions = c(0.01,0.13,0.35))
# save
ggsave(filename=paste(output.path,cohort,"_RFI_forest_E.pdf",sep=""),
       width = 560,
       height = 480,
       units = "mm")

#######################################################################
#######################################################################
