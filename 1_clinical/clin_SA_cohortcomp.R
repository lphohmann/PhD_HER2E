# Script: SA comparison in different SCANB data releases

# TODO: 
# - map release data 1 and 2 to NPJ file (add outcome data as extra columns)
# - select only NPJ follow up TRUE cases for SA
# - do DRFI/OS/RFI for NPJ and OS/RFI for the review data 1&2
# - 

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
# load data and combine into one file
#######################################################################

# SCANB review release 1
load("data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_with_ExternalReview_Bosch_data.RData")
clin.rev1 <- pam50.frame
# SCANB review release 2
clin.rev2 <- read.table("data/SCANB/1_clinical/raw/Project51_ERp_AB_transformed_SPSSselection.txt",sep="\t",fill=TRUE, header = TRUE) %>% filter(str_detect(Case_selected, "^C")) # weird import bug which imported cells with invis characters

# SCANB NPJ release (full)
clin.rel4 <- as.data.frame(read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))
#s <- clin.rel4 %>% filter(Case=="C006808") # checking a sample that was in review file but not follow up NPJ cohort

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
    dplyr::rename(OS_all = OS_days,
                  OSbin_all = OS_event,
                  RFI_all = RFi_days,
                  RFIbin_all = RFi_event,
                  DRFI_all = DRFi_days,
                  DRFIbin_all = DRFi_event,
                  PAM50 = NCN.PAM50) %>% 
    dplyr::select(Case,ER,HER2, 
                  TreatGroup, PAM50,
                  Follow.up.cohort,
                  OS_all, OSbin_all, 
                  RFI_all, RFIbin_all, 
                  DRFI_all, DRFIbin_all) %>% 
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
# 4. Plot for Endo treatment group
#######################################################################

# ET OS rel4
KMplot(group.cohort.version = "ET (cohort: rel4)",
       OMstring = "Overall survival",
       OM = E.data$OS_all,
       OMbin = E.data$OSbin_all,
       sdata = E.data)

# ET OS rev1
KMplot(group.cohort.version = "ET (cohort: rev2)",
       OMstring = "Overall survival",
       OM = E.data$OS_rev2,
       OMbin = E.data$OSbin_rev2,
       sdata = E.data)

# ET OS rev2
KMplot(group.cohort.version = "ET (cohort: rev2)",
       OMstring = "Overall survival",
       OM = E.data$OS_rev2,
       OMbin = E.data$OSbin_rev2,
       sdata = E.data)


# save
ggsave(filename=paste(output.path,cohort,"_RFI_KM_EC.pdf",sep=""), #_basal
       width = 325,
       height = 210,
       units = "mm")


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
