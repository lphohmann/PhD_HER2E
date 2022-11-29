# Script: Survival analyses in HER2E-FGFR4 expression groups in SCANB

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" 

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

#######################################################################
# Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# load data
clin.rel4 <- as.data.frame(read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))

clin.rel4 <- clin.rel4 %>% 
    filter(Follow.up.cohort==TRUE) %>% 
    filter(HER2 == "Negative") %>% 
    filter(ER == "Positive") %>% 
    filter(NCN.PAM50 == "Her2") %>% 
    dplyr::rename(OS_rel4 = OS_days,
                  OSbin_rel4 = OS_event,
                  RFI_rel4 = RFi_days,
                  RFIbin_rel4 = RFi_event,
                  DRFI_rel4 = DRFi_days,
                  DRFIbin_rel4 = DRFi_event,
                  PAM50 = NCN.PAM50,
                  sampleID = GEX.assay) 

# load FGFR data
load(file = "data/SCANB/2_transcriptomic/processed/fgfr4.gex.RData") #1.66 mean her2e

fgfr4.gex <- fgfr4.gex %>% filter(sampleID %in% clin.rel4$sampleID) %>% dplyr::select(-c(PAM50))

clin.rel4 <- merge(clin.rel4,fgfr4.gex,by="sampleID")

#clin.rel4< - clin.rel4 %>% 
    #mutate(FGFR4= ifelse(FGFR4>=1.66,"high","low"))

clin.rel4 <- clin.rel4  %>%
  mutate(FGFR4_expr = ntile(FGFR4, 3)) %>%
  mutate(FGFR4_expr = if_else(FGFR4_expr == 1, 'Low', if_else(FGFR4_expr == 2, 'Medium', 'High'))) %>%
  arrange(FGFR4)

# plot


pdf(file = paste(output.path,cohort,"_","HER2n_FGFR4_KMplots.pdf", sep=""), onefile = TRUE, width = 21, height = 14.8) 

# OS
# surv object
data.surv <- Surv(clin.rel4$OS_rel4, clin.rel4$OSbin_rel4) 
    
# fit
fit <- survminer::surv_fit(data.surv~FGFR4_expr, data=clin.rel4, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
        fit,
        censor.size = 8,
        censor.shape = "|",
        size = 5,
        risk.table = FALSE,       
        pval = TRUE,
        pval.size = 8,
        pval.coord = c(0,0.1),
        conf.int = FALSE,         
        xlim = c(0,max(clin.rel4$OS_rel4[is.finite(clin.rel4$OS_rel4)])),         
        xlab = "OS days",
        ylab = "OS event probability", # ggf just label as "event probability"
        ylim = c(0,1),
        palette = c("#d334eb", "#2176d5", "#34c6eb"), 
        legend = c(0.9,0.96),
        ggtheme = theme(legend.title = element_text(size=25), #20
                        legend.key.size = unit(0.5,"cm"), 
                        legend.text = element_text(size = 25), #20
                        axis.text.x = element_text(size = 25), #20
                        axis.title.x = element_text(size = 30), #25
                        axis.text.y = element_text(size = 25), #20
                        axis.title.y = element_text(size = 30),
                        plot.title = element_text(size=30)),
        title= "Overall survival in FGFR4 expression groups",
        legend.title = "FGFR4",
        legend.labs = c(paste("High"," (",table(clin.rel4[!is.na(clin.rel4$OS_rel4),]$FGFR4_expr)[1],")",sep = ""),
                        paste("Low"," (",table(clin.rel4[!is.na(clin.rel4$OS_rel4),]$FGFR4_expr)[2],")",sep = ""),
                        paste("Medium"," (",table(clin.rel4[!is.na(clin.rel4$OS_rel4),]$FGFR4_expr)[3],")",sep = "")),
        break.x.by = 500, # break X axis in time intervals of x (what is nicest here? maybe 365)
        break.y.by = 0.1)
    
    print(plot)
    
# RFI
# surv object
data.surv <- Surv(clin.rel4$RFI_rel4, clin.rel4$RFIbin_rel4) 
    
# fit
fit <- survminer::surv_fit(data.surv~FGFR4_expr, data=clin.rel4, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
      fit,
      censor.size = 8,
      censor.shape = "|",
      size = 5,
      risk.table = FALSE,       
      pval = TRUE,
      pval.size = 8,
      pval.coord = c(0,0.1),
      conf.int = FALSE,         
      xlim = c(0,max(clin.rel4$RFI_rel4[is.finite(clin.rel4$RFI_rel4)])),         
      xlab = "RFI days",
      ylab = "RFI event probability", # ggf just label as "event probability"
      ylim = c(0,1),
      palette = c("#d334eb", "#2176d5", "#34c6eb"), 
      legend = c(0.9,0.96),
      ggtheme = theme(legend.title = element_text(size=25), #20
                      legend.key.size = unit(0.5,"cm"), 
                      legend.text = element_text(size = 25), #20
                      axis.text.x = element_text(size = 25), #20
                      axis.title.x = element_text(size = 30), #25
                      axis.text.y = element_text(size = 25), #20
                      axis.title.y = element_text(size = 30),
                      plot.title = element_text(size=30)),
      title= "RFI in FGFR4 expression groups",
      legend.title = "FGFR4",
      legend.labs = c(paste("High"," (",table(clin.rel4[!is.na(clin.rel4$RFI_rel4),]$FGFR4_expr)[1],")",sep = ""),
                      paste("Low"," (",table(clin.rel4[!is.na(clin.rel4$RFI_rel4),]$FGFR4_expr)[2],")",sep = ""),
                      paste("Medium"," (",table(clin.rel4[!is.na(clin.rel4$RFI_rel4),]$FGFR4_expr)[3],")",sep = "")),
      break.x.by = 500, # break X axis in time intervals of x (what is nicest here? maybe 365)
      break.y.by = 0.1)
    
    print(plot)
    

dev.off()

