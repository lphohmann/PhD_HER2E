# Script: Plot Pam50 correlations for Her2 samples  

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/2_transcriptomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots
plot.file <- paste(output.path,cohort,"_HER2n_pam50corr.pdf",sep = "")

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# load data
corr.data <- loadRData("./data/SCANB/1_clinical/raw/SampleSet_WhoAmI_PAM50_n6233_Rel4_with_PAM50correlations_NoUnclassified.RData") %>% dplyr::rename(sampleID=rba)

# filter against my summary object scanb
anno <- loadRData(file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>%
  dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50) %>% 
  filter(PAM50=="Her2")

corr.data <- corr.data %>% filter(sampleID %in% anno$sampleID)

setdiff(anno$sampleID,corr.data$sampleID) # 5 samples missing
#View(head(corr.data))

#######################################################################
# Plot
#######################################################################

# which pam50 subtypes to check correlations for
pam50.subtypes <- c("meanBasal","meanLumA","meanLumB")

for (i in 1:3) {
  
  # select data 
  plot.data <- corr.data[c("majoritySecondBestClass","meanHer2",pam50.subtypes[i])]

  # plot
  plot.list <- append(plot.list,list(
    ggplot(plot.data,aes(y=meanHer2, x=!!sym(pam50.subtypes[i]))) + #,color=majoritySecondBestClass
      geom_point(aes(size=2.5)) +
      geom_abline() +
      ylim(c(0,1)) +
      xlim(c(-1,1)) +
      ylab("Her2 correlation") +
      xlab(paste(gsub('^.{4}', '', pam50.subtypes[i])," correlation",sep="")) +
      ggtitle(paste(
        "PAM50 centroid correlations (PAM50 class: HER2E (n=",nrow(plot.data),")",sep="")) +
      theme(plot.title = element_text(size = 25),
            axis.text.x = element_text(size = 20),
            axis.title.x = element_text(size = 25),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 25),
            legend.position = "none")))
}

# save plots
pdf(file = plot.file, 
    onefile = TRUE, width = 10, height = 10) 

for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}

dev.off()

