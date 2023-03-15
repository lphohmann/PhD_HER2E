# Script: Driver mutations in the HER2E subtype (WGS based)

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/3_genomic//"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_driverWF.pdf",sep = "")

#packages
source("scripts/3_genomic/src/gen_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(GenVisR)
library(reshape2)
library(ggstatsplot)
library(ggsignif)

#######################################################################
#######################################################################

# load data
sub.drivers <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "SubsDrivers"))) %>% 
  dplyr::select(TumorID_simple,VD_Gene,VC) %>% 
  dplyr::rename(sample=TumorID_simple,gene=VD_Gene,variant_class=VC) %>% 
  mutate(variant_class = paste("sub_",variant_class, sep = ""))
indel.drivers <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "IndelDrivers")))%>% 
  dplyr::select(TumorID_simple,VD_Gene,VC) %>% 
  dplyr::rename(sample=TumorID_simple,gene=VD_Gene,variant_class=VC) %>% 
  mutate(variant_class = paste("indel_",variant_class, sep = ""))

indel.all <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "AllCodingIndels")))
sub.all <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "AllCodingSubs")))

# no rearr and cn drivers in new file
# # rearrangement
# rearr.drivers <- rearr.drivers %>% select(sample,gene2,svclass) %>% 
#     dplyr::rename(gene=gene2,variant_class=svclass) %>% mutate(variant_class = paste("rearr_",variant_class, sep = ""))
# unique(rearr.drivers$variant_class) # "translocation" "deletion" "inversion" "tandem-duplication"
# # copy number
# cn.drivers <- cn.drivers %>% filter(Type == "Amplified") %>% 
#     select(Sample...2,`Gene Symbol`,Type) %>% 
#     dplyr::rename(sample=Sample...2,gene=`Gene Symbol`,variant_class=Type) %>% mutate(variant_class = "Amplification") %>% mutate(variant_class = paste("CN_",variant_class, sep = ""))
# unique(cn.drivers$variant_class) # "Amplified"

#####
# combine
mut.drivers <- rbind(sub.drivers,indel.drivers) #cn.drivers,

length(unique(mut.drivers$gene)) #26


# alternative to waterfall plot
as.character(unique(mut.drivers$variant_class))

mutation_priority <- as.character(unique(mut.drivers$variant_class))

mutation_priority <- c("CN_Amplification","indel_frameshift","indel_inframe","indel_ess_splice","sub_nonsense","sub_missense")


custom_pallete <- c("#0e0421", "#d4136d", "#12e0dd", "#c70c0c", "#2a18cc", "#0b9c32")

#","#a903fc")

# plot 
pdf(file = "~/Desktop/MTP_project/Data/SCAN_B/WGS/waterfall_driver_her2.pdf",width=20, height=10)
waterfall(mut.drivers,fileType = "Custom", variant_class_order=mutation_priority, mainRecurCutoff = 0, maxGenes = 26, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15,mainPalette = custom_pallete)
dev.off()

################################################################################
# Driver waterfall plot etc. (initial apprach only based on subs)
################################################################################

# load the cn data foe one sample
mut.data <- as.data.frame((read_excel("Data/SCAN_B/WGS/CodingAndDriversForJohan_21July22.xlsx")))
View(mut.data)
str(mut.data)
driver.data <- mut.data %>% filter(grepl("yes",mut.data$driver)) %>% 
    select(Sample,VD_Gene,VC) %>% 
    dplyr::rename(sample=Sample,gene=VD_Gene,variant_class=VC)
View(driver.data)

# plot

# Create a vector to save mutation priority order for plotting: CHANGE THIS
mutation_priority <- as.character(unique(driver.data$variant_class))

mutation_priority <- c("nonsense","missense")
# plotting parameters
# 1. mainRecurCutoff accepts a numeric value between 0 and 1, and will only plot genes with mutations in x proportion of samples.
# 2. if there are specific genes of interest those can be specified directly via the plotGenes parameter. Input to plotGenes should be a character vector of a list of genes that are desireable to be shown and is case sensitive. 
# 3. plot only specific samples. This can be achieved via the parameter plotSamples
# 4. the maxGenes parameter will only plot the top x genes and takes an integer value. This is usefull for example if when using the mainRecurCutoff parameter a vector of genes have values at x cutoff and all of them are not desired. 
# 5. the rmvSilent parameter will remove all silent mutations from the data.


# plot for 
# her2 driver
pdf(file = "~/Desktop/MTP_project/Data/SCAN_B/WGS/waterfall_driver_her2.pdf",width=20, height=10)
waterfall(driver.data,fileType = "Custom", variant_class_order=mutation_priority, mainRecurCutoff = 0, maxGenes = 27, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15)
dev.off()

################################################################################
# her2 all wf plot
################################################################################
mut.all <- mut.data %>%
    select(Sample,VD_Gene,VC) %>% 
    dplyr::rename(sample=Sample,gene=VD_Gene,variant_class=VC)

# Create a vector to save mutation priority order for plotting: CHANGE THIS
mutation_priority <- as.character(unique(mut.all$variant_class))

mutation_priority <- c("nonsense","start_lost","stop_lost","missense","ess_splice","5prime_UTR_ess_splice","splice_region","silent")
# make a custom colour pallete
custom_pallete <- c("#0e0421", "#d4136d", "#12e0dd", "#c70c0c", "#2a18cc", "#c7c41e",
                    "#37e019","#a903fc")
# plot for each pam50 subtype
# her2 driver
pdf(file = "~/Desktop/MTP_project/Data/SCAN_B/WGS/waterfall_all_her2.pdf",width=20, height=10)
waterfall(mut.all,fileType = "Custom", variant_class_order=mutation_priority, mainRecurCutoff = 0, maxGenes = 20, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15,rmvSilent = TRUE,mainPalette = custom_pallete)
dev.off()
