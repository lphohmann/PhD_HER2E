# Script: Driver mutations in the HER2E subtype (WGS based) - waterfall plot

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
#library(data.table)
library(grid)

#######################################################################
#######################################################################

# load data
sub.drivers <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "SubsDrivers"))) %>% 
  dplyr::select(TumorID_simple,VD_Gene,VC) %>% 
  dplyr::rename(sample=TumorID_simple,gene=VD_Gene,variant_class=VC) %>% 
  mutate(variant_class = paste("sub_",variant_class, sep = ""))

indel.drivers <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "IndelDrivers"))) %>% 
  dplyr::select(TumorID_simple,VD_Gene,VC) %>% 
  dplyr::rename(sample=TumorID_simple,gene=VD_Gene,variant_class=VC) %>% 
  mutate(variant_class = paste("indel_",variant_class, sep = ""))

indel.all <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "AllCodingIndels"))) %>% 
  dplyr::select(Sample,VD_Gene,VC) %>% 
  dplyr::rename(sample=Sample,gene=VD_Gene,variant_class=VC) %>% 
  mutate(variant_class = paste("indel_",variant_class, sep = ""))

sub.all <- as.data.frame((read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "AllCodingSubs"))) %>% 
  dplyr::select(Sample,VD_Gene,VC) %>% 
  dplyr::rename(sample=Sample,gene=VD_Gene,variant_class=VC) %>% 
  mutate(variant_class = paste("sub_",variant_class, sep = ""))

# cn drivers
driver.genes <- unique(c(indel.drivers$gene,sub.drivers$gene))
amp.dat <- loadRData("data/SCANB/4_CN/processed/CN_amp_genpos_genmap.RData")
cn.drivers <- amp.dat[lengths(amp.dat$Gene_symbol) > 0,] %>% 
  filter(!is.na(Gene_symbol)) %>% 
  filter(as.character(Gene_symbol) %in% driver.genes) %>% 
  filter(if_any(starts_with("S"), ~ . > 0))

# final matrix to be filled
cn.drivers.long <- data.frame()

# get it down to gene level
for (sample in colnames(cn.drivers)[grepl("S",colnames(cn.drivers))]) {
  # get sample data
  sample.data <- cn.drivers[c(sample,"ProbeID","Gene_symbol")] #%>% filter(.[[sample]] > 0)
  
  # jump to next iteration if the sample has no driver amps
  if (sort(unique(sample.data[[sample]]))[1]==0 & 
      length(unique(sample.data[[sample]]))==1) { 
    #print(paste("No amps in this sample: ",sample,sep=""))
    next
  }
  #print(paste("AMP IN THIS SAMPLE: ",sample,sep=""))
  for (gene in unique(sample.data$Gene_symbol)) {
     # probes for a gene
     gene.probes <- sample.data %>% 
       filter(Gene_symbol == !!gene) %>% 
       pull("ProbeID")
     probe.count <- sample.data %>% 
       filter(ProbeID %in% gene.probes) %>% 
       filter(.[[sample]] > 0) %>% 
       pull("ProbeID") %>% length(.)
     if (probe.count >= length(gene.probes)/2) {
       cn.drivers.long <- rbind(cn.drivers.long,c(sample,gene,"amplified"))
       
     }
   }
}

# set col names
names(cn.drivers.long) <- c("sample", "gene", "variant_class")

#mutate(variant_class = paste("CN_",variant_class, sep = ""))

################################################################################
# Plotting parameters
################################################################################

# plotting parameters
# 1. mainRecurCutoff accepts a numeric value between 0 and 1, and will only plot genes with mutations in x proportion of samples.
# 2. if there are specific genes of interest those can be specified directly via the plotGenes parameter. Input to plotGenes should be a character vector of a list of genes that are desireable to be shown and is case sensitive. 
# 3. plot only specific samples. This can be achieved via the parameter plotSamples
# 4. the maxGenes parameter will only plot the top x genes and takes an integer value. This is usefull for example if when using the mainRecurCutoff parameter a vector of genes have values at x cutoff and all of them are not desired. 
# 5. the rmvSilent parameter will remove all silent mutations from the data.

# combine
#mut.drivers <- rbind(sub.drivers,indel.drivers) #cn.drivers,
#mut.all <- rbind(sub.all,indel.all)

colors <- c("#0e0421", "#d4136d", "#12e0dd", "#c70c0c", 
                    "#2a18cc", "#0b9c32")

################################################################################
# Plots
################################################################################

datasets <- list(sub.drivers = sub.drivers,
                 indel.drivers = indel.drivers,
                 sub.all = sub.all,
                 indel.all = indel.all)

for (i in names(datasets)) { #1:4
  
  # data
  data <- datasets[[i]] 
  mutation.priority <- as.character(unique(data$variant_class))
  custom.pallete <- colors[1:length(mutation.priority)]
  
  # title
  layer <- list(ggtitle(i))
  
  # plot # idea include all sample but only plot the 25 samples because toherwise the % mutatnt sidebar is not correct in relation to all 30 samples
  plot <- waterfall(data, 
                    fileType = "Custom", 
                    variant_class_order = mutation.priority,
                    mainGrid = TRUE,
                    plotMutBurden = TRUE,
                    mainPalette = custom.pallete,
                    main_geneLabSize = 15,
                    mainRecurCutoff = 0,
                    maxGenes = 20,
                    mainDropMut = TRUE, # drop unused mutation types from legend
                    #rmvSilent = TRUE,
                    out= "grob",
                    mutBurdenLayer = layer)
  #plotSamples = c()
  
  grid.draw(plot)
  
  # append to list
  plot.list <- append(plot.list,list(plot))
  
}

################################################################################

# # data
# data <- list(sub.drivers,indel.drivers,sub.all,indel.all)[[3]] 
# mutation.priority <- as.character(unique(data$variant_class))
# custom.pallete <- colors[1:length(mutation.priority)]
#   
# # plot # idea include all sample but only plot the 25 samples because toherwise the % mutatnt sidebar is not correct in relation to all 30 samples
# plot <- waterfall(data, 
#                   fileType = "Custom", 
#                   variant_class_order = mutation.priority,
#                   mainGrid = TRUE,
#                   plotMutBurden = FALSE,
#                   mainPalette = custom.pallete,
#                   main_geneLabSize = 15,
#                   mainRecurCutoff = 0,
#                   maxGenes = 20,
#                   mainDropMut = TRUE, # drop unused mutation types from legend
#                   #rmvSilent = TRUE,
#                   out= "grob")
# #plotSamples = c()
#   
# grid.draw(plot)
#   
# # append to list
# plot.list <- append(plot.list,list(plot))
  

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE, height = 10, width = 15)

for (i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
  
  #print(plot.list[[i]])
}

dev.off()

