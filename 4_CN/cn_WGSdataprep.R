# Script: Preparing analyses based on WGS data

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB 

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_WGS.pdf",sep = "")

#packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)
library(GenVisR)
library(reshape2)
library(ggstatsplot)
#library(data.table)

#######################################################################
#######################################################################

# load data
dat1 <- read_table("./data/SCANB/3_genomic/raw/ascat_easysegments/S000763_easySegments.txt")
dat2 <- read_table("./data/SCANB/3_genomic/raw/ascat_TC_ploidy/S000763_aberrantcellfraction_ploidy.txt")
View(dat1)
View(dat2)

# compile all files into one R object
