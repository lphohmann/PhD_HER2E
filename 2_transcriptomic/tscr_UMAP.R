# Script: Generate UMAP

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "METABRIC" # SCANB or METABRIC

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/2_transcriptomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_heatmaps.pdf",sep = "")

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#library(Hmisc)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)
library(amap)

#######################################################################
# testing
#######################################################################

library(ggplot2)
library(umap)

dat <- iris
w<- dat[!duplicated(dat), -5]

umap <- umap(w) # sample=row, vars = columns

df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 Species = dat[!duplicated(dat), 5])

ggplot(df, aes(x, y, colour = Species)) +
  geom_point()

# my data

