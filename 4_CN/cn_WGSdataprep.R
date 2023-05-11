# Script: Preparing analyses based on WGS data

#TODO: 
# plot to test if the probe states and segment states match
# do the positions of the probes make sense (is it per chr or genome)

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
library(reshape2)
#library(data.table)

#######################################################################
#######################################################################

# compile all files into one R object
# ascat data files
temp <- list.files(
  path="./data/SCANB/3_genomic/raw/ascat_easysegments/",
  pattern="*.txt", full.names = TRUE)
segment.files <- lapply(temp, read_table)
names(segment.files) <- lapply(temp, function(x) {
  gsub("_easySegments.txt","",
       gsub("./data/SCANB/3_genomic/raw/ascat_easysegments//","",x))})
str(segment.files)

# ploidy file
ploidy.df <- list.files(
  path="./data/SCANB/3_genomic/raw/ascat_TC_ploidy/",
  pattern = "*.txt",full.names = TRUE) %>% 
  map_df(~read_table(.)) %>% 
  as.data.frame()
str(ploidy.df)

# probe anno
probe.df <- read.delim("./data/SCANB/3_genomic/raw/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv",sep="\t")[1:3] %>% mutate(Chr = gsub("chr","",Chr)) %>% dplyr::rename("ProbeID"=X)
str(probe.df)

#######################################################################
# functions (move to CN functions later)
#######################################################################

# function that that creates a matrix (probeID | Chr | Position | CNstate_sample) for a all samples
make_finmatrix <- function(probe.df, segment.files, ploidy.df, matrix.type) {
  
  fin.list <- list()
  
  # for each file (sample)
  for(i in c(1:length(segment.files))) {
    
    # sample data
    sample.ID <- names(segment.files[i]) 
    sample.data <- segment.files[[sample.ID]] 
    sample.ploidy <- ploidy.df[which(ploidy.df$Sample==sample.ID),]$ploidy
    
    # create a matrix (probeID | Chr | Position | CNstate_sample) for a single sample
    sample.matrix <- make_sampmatrix(probe.df, sample.data, sample.ID, matrix.type, sample.ploidy)
    
    # list with sample matrices
    fin.list <- append(fin.list, list(sample.matrix))
  }
            
  # merge fin.list into final matrix
  
  # reduce(lambda x, y: pd.merge(x, y, on = ['ProbeID','Chr','Position']), fin_list)
  #fin.matrix <- 
  
  return(fin.matrix)
}

#######################################################################

## function that that creates a matrix (probeID | Chr | Position | CNstate_sample) for a single sample
make_sampmatrix <- function(probe.df, sample.data, sample.ID, matrix.type, sample.ploidy) {
  
  chr.set <- unique(sample.data$chr)
  res.list <- list()
  
  for(j in c(1:length(chr.set))) { 
    # chromosome
    chr <- chr.set[j] #j
  # relevant probes
    chr.probes <- probe.df[probe.df$Chr == chr,]
    # relevant segments
    chr.segments <- sample.data[sample.data$chr == chr,]
    # get chr matrix
    chr.matrix <- make_chrmatrix(
      chr.probes, chr.segments, matrix.type, sample.ploidy)
    
    # save in list
    res.list <- append(res.list,list(chr.matrix))
  }

  #make one dataframe out of the list
  sample.matrix <- do.call(rbind, res.list)

  # rename cn_state to the sample id
  sample.matrix <- sample.matrix %>% dplyr::rename(!!sample.ID := CN_state)
                      
  return(sample.matrix) 
}

#######################################################################

# function that creates a CN matrix (probeID | Chr | Position | CNstate_sample) for a single chromosome with a method defining what type of matrix is created (total, major, minor)
make_chrmatrix <- function(chr.probes, chr.segments, matrix.type, sample.ploidy) {
  
  # add a column to probe file that is to be filled with the CN state for each probe
  chr.matrix <- chr.probes
  chr.matrix$CN_state <- NA 
  
  # for each chromosome segment
  for(k in c(1:nrow(chr.segments))) {
    
    # get sample info
    seg.data <- chr.segments[k,] 
    seg.end <- seg.data[["endpos"]]
    seg.start <- seg.data[["startpos"]]
    
    # define calculations for each matrix type
    # major CN state
    if (matrix.type == "major") {
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- seg.data[["nMajor"]]
    # minor CN state
      } else if (matrix.type == "minor") {
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- seg.data[["nMinor"]]
    # total CN state  
      } else if (matrix.type == "total") {
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- seg.data[["nTot"]]
    # gain or loss CN state
      } else if (matrix.type == "gainloss") {
        if (seg.data[["nTot"]]  >= (sample.ploidy + 0.6)) {
          state.gl <- 1
          } else if (seg.data[["nTot"]] <= (sample.ploidy - 0.6)) {
            state.gl <- -1
          } else {
            state.gl <- 0
        }
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- state.gl
    # amplification CN state
      } else if (matrix.type == "amplification") {
        if (seg.data[["nTot"]]  >= (sample.ploidy*4)) { # focal amp
          state.amp <- 2
          } else if (seg.data[["nTot"]] >= (sample.ploidy*2)) { # amp
            state.amp <- 1
          } else {
            state.amp <- 0 # no amp
        }
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- state.amp
    # LOH CN state
      } else if (matrix.type == "LOH") {
        if (seg.data[["nMinor"]]  == 0) { # LOH
          state.loh <- 1
          } else {
            state.loh <- 0 # no LOH
          }
        # assign the CN_state to the probes that correspond to that segment
        chr.matrix[which(chr.matrix$Position >= seg.start & 
                           chr.matrix$Position <= seg.end),]$CN_state <- state.loh
      } else { print("Check matrix type input") }
  } 
  return(chr.matrix)
}

################ test


res <- make_chrmatrix(chr.probes, chr.segments, matrix.type="major", sample.ploidy)
str(res)

res2 <- make_sampmatrix(probe.df, sample.data, sample.ID, matrix.type="major", sample.ploidy)

str(res2) # cnstate is missing whyyyy?
unique(res2$Chr)
str(res2)
View(sample.data)
