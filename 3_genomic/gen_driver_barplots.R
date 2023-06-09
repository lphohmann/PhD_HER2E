# Script: Driver mutations barplots for SCANB and BASIS

#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "COMBINED" 

# set/create output directory for plots
output.path <- "output/plots/3_genomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/3_genomic/processed/",sep="")
dir.create(data.path)

# plot
plot.list <- list() # object to store plots; note: if the output is not in string format use capture.output()
plot.file <- paste(output.path,cohort,"_HER2n_driverBarplots.pdf",sep = "")

#packages
source("scripts/3_genomic/src/gen_functions.R")
source("scripts/4_CN/src/cn_functions.R")
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

# scanb mut data
scanb.dmut <- loadRData("data/SCANB/3_genomic/processed/driver_mutations_all.RData") %>% 
  mutate(PAM50 = "Her2e") %>% 
  dplyr::rename(Sample=sample, Gene=gene, Effect=variant_class) %>% 
  # split the mutation_type column in two columns (add effect column)
  mutate(Mutation_Type = sub("_$","",str_extract(Effect,"^.*?_"))) %>% 
  mutate(Effect = gsub("^.*?_","",Effect)) %>% 
  mutate(Mutation_Type = if_else(Effect=="amplified","CN",Mutation_Type)) 

# basis mut data
basis.dmut <- as.data.frame(read_excel("data/BASIS/3_genomic/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx", sheet = "COMBINED_EVENTS")) %>% 
  filter(Mutation_Type %in% c("Sub","Ins","Del","CopyNumber")) %>% 
  # split into two rows where #Chr8:(ZNF703/FGFR1)
  mutate(Gene = gsub("^Chr8:","",Gene)) %>% #Chr8:(ZNF703/FGFR1)
  mutate(Gene = gsub("[()]", "", Gene)) %>%
  mutate(Gene = strsplit(Gene, "/")) %>% 
  unnest(Gene) %>% 
  mutate(Mutation_Type = case_when(Mutation_Type=="Sub" ~ "sub",
                                    Mutation_Type=="Ins" ~ "indel",
                                    Mutation_Type=="Del" ~ "indel",
                                    Mutation_Type=="CopyNumber" ~ "CN")) %>% 
  mutate(Effect = if_else(Effect=="amp","amplified",Effect)) %>% 
  dplyr::select(-c(Tumour_name))

# get pam50 annotation
basis.anno <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData") %>% 
  filter(ClinicalGroup == "ERposHER2neg" & PAM50_AIMS %in% c("LumA","LumB")) %>% 
  dplyr::select(c("sample_name","PAM50_AIMS")) %>% 
  dplyr::rename(Sample=sample_name,PAM50=PAM50_AIMS)


basis.dmut <- merge(basis.anno,basis.dmut,by="Sample")

# combined
all.dmut <- rbind(basis.dmut,scanb.dmut)

#######################################################################
# plot genes of interest
#######################################################################
# PAM50, mutation type, count - but this would count multiple mutations per sample
count.type <- function(data,gene) {
  data[data$Gene==gene,] %>% 
    count(PAM50, Mutation_Type)
}

# 1 mut per sample
count.sample <- function(data,gene) {
  data <- data %>% 
    distinct(Sample,Gene, .keep_all = TRUE)
  data[data$Gene==gene,] %>% 
    count(PAM50)
}

gene.vec <-c("ERBB2","ESR1","TP53","PIK3CA","FGFR4")

for (g in gene.vec) {
  
  # Plots overall mutations regardless of multiple per sample
  p1 <- ggplot(count.type(all.dmut, gene=g), aes(fill=Mutation_Type,x=PAM50,y=n))+
    geom_bar(position="stack", stat="identity")
  
  # Plots freq samples mutated in PAM50 groups
  p2 <- ggplot(count.sample(all.dmut, gene=g), aes(x=PAM50,y=n))+
    geom_bar(position="stack", stat="identity")
  
  plot.list <- append(plot.list,list(p1)) 
  plot.list <- append(plot.list,list(p2)) 
  
}

# how to save multiple plots per page, also check firs tthat the loop above gives the right plots

# save plots
pdf(file = plot.file, onefile = TRUE) #, height = 10, width = 15)

# flatten plot list
plot.list <- do.call(c, plot.list)

for (i in 1:length(plot.list)) {
  #print(plot.list[[i]])
  do.call("grid.arrange", print(plot.list[[i]]))
}

dev.off()

library(gridExtra)

pdf("plots.pdf", onefile = TRUE)
for (i in seq(length(p))) {
  do.call("grid.arrange", p[[i]])  
}
dev.off()