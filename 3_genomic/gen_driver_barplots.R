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
library(gridExtra)

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
#View(all.dmut)

#######################################################################
# plot genes of interest: mutation frequencies
#######################################################################
# PAM50, mutation type, count - but this would count multiple mutations per sample
count.type <- function(data,gene) {
  data[data$Gene==gene,] %>% 
    count(PAM50, Mutation_Type) %>% 
    mutate(Freq=case_when(PAM50=="Her2e" ~ n/)) # check
}

# 1 mut per sample
count.sample <- function(data,gene) {
  data <- data %>% 
    distinct(Sample,Gene, .keep_all = TRUE)
  data[data$Gene==gene,] %>% 
    count(PAM50)
}

gene.vec <-c("ERBB2","ESR1","TP53","PIK3CA","FGFR4")
g <- gene.vec[3]
count.type(all.dmut, gene=g)
count.sample(all.dmut, gene=g)

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

# save plots
pdf(file = plot.file, onefile = TRUE) #, height = 10, width = 15)


# number of plots
n <- length(plot.list)
# plots per page
pp <- 6 
  
for(i in seq(ceiling(n/pp))) { # pp plots per page
  
  # for the last page
  if(n > pp) { 
    k <- pp
    } else {k <- n} # print only the remaining plots
  n <- n - pp
  
  # subset the plot list
  plots.page.list <- plot.list[1:k]
  # remove these plots from the original list
  plot.list[!(plot.list %in% plots.page.list)]
  # 
  print(do.call(gridExtra::grid.arrange, plots.page.list))
}

dev.off()