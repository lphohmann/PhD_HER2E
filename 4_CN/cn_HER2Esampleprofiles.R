# Script: Plot total cn profile for each HER2E sample and calc. % of genome altered in HER2E samples

# TODO:
# - 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run (for plot naming)
cohort <- "SCANB" 

# set/create output directory for plots
output.path <- "output/plots/4_CN/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/4_CN/processed/",sep="")
dir.create(data.path)

# packages
source("scripts/4_CN/src/cn_functions.R")
library(ggplot2)
library(tidyverse)
library(readxl)

################################################################################
# loading data
################################################################################

total.cn.scanb <- as.data.frame((read.table(file = 'data/SCANB/4_CN/raw/total_matrix.tsv', sep = '\t', header = TRUE))) 

################################################################################
# plot the total CN for each sample
################################################################################

results <- processing(total.cn.scanb)
total.cn.scanb <- results[1] %>% as.data.frame()
chr.lengths <- results[2] %>% as.data.frame()

save(total.cn.scanb, file = paste(data.path,"CN_total",sep=""))
#load(file = paste(data.path,"CN_total",sep=""))

pdf(file = paste(output.path,cohort,"_sampleprofiles.pdf", sep =""),onefile = TRUE, height = 10.5, width = 14.8)

pb = txtProgressBar(min = 0, max = ncol(total.cn.scanb)-4, initial = 0, style = 3) 
# for each sample
for(i in c(1:(ncol(total.cn.scanb)-4))) { 
    setTxtProgressBar(pb,i)
    # loop iteration example
    sampleID <- colnames(total.cn.scanb)[i+4] #gsub("^.{0,3}", "", colnames(total.cn.scanb)[i+4])
    #print(sampleID)
    sample.data <- as.data.frame(total.cn.scanb[,c(1:4,i+4)]) 
    #print(head(sample.data))
    plot <- ggplot(sample.data, aes(x=sample.data$genome_pos, y=sample.data[,sampleID])) + #, group=1 #scale_y_continuous(limits=c(-2, 2), breaks=c(-2, -1, 0, 1, 2)) +
        geom_step() + #geom_path()
        geom_point() +
        ggtitle(paste("Copy number profile: ",gsub("^.{0,3}", "", sampleID),sep="")) +
        geom_vline(xintercept = chr.lengths$genome, linetype="dotted") +
        ylab("copy number state") +
        xlab("genome position (chromosome)") +
        scale_x_continuous(breaks=chr.lengths$genome,
                           labels=chr.lengths$Chr)+
        theme(plot.title = element_text(size = 30),
              axis.text.x = element_text(size = 20),
              axis.title.x = element_text(size = 25),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 25),
              legend.position = "none")
    print(plot)
}
dev.off()



##########################################################################
# genome altered (%)
##########################################################################

str(gainloss.cn.scanb)

cn.data <- gainloss.cn.scanb %>% select(-c("ProbeID","Chr","Position"))

# alterations per sample 
res <- colSums(cn.data != 0) %>% as.data.frame() %>% rownames_to_column(var = "sampleID") %>% dplyr::rename(Nalt = 2)

# as percent
res$Nalt <- (res$Nalt/nrow(cn.data))*100

boxplot(res$Nalt)
median(res$Nalt) # 43.05 % of genome is altered
