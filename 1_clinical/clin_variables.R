# Script: Create table with clinicopathological variables comparisons between subtypes

# if run for SCANB, output: 
# - table whole rel 4 + h2low +% review
# - table review 1 

# for METABRIC:
# - table + IC10

# TODO: - create two output tables, 1 for whole rel4 (+ % review and HER2low) and 1 for review 1 cases
# - for metabric include ic10 class distribution
# make table 1 for scanb and mb in one 

# - check that percentage is logical and pvalue
# - exclude cases that dont have all variables?

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # Metabric or SCANB 

# set/create output directory for plots
output.path <- "output/plots/1_clinical/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/1_clinical/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/1_clinical/src/clin_functions.R")
library(tidyverse)
library(matrixStats)
library(DescTools)
library(openxlsx)

#######################################################################
# Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for metabric cohort
if (cohort=="Metabric") {
    
    # load data
    # extract relevant variables
    anno <- loadRData("./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData") %>% 
        mutate(HER2_Low = case_when(
          HER2_IHC_status==1 ~ 1,
          HER2_IHC_status==2 & HER2_SNP6_state!="GAIN" ~ 1)) %>% 
        mutate(LN = ifelse(lymph_nodes_positive > 0, 1, 0)) %>%  
        dplyr::rename(sampleID=METABRIC_ID,
                      NHG=Grade,
                      Size=TumSize,
                      IC10=IC10_Rueda) %>% 
        dplyr::select(sampleID, PAM50, Age, Size, NHG, LN, HER2_Low, IC10) # no pr status available
    
    # replace na values with 0
    anno$HER2_Low[is.na(anno$HER2_Low)] <- 0
    
    # 
    anno$IC10 <- as.factor(anno$IC10)
    
# for SCANB cohort
} else if (cohort=="SCANB") {
    
    # select relevant variables
    anno <- loadRData("./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>% 
        mutate(LN_Bosch = case_when(LNstatus_Bosch == "N0" ~ 0, LNstatus_Bosch == "N+" ~ 1)) %>% 
      dplyr::select(
        GEX.assay, 
        # table 1: N, Age, size, PR, NHG, LN, % review, H2low
        NCN.PAM50, Age, PR, Size.mm, NHG, LN, HER2_Low,
        # table 2: N, Age, size, PR, NHG, LN, H2low
        Bosch_RS1, TumSize_Bosch, NHG_Bosch, LN_Bosch
        ) %>% 
      dplyr::rename(
            sampleID=GEX.assay,
            PAM50=NCN.PAM50,
            Size=Size.mm,
            Size_Bosch = TumSize_Bosch)
    
    # do data type correction for scanb exclusive vars 
    #anno$PR <- as.factor(anno$PR)
    anno$Size_Bosch <- as.numeric(anno$Size_Bosch)
    anno$NHG_Bosch <- as.factor(anno$NHG_Bosch)
    anno$LN_Bosch <- as.factor(anno$LN_Bosch)
    anno$PR[anno$PR == ""] <- NA
}

# get correct data types
anno$PAM50 <- as.factor(anno$PAM50)
anno$NHG[anno$NHG == ""] <- NA
anno$NHG <- as.factor(anno$NHG)
anno$Size <- as.numeric(anno$Size)
anno$Age <- as.numeric(anno$Age)
anno$LN <- as.factor(anno$LN)
anno$HER2_Low <- as.factor(anno$HER2_Low)

# check if i can do that differently
# sample names
Her2.ids <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA.ids <- anno %>% filter(PAM50=="LumA") %>% pull(sampleID)
LumB.ids <- anno %>% filter(PAM50=="LumB") %>% pull(sampleID)
anno <- anno %>% column_to_rownames(var = "sampleID")

#storing format: excel document with sheets for each variable
# keep column names consitent so later I just need to copy & paste tables together
col.names <- c("Variable","HER2E(ref)","HER2E.%","LUMA","LUMA.%","LUMA.pval","LUMB","LUMB.%","LUMB.pval")

########################################################################
# 1. data vars (not review)
# table 1: N, Age, size, PR, NHG, LN, % review, H2low

#######################################################################
# PAM50 count
#######################################################################

# matrix
pam50.count <- as.data.frame(matrix(ncol=length(col.names)))
names(pam50.count) <- col.names
pam50.count$Variable <- "PAM50"

t <- table(anno$PAM50)

# add count data
pam50.count$`HER2E(ref)`[1] <- t[1]
pam50.count$LUMA[1] <- t[2]
pam50.count$LUMB[1] <- t[3]

# add % data
pam50.count$`HER2E.%`[1] <- (t[1]/length(anno$PAM50))*100
pam50.count$`LUMA.%`[1] <- (t[2]/length(anno$PAM50))*100
pam50.count$`LUMB.%`[1] <- (t[3]/length(anno$PAM50))*100

pam50.count

#######################################################################
# Age comparison
#######################################################################

age.data <- anno %>% filter(!is.na(Age))
table(anno$PAM50,is.na(anno$Age))

# matrix
age.df <- as.data.frame(matrix(ncol=length(col.names)))
names(age.df) <- col.names
age.df$Variable <- "Age"

# for Her2 vs LumA
# store result
age.df$`HER2E(ref)` <-  mean(age.data[which(age.data$PAM50=="Her2"),]$Age)
age.df$LUMA <-          mean(age.data[which(age.data$PAM50=="LumA"),]$Age)
age.df$LUMA.pval <-     TwoGroups.ttest(age.data[Her2.ids,"Age"],age.data[LumA.ids,"Age"])$p.value

# store result
age.df$LUMB <-          mean(age.data[which(age.data$PAM50=="LumB"),]$Age)
age.df$LUMB.pval <-     TwoGroups.ttest(age.data[Her2.ids,"Age"],age.data[LumB.ids,"Age"])$p.value

age.df

#######################################################################
# 4. Size
#######################################################################

size.data <- anno %>% filter(!is.na(Size))
table(anno$PAM50,is.na(anno$Size))

# matrix
size.df <- as.data.frame(matrix(ncol=length(col.names)))
names(size.df) <- col.names
size.df$Variable <- "Size"

# for Her2 vs LumA
# store result
size.df$`HER2E(ref)` <-  mean(
  size.data[which(size.data$PAM50=="Her2" & !is.na(size.data$Size)),]$Size)
size.df$LUMA <-          mean(size.data[which(size.data$PAM50=="LumA" & !is.na(size.data$Size)),]$Size)
size.df$LUMA.pval <-     TwoGroups.ttest(size.data[Her2.ids,"Size"],size.data[LumA.ids,"Size"])$p.value

# for Her2 vs LumB
# store result
size.df$LUMB <-          mean(size.data[which(size.data$PAM50=="LumB" & !is.na(size.data$Size)),]$Size)
size.df$LUMB.pval <-     TwoGroups.ttest(size.data[Her2.ids,"Size"],size.data[LumB.ids,"Size"])$p.value

size.df

#######################################################################
# 4. PR status (2x2 contingency table)
#######################################################################

pr.data <- anno %>% filter(!is.na(PR))

# matrix
pr.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
names(pr.df) <- col.names
pr.df$Variable <- c("Negative","Positive")

# for Her2 vs LumA
luma.ct <- table(pr.data$PAM50,pr.data$PR)[c("Her2","LumA"),]
pr.df$LUMA.pval <- fisher.test(luma.ct)$p.value

# for Her2 vs LumB
lumb.ct <- table(pr.data$PAM50,pr.data$PR)[c("Her2","LumB"),]
pr.df$LUMB.pval <- fisher.test(lumb.ct)$p.value

groups <- c("Her2","LumA","LumB")
cols <- c("HER2E","LUMA","LUMB")
for (i in 1:3) {
  type <- groups[i]
  col <- cols[i]
  type.dat <- pr.data[which(pr.data$PAM50==type),]
  type.dat$PAM50 <- droplevels(type.dat$PAM50)
  type.counts <- table(type.dat$PAM50,type.dat$PR)
  
  # count column
  if (col=="HER2E") {
    pr.df[[paste(col,"(ref)",sep="")]][1] <- type.counts[1]
    pr.df[[paste(col,"(ref)",sep="")]][2] <- type.counts[2]
  } else {
    pr.df[[col]][1] <- type.counts[1]
    pr.df[[col]][2] <- type.counts[2]
  }
  # % column
  pr.df[[paste(col,".%",sep="")]][1] <- round(type.counts[1]/sum(type.counts)*100)
  pr.df[[paste(col,".%",sep="")]][2] <- round(type.counts[2]/sum(type.counts)*100)
} 

pr.df

#######################################################################
# 4. NHG (2x3 contingency table)
#######################################################################

nhg.data <- anno %>% filter(!is.na(NHG))
table(anno$PAM50,is.na(anno$NHG))

# matrix
grade.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 3))
names(grade.df) <- col.names
grade.df$Variable <- c("NHG1","NHG2","NHG3")

# for Her2 vs LumA
luma.ct <- table(nhg.data$PAM50,nhg.data$NHG)[c("Her2","LumA"),]
grade.df$LUMA.pval <- fisher.test(luma.ct)$p.value

#GTest(table(a_nhg$PAM50,a_nhg$NHG))

# for Her2 vs LumB
lumb.ct <- table(nhg.data$PAM50,nhg.data$NHG)[c("Her2","LumB"),]
grade.df$LUMB.pval <- fisher.test(lumb.ct)$p.value

groups <- c("Her2","LumA","LumB")
cols <- c("HER2E","LUMA","LUMB")
for (i in 1:3) {
  type <- groups[i]
  col <- cols[i]
  type.dat <- nhg.data[which(nhg.data$PAM50==type),]
  type.dat$PAM50 <- droplevels(type.dat$PAM50)
  type.counts <- table(type.dat$PAM50,type.dat$NHG)

  # count column
  if (col=="HER2E") {
    grade.df[[paste(col,"(ref)",sep="")]][1] <- type.counts[1]
    grade.df[[paste(col,"(ref)",sep="")]][2] <- type.counts[2]
    grade.df[[paste(col,"(ref)",sep="")]][3] <- type.counts[3]
  } else {
  grade.df[[col]][1] <- type.counts[1]
  grade.df[[col]][2] <- type.counts[2]
  grade.df[[col]][3] <- type.counts[3]
  }
  # % column
  grade.df[[paste(col,".%",sep="")]][1] <- round(type.counts[1]/sum(type.counts)*100)
  grade.df[[paste(col,".%",sep="")]][2] <- round(type.counts[2]/sum(type.counts)*100)
  grade.df[[paste(col,".%",sep="")]][3] <- round(type.counts[3]/sum(type.counts)*100)
} 

#######################################################################
# 4. LN (2x2 contingency table)
#######################################################################

ln.data <- anno %>% filter(!is.na(LN))
table(anno$PAM50,is.na(anno$LN))

# matrix
ln.df <- as.data.frame(matrix(ncol=length(col.names),nrow = 2))
names(ln.df) <- col.names
ln.df$Variable <- c("N0","N+")

# for Her2 vs LumA
luma.ct <- table(ln.data$PAM50,ln.data$LN)[c("Her2","LumA"),]
ln.df$LUMA.pval <- fisher.test(luma.ct)$p.value

# for Her2 vs LumB
lumb.ct <- table(ln.data$PAM50,ln.data$LN)[c("Her2","LumB"),]
ln.df$LUMB.pval <- fisher.test(lumb.ct)$p.value

groups <- c("Her2","LumA","LumB")
cols <- c("HER2E","LUMA","LUMB")
for (i in 1:3) {
  type <- groups[i]
  col <- cols[i]
  type.dat <- ln.data[which(ln.data$PAM50==type),]
  type.dat$PAM50 <- droplevels(type.dat$PAM50)
  type.counts <- table(type.dat$PAM50,type.dat$LN)
  
  # count column
  if (col=="HER2E") {
    ln.df[[paste(col,"(ref)",sep="")]][1] <- type.counts[1]
    ln.df[[paste(col,"(ref)",sep="")]][2] <- type.counts[2]
  } else {
    ln.df[[col]][1] <- type.counts[1]
    ln.df[[col]][2] <- type.counts[2]
  }
  # % column
  ln.df[[paste(col,".%",sep="")]][1] <- round(type.counts[1]/sum(type.counts)*100)
  ln.df[[paste(col,".%",sep="")]][2] <- round(type.counts[2]/sum(type.counts)*100)
} 

ln.df

#######################################################################
# export to excel
#######################################################################

# save
write.xlsx(final_matrix, file = paste(data.path,"clin_variables.xlsx",sep=""),overwrite = TRUE)

########################################################################
# 2. review data vars

#######################################################################
# HER2low frequency (2x2 contingency table)
#######################################################################

# for Her2 vs LumA
a_hl <- subset(anno,PAM50 %in% c("Her2","LumA")) %>% droplevels()
a_cont_table <- table(a_hl$PAM50,a_hl$HER2_Low)
a_result <- fisher.test(a_cont_table)

# for Her2 vs LumB
b_hl <- subset(anno,PAM50 %in% c("Her2","LumB")) %>% droplevels()
b_cont_table <- table(b_hl$PAM50,b_hl$HER2_Low)
b_result <- fisher.test(b_cont_table) 

H_tot <- a_cont_table[1] + a_cont_table[3]
A_tot <- a_cont_table[2] + a_cont_table[4]
B_tot <- b_cont_table[2] + b_cont_table[4]

#her2e
final_matrix$`HER2E(ref)`[2] <- a_cont_table[3]
final_matrix$`HER2E.%`[2] <- round(a_cont_table[3]/H_tot*100)

#luma
final_matrix$LUMA[2] <- a_cont_table[4]
final_matrix$`LUMA.%`[2] <- round(a_cont_table[4]/A_tot*100)
final_matrix$LUMA.pval[2] <- a_result$p.value

#lumb
final_matrix$LUMB[2] <- b_cont_table[4]
final_matrix$`LUMB.%`[2] <- round(b_cont_table[4]/B_tot*100)
final_matrix$LUMB.pval[2] <- b_result$p.value
    
final_matrix



