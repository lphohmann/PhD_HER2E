# Script: Create table with clinicopathological variables comparisons between subtypes
# cohorts metabric or scanb review 1

# TODO: - create final matrix directly
# - check that percentage is logical and pvalue

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
    load("data/METABRIC/1_clinical/raw/Merged_annotations.RData")
    # extract relevant variables
    anno <- anno %>% 
        mutate(HER2_Low = case_when(HER2_IHC_status==1 ~ 1,
                                    HER2_IHC_status==2 & HER2_SNP6_state!="GAIN" ~ 1)) %>% 
        mutate(LN = ifelse(lymph_nodes_positive > 0, 1, 0)) %>%  
        filter(grepl('ERpHER2n', ClinGroup)) %>% 
        dplyr::rename(sampleID=METABRIC_ID,
                      NHG=Grade,
                      Size=TumSize) %>% 
        dplyr::select(sampleID,PAM50,Size,Age,NHG,LN,HER2_Low)
    
    # replace na values with 0
    anno$HER2_Low[is.na(anno$HER2_Low)] <- 0
    
# for SCANB cohort
} else if (cohort=="SCANB") {
    
    # load data
    load("data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_with_ExternalReview_Bosch_data.RData")
    # extract relevant variables
    anno <- pam50.frame %>% 
        filter(fuV8 == 1) %>% 
        filter(!is.na(treatment_Bosch)) %>% # only include the review data
        filter(ERpHER2n_Bosch==1) %>%
        mutate(LN = case_when(LNstatus_Bosch == "N0" ~ 0, LNstatus_Bosch == "N+" ~ 1)) %>% 
        dplyr::rename(
            sampleID=rba_rel4,
            PAM50=PAM50_NCN_rel4,
            NHG=NHG_Bosch,
            Size=TumSize_Bosch) %>%
        dplyr::select(
            sampleID, PAM50, Size, Age, NHG, LN, HER2_Low)
}

# get correct data types
anno$PAM50 <- as.factor(anno$PAM50)
anno$NHG[anno$NHG == ""] <- NA
anno$NHG <- as.factor(anno$NHG)
anno$Size <- as.numeric(anno$Size)
anno$Age <- as.numeric(anno$Age)
anno$LN <- as.factor(anno$LN)
anno$HER2_Low <- as.factor(anno$HER2_Low)


# filter to only include subjects that are PAM50 == Her2 | LumA | LumB
anno <- anno %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) # basal
anno$PAM50 <- droplevels(anno$PAM50) # drop empty levels
barplot(table(anno$PAM50))

# check if i can do that differently
# sample names
Her2_samp <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA_samp <- anno %>% filter(PAM50=="LumA") %>% pull(sampleID)
LumB_samp <- anno %>% filter(PAM50=="LumB") %>% pull(sampleID)
anno <- anno %>% column_to_rownames(var = "sampleID")

#######################################################################
# Age
#######################################################################

# initialize storing object
age_matrix <- data.frame(matrix(ncol = 0, nrow = 1))

# for Her2 vs LumA
# equal variance check
if (var.test(unlist(anno[Her2_samp,"Age"]),unlist(anno[LumA_samp,"Age"]), 
             alternative = "two.sided")$p.value <= 0.05) {
    result <- t.test(anno[Her2_samp,"Age"],anno[LumA_samp,"Age"], 
           var.equal = FALSE)
} else {
    result <- t.test(anno[Her2_samp,"Age"],anno[LumA_samp,"Age"], 
           var.equal = TRUE)  
}

# store result
A_pvalue <- result$p.value
H_mean <- result$estimate[1]
A_mean <- result$estimate[2]

# for Her2 vs LumB
# equal variance check
if (var.test(unlist(anno[Her2_samp,"Age"]),unlist(anno[LumB_samp,"Age"]), 
             alternative = "two.sided")$p.value <= 0.05) {
    result <- t.test(anno[Her2_samp,"Age"],anno[LumB_samp,"Age"], 
                     var.equal = FALSE)
} else {
    result <- t.test(anno[Her2_samp,"Age"],anno[LumB_samp,"Age"], 
                     var.equal = TRUE)  
}

# store result
B_pvalue <- result$p.value
B_mean <- result$estimate[2]

# add to final result df
age_matrix$H_mean <- H_mean
age_matrix$A_mean <- A_mean
age_matrix$A_pvalue <- A_pvalue
age_matrix$B_mean <- B_mean
age_matrix$B_pvalue <- B_pvalue

row.names(age_matrix) <- c("Age")

#View(age_matrix)
#######################################################################
# 4. Size
#######################################################################

# initialize storing object
size_matrix <- data.frame(matrix(ncol = 0, nrow = 1))

# for Her2 vs LumA
# equal variance check
if (var.test(unlist(anno[Her2_samp,"Size"]),unlist(anno[LumA_samp,"Size"]), 
             alternative = "two.sided")$p.value <= 0.05) {
    result <- t.test(anno[Her2_samp,"Size"],anno[LumA_samp,"Size"], 
                     var.equal = FALSE)
} else {
    result <- t.test(anno[Her2_samp,"Size"],anno[LumA_samp,"Size"], 
                     var.equal = TRUE)  
}

# store result
A_pvalue <- result$p.value
H_mean <- result$estimate[1]
A_mean <- result$estimate[2]

# for Her2 vs LumB
# equal variance check
if (var.test(unlist(anno[Her2_samp,"Size"]),unlist(anno[LumB_samp,"Size"]), 
             alternative = "two.sided")$p.value <= 0.05) {
    result <- t.test(anno[Her2_samp,"Size"],anno[LumB_samp,"Size"], 
                     var.equal = FALSE)
} else {
    result <- t.test(anno[Her2_samp,"Size"],anno[LumB_samp,"Size"], 
                     var.equal = TRUE)  
}

# store result
B_pvalue <- result$p.value
B_mean <- result$estimate[2]

# add to final result df
size_matrix$H_mean <- H_mean
size_matrix$A_mean <- A_mean
size_matrix$A_pvalue <- A_pvalue
size_matrix$B_mean <- B_mean
size_matrix$B_pvalue <- B_pvalue

row.names(size_matrix) <- c("Size")
#View(size_matrix)

#######################################################################
# 4. NHG (2x3 contingency table)
#######################################################################

nhg_matrix <- data.frame(matrix(ncol = 0, nrow = 1))

# for Her2 vs LumA
a_nhg <- subset(anno,PAM50 %in% c("Her2","LumA")) %>% droplevels()
a_cont_table <- table(a_nhg$PAM50,a_nhg$NHG)
#GTest(table(a_nhg$PAM50,a_nhg$NHG))
a_result <- fisher.test(a_cont_table) # maybe use fishers due to expected n below 5

        
# for Her2 vs LumB
b_nhg <- subset(anno,PAM50 %in% c("Her2","LumB")) %>% droplevels()
b_cont_table <- table(b_nhg$PAM50,b_nhg$NHG)

b_result <- fisher.test(b_cont_table) # maybe use fishers due to expected n below 5

H_tot <- a_cont_table[1] + a_cont_table[3] + a_cont_table[5]
A_tot <- a_cont_table[2] + a_cont_table[4] + a_cont_table[6]
B_tot <- b_cont_table[2] + b_cont_table[4] + b_cont_table[6]

# final results
nhg_matrix$A_nhg_pvalue <- a_result$p.value
nhg_matrix$B_nhg_pvalue <- b_result$p.value

nhg_matrix$H_nhg1_count <- a_cont_table[1]
nhg_matrix$A_nhg1_count <- a_cont_table[2]
nhg_matrix$B_nhg1_count <- b_cont_table[2]

nhg_matrix$H_nhg1_perc <- round(a_cont_table[1]/H_tot*100)
nhg_matrix$A_nhg1_perc <- round(a_cont_table[2]/A_tot*100)
nhg_matrix$B_nhg1_perc <- round(b_cont_table[2]/B_tot*100)

nhg_matrix$H_nhg2_count <- a_cont_table[3]
nhg_matrix$A_nhg2_count <- a_cont_table[4]
nhg_matrix$B_nhg2_count <- b_cont_table[4]

nhg_matrix$H_nhg2_perc <- round(a_cont_table[3]/H_tot*100)
nhg_matrix$A_nhg2_perc <- round(a_cont_table[4]/A_tot*100)
nhg_matrix$B_nhg2_perc <- round(b_cont_table[4]/B_tot*100)

nhg_matrix$H_nhg3_count <- a_cont_table[5]
nhg_matrix$A_nhg3_count <- a_cont_table[6]
nhg_matrix$B_nhg3_count <- b_cont_table[6]

nhg_matrix$H_nhg3_perc <- round(a_cont_table[5]/H_tot*100)
nhg_matrix$A_nhg3_perc <- round(a_cont_table[6]/A_tot*100)
nhg_matrix$B_nhg3_perc <- round(b_cont_table[6]/B_tot*100)

#View(nhg_matrix)

#######################################################################
# 4. LN (2x2 contingency table)
#######################################################################

ln_matrix <- data.frame(matrix(ncol = 0, nrow = 1))

# for Her2 vs LumA
a_ln <- subset(anno,PAM50 %in% c("Her2","LumA")) %>% droplevels()
a_cont_table <- table(a_ln$PAM50,a_ln$LN)
a_result <- fisher.test(a_cont_table)

# for Her2 vs LumB
b_ln <- subset(anno,PAM50 %in% c("Her2","LumB")) %>% droplevels()
b_cont_table <- table(b_ln$PAM50,b_ln$LN)
b_result <- fisher.test(b_cont_table) 

H_tot <- a_cont_table[1] + a_cont_table[3]
A_tot <- a_cont_table[2] + a_cont_table[4]
B_tot <- b_cont_table[2] + b_cont_table[4]

# add to final matrix

#lumb
final_matrix$LUMB[9] <- b_cont_table[2]
final_matrix$LUMB[10] <- b_cont_table[4]
final_matrix$`LUMB.%`[9] <- round(b_cont_table[2]/B_tot*100)
final_matrix$`LUMB.%`[10] <-round(b_cont_table[4]/B_tot*100)
final_matrix$LUMB.pval[8] <- b_result$p.value
    
#luma
final_matrix$LUMA[9] <- a_cont_table[2]
final_matrix$LUMA[10] <- a_cont_table[4]
final_matrix$`LUMA.%`[9] <- round(a_cont_table[1]/H_tot*100)
final_matrix$`LUMA.%`[10] <- round(a_cont_table[4]/A_tot*100)
final_matrix$LUMA.pval[8] <- a_result$p.value

#her2e
final_matrix$`HER2E(ref)`[9] <- a_cont_table[1]
final_matrix$`HER2E(ref)`[10] <- a_cont_table[3]
final_matrix$`HER2E.%`[9] <- round(a_cont_table[1]/H_tot*100)
final_matrix$`HER2E.%`[10] <- round(a_cont_table[3]/H_tot*100)

final_matrix    # da stimmt was nicht
#View(ln_matrix)

#######################################################################
# 4. export to excel
#######################################################################

# get desired format
final_matrix <- data.frame(matrix(ncol = 9, nrow = 10))
colnames(final_matrix) <- c("Variable","HER2E(ref)","HER2E.%","LUMA","LUMA.%","LUMA.pval","LUMB","LUMB.%","LUMB.pval")
final_matrix$Variable <- c("N","Age in years (mean)","Tumor Size in mm (mean)","Tumor grade (count):","grade 1","grade 2", "grade 3","Lymph node status (count):","N0","N+")
final_matrix
# add count data

final_matrix$`HER2E(ref)`[1] <- table(anno$PAM50)[1]
final_matrix$LUMA[1] <- table(anno$PAM50)[2]
final_matrix$LUMB[1] <- table(anno$PAM50)[3]
final_matrix    



# save
datasheets <- list("Final" = final_matrix, "Age" = age_matrix, "Size" = size_matrix, "NHG" = nhg_matrix, "LN" = ln_matrix)
write.xlsx(datasheets, file = paste(data.path,"clin_variables.xlsx",sep=""),overwrite = TRUE)

##################################

# plot exploration
ggplot(anno, aes(x=as.factor(PAM50),y=Age)) +
    geom_boxplot(fill="slateblue",alpha=0.2) +
    xlab("PAM50 subtype") +
    ylab("Age") +
    ggtitle("one of many") +
    geom_text(data=as.data.frame(dplyr::count(anno, PAM50)), aes(y = 0, label = paste("n=",n,sep = "")))


###
#######################################################################
# 4. HER2low frequency (2x2 contingency table)
#######################################################################

hl_matrix <- data.frame(matrix(ncol = 0, nrow = 1))


# for Her2 vs LumA
a_hl <- subset(anno,PAM50 %in% c("Her2","LumA")) %>% droplevels()
a_cont_table <- table(a_hl$PAM50,a_hl$HER2_Low)
a_cont_table
a_result <- fisher.test(a_cont_table)

# for Her2 vs LumB
b_hl <- subset(anno,PAM50 %in% c("Her2","LumB")) %>% droplevels()
b_cont_table <- table(b_hl$PAM50,b_hl$HER2_Low)
b_cont_table
b_result <- fisher.test(b_cont_table) 

H_tot <- a_cont_table[1] + a_cont_table[3]
A_tot <- a_cont_table[2] + a_cont_table[4]
B_tot <- b_cont_table[2] + b_cont_table[4]

# final results
hl_matrix$A_hl_pvalue <- a_result$p.value
hl_matrix$B_hl_pvalue <- b_result$p.value

hl_matrix$H_H0_count <- a_cont_table[1]
hl_matrix$A_H0_count <- a_cont_table[2]
hl_matrix$B_H0_count <- b_cont_table[2]

hl_matrix$H_H0_perc <- round(a_cont_table[1]/H_tot*100)
hl_matrix$A_H0_perc <- round(a_cont_table[2]/A_tot*100)
hl_matrix$B_H0_perc <- round(b_cont_table[2]/B_tot*100)

hl_matrix$H_Hlow_count <- a_cont_table[3]
hl_matrix$A_Hlow_count <- a_cont_table[4]
hl_matrix$B_Hlow_count <- b_cont_table[4]

hl_matrix$H_Hlow_perc <- round(a_cont_table[3]/H_tot*100)
hl_matrix$A_Hlow_perc <- round(a_cont_table[4]/A_tot*100)
hl_matrix$B_Hlow_perc <- round(b_cont_table[4]/B_tot*100)

View(hl_matrix)

write.xlsx(hl_matrix, file = paste("~/Desktop/MTP_project/Data/",cohort,"/her2low_matrix.xlsx",sep=''),overwrite = TRUE)


# check treatments
anno <- as.data.frame(pam50.frame) %>% 
    filter(PAM50_NCN_ProSigna_rel4 %in% c("LumA", "LumB", "Her2")) %>% 
    mutate(LN = case_when(LNstatus_Bosch == "N0" ~ 0, LNstatus_Bosch == "N+" ~ 1)) %>%
    dplyr::rename(sampleID=rba_rel4,PAM50=PAM50_NCN_ProSigna_rel4,NHG=NHG_Bosch,Size=TumSize_Bosch) %>% 
    dplyr::select(sampleID,PAM50,Size, Age, NHG, LN, HER2_Low,AI_Bosch,Tamoxifen_Bosch,Chemo_Regime_Bosch,Chemo_Bosch) 


# sample namesx
Her2_samp <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA_samp <- anno %>% filter(PAM50=="LumA") %>% pull(sampleID)
LumB_samp <- anno %>% filter(PAM50=="LumB") %>% pull(sampleID)
anno <- anno %>% column_to_rownames(var = "sampleID")

table(anno[Her2_samp,]$Tamoxifen_Bosch) # 22 yes (59%) 14 no
table(anno[Her2_samp,]$AI_Bosch) # 16 yes (44%), 20 n0
table(anno[Her2_samp,]$AI_Bosch,anno[Her2_samp,]$Tamoxifen_Bosch)

x <- anno %>% rownames_to_column(var="ids")
h <- x %>%
    filter(ids %in% Her2_samp) %>%
    filter(Chemo_Bosch==TRUE)
table(h$AI_Bosch)
table(h$Tamoxifen_Bosch)
