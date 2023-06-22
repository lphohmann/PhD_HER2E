# Script: checking stuff
#TODO: 

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

#packages
source("./scripts/1_clinical/src/clin_functions.R")
library(ggplot2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)

##########################################################################
# create summary objects
##########################################################################

# 1. SCANB ; ERpHER2n ; LUMA ; LUMB ; HER2E
SCANB.ERpHER2n.data <- loadRData(
  "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData") %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(fuV8==TRUE) %>% 
  filter(ER=="Positive") %>% 
  filter(HER2=="Negative") %>% 
  filter(NCN.PAM50 %in% c("LumA","LumB","Her2"))
# add the metagene scores
SCANB.mg.scores <- loadRData(
  "./data/SCANB/2_transcriptomic/processed/mg_anno.RData")[[1]] %>% 
  dplyr::select(-c(PAM50)) %>% 
  dplyr::rename(GEX.assay = "sampleID")
SCANB.ERpHER2n.mg.data <- merge(SCANB.ERpHER2n.data,SCANB.mg.scores,by="GEX.assay",all=TRUE)

save(SCANB.ERpHER2n.mg.data,file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData")

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# 2. SCANB ; ERpHER2p ; ERpHER2nHER2E
SCANB.ERpHER2p.data <- loadRData(
  "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData") %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(fuV8==TRUE) %>% 
  filter(ER=="Positive") %>% 
  mutate(Group = case_when(
    HER2 == "Negative" & NCN.PAM50 == "Her2" ~ "HER2n_HER2E",
    HER2 == "Positive" & NCN.PAM50 == "Her2" ~ "HER2p_HER2E",
    HER2 == "Positive" & NCN.PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
  filter(Group %in% 
           c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E"))
# add metagene scores
SCANB.mg.scores <- loadRData(
  "./data/SCANB/2_transcriptomic/processed/mg_anno_HER2p.RData")[[1]] %>% 
  dplyr::select(-c(Group)) %>% 
  dplyr::rename(GEX.assay = "sampleID")
SCANB.ERpHER2p.mg.data <- merge(SCANB.ERpHER2p.data,SCANB.mg.scores,by="GEX.assay",all=TRUE)

save(SCANB.ERpHER2p.mg.data,file="./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2p.RData")

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# 3. METARBRIC ; ERpHER2n ; LUMA ; LUMB ; HER2E
#just cheking definition of ERpHER2n
#gg <- anno[which(anno$ER_IHC_status=="pos" & anno$HER2_SNP6_state!="GAIN" & anno$HER2_SNP6_state!="UNDEF"),]
#setequal(gg$METABRIC_ID,anno[which(anno$ClinGroup == "ERpHER2n"),]$METABRIC_ID)

METABRIC.ERpHER2n.data <- loadRData(
  "data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>% 
  filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
  filter(grepl('ERpHER2n', ClinGroup))
METABRIC.mg.scores <- loadRData(
  "./data/METABRIC/2_transcriptomic/processed/mg_anno.RData")[[1]] %>% 
  dplyr::select(-c(PAM50)) %>% 
  dplyr::rename(METABRIC_ID = "sampleID")
METABRIC.ERpHER2n.mg.data <- merge(METABRIC.ERpHER2n.data,METABRIC.mg.scores,by="METABRIC_ID",all=TRUE)

save(METABRIC.ERpHER2n.mg.data,file="./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData")

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

# 4. METABRIC ; ERpHER2p ; ERpHER2nHER2E
# extract relevant variables
METABRIC.ERpHER2p.data <- loadRData(
  "data/METABRIC/1_clinical/raw/Merged_annotations.RData") %>% 
  mutate(ER = if_else(
    ER_IHC_status=="pos","Positive","Negative")) %>% 
  mutate(HER2 = case_when(
    HER2_SNP6_state == "GAIN" ~ "Positive",
    HER2_SNP6_state == "UNDEF" ~ "Undefined",
    HER2_SNP6_state %in% c("LOSS","NEUT") ~ "Negative")) %>% 
  filter(ER == "Positive") %>% 
  dplyr::select(-c("group")) %>% 
  mutate(Group = case_when(
    HER2 == "Negative" & PAM50 == "Her2" ~ "HER2n_HER2E",
    HER2 == "Positive" & PAM50 == "Her2" ~ "HER2p_HER2E",
    HER2 == "Positive" & PAM50 != "Her2" ~ "HER2p_nonHER2E")) %>% 
  filter(Group %in% 
           c("HER2n_HER2E","HER2p_HER2E","HER2p_nonHER2E"))
METABRIC.mg.scores <- loadRData(
  "./data/METABRIC/2_transcriptomic/processed/mg_anno_HER2p.RData")[[1]] %>% 
  dplyr::select(-c(Group)) %>% 
  dplyr::rename(METABRIC_ID = "sampleID")
METABRIC.ERpHER2p.mg.data <- merge(METABRIC.ERpHER2p.data,METABRIC.mg.scores,by="METABRIC_ID",all=TRUE)

save(METABRIC.ERpHER2p.mg.data,file="./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2p.RData")


#
load("./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData")
load("./data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2p.RData")
load("./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2n.RData")
load("./data/METABRIC/1_clinical/processed/Merged_annotations_ERpHER2p.RData")

colnames(SCANB.ERpHER2n.mg.data)
colnames(SCANB.ERpHER2p.mg.data)
colnames(METABRIC.ERpHER2n.mg.data)
colnames(METABRIC.ERpHER2p.mg.data)
table(SCANB.ERpHER2n.mg.data$NCN.PAM50)
table(SCANB.ERpHER2p.mg.data$Group)
table(METABRIC.ERpHER2n.mg.data$PAM50)
table(METABRIC.ERpHER2p.mg.data$Group)

##########################################################################
# create clinical summary object SCANB
##########################################################################

test <- loadRData("./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData") %>% filter(Follow.up.cohort==TRUE)
View(head(test))

# make clin summary object
clin.data <- loadRData("./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData") %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(fuV8==TRUE) %>% 
  dplyr::select(c("GEX.assay","Sample","Follow.up.cohort",
                  "ER","PR","HER2","NCN.PAM50","Age","LN","NHG","Size.mm",
                  "TreatGroup","ClinGroup",
                  "OS","OSbin","DRFI","DRFIbin","RFI","RFIbin",
                  "Bosch_RS1","DRFI_bin_Bosch","DRFI_Bosch","RFI_bin_Bosch","RFI_Bosch",
                  "IDFS_bin_Bosch","IDFS_Bosch","ERpHER2n_Bosch",
                  "Chemo_Bosch","ET_Bosch","AI_Bosch","Tamoxifen_Bosch","TumSize_Bosch","NHG_Bosch","LNstatus_Bosch","Ki67_Bosch","Ki67_Bosch_RS2","HER2_Low"))
#View(clin.data)
#save(clin.data,file = "./data/SCANB/1_clinical/processed/SCANB_clinData.RData")

##########################################################################
# matching wgs data to annotation data
##########################################################################

# load annotation data 1
data.ids <- as.data.frame(
    read_excel("data/SCANB/3_genomic/raw/JVCimpression2022-11-09_ERpos_WGS_Batch1-3.xlsx")) %>% 
    dplyr::select(c(SENT.TUMOR, TUMOR.alias)) %>% 
    dplyr::rename(sampleID = SENT.TUMOR) 

# load annotation data 2
data <- as.data.frame(
    read_excel("data/SCANB/3_genomic/raw/SCANB_ERpos_Summary_11Nov22.xlsx",sheet=2)) %>% 
    dplyr::select(c(Batch, Tumour,`Final QC`)) %>% 
    mutate(sampleID = ifelse(
        Batch == 3, data.ids[data.ids$TUMOR.alias==Tumour,]$sampleID,Tumour))

data$sampleID <- substr(data$sampleID,1,7)
test <- loadRData("./data/SCANB/1_clinical/raw/fu.all_v1.1.Rdata")
View(head(test))

# load clin annotation data
clin.rel4 <- as.data.frame(
    read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))

View(clin.rel4)
# check basal n
ERpHER2n.Basal <- clin.rel4 %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(NCN.PAM50 == "Basal") %>% 
  filter(ER=="Positive") %>% 
  filter(HER2=="Negative") %>% 
  dplyr::select(Sample, ER, HER2, NCN.PAM50)
save(ERpHER2n.Basal,file="./data/SCANB/1_clinical/processed/ERpHER2nBasal_samples.RData")
rm(list=ls())
load("./data/SCANB/1_clinical/processed/ERpHER2nBasal_samples.RData")

# select subgroup data
anno <- clin.rel4 %>% 
    filter(Follow.up.cohort==TRUE) %>% 
    filter(NCN.PAM50 == "Her2") %>% 
    filter(ER=="Positive") %>% 
    filter(HER2=="Negative") %>% 
    dplyr::rename(sampleID = Sample, PAM50 = NCN.PAM50) %>% 
    dplyr::select(sampleID, ER, HER2, PAM50)

her2e.wgs.samples.df <- merge(anno, data, by = "sampleID") %>% 
    dplyr::rename(rel4.Sample = sampleID, 
                  rel4.ER = ER,
                  rel4.HER2 = HER2,
                  rel4.PAM50 = PAM50)

save(her2e.wgs.samples.df, file = "./data/SCANB/3_genomic/raw/her2samples_wgs_status.RData")
load(file = "./data/SCANB/3_genomic/raw/her2samples_wgs_status.RData")
View(her2e.wgs.samples.df)




### 
#wgs data subset check
################# are all samples really her2e
# only include samples with a pass in the final QC

# all her2e samples
ERpHER2n.dat <- loadRData("data/SCANB/1_clinical/processed/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data_ERpHER2n.RData") %>% 
  dplyr::select(c(Sample,ER,HER2,NCN.PAM50))
her2e.samples <- ERpHER2n.dat %>%
  filter(NCN.PAM50=="Her2") %>% 
  pull(Sample)
View(ERpHER2n.dat)
# qc annotation
wgs.df <- read_excel("data/SCANB/3_genomic/raw/SCANB_ERpos_Summary_11Nov22.xlsx", sheet = "Batch1+2+3") %>% 
  dplyr::select(Tumour,Batch,`Final QC`) %>% 
  filter(`Final QC`!="MISSING")#"pass")
View(wgs.df)

# get normal IDs before merge
id.key <- read_excel("data/SCANB/3_genomic/raw/JVCimpression2022-11-09_ERpos_WGS_Batch1-3.xlsx") %>% 
  dplyr::select(c(SENT.TUMOR,SENT.TUMOR.aliquot,TUMOR.alias)) %>% 
  mutate(SENT.TUMOR.aliquot = gsub("\\.","_",SENT.TUMOR.aliquot))

# 1. convert epb IDs to normal sample IDs
#userdata$ID <- userids$ID[match(userdata$ID, userids$USER)]
wgs.df$Sample <- id.key$SENT.TUMOR[match(wgs.df$Tumour, id.key$SENT.TUMOR.aliquot)]
# 2. convert the other IDs to normal sample IDs
wgs.df$Sample2 <- id.key$SENT.TUMOR[match(wgs.df$Tumour, id.key$TUMOR.alias)]
wgs.df$Sample <- ifelse(is.na(wgs.df$Sample), wgs.df$Sample2, wgs.df$Sample)
wgs.df$Sample2 <- NULL
#View(wgs.df)

wgs.anno <- merge(wgs.df,ERpHER2n.dat,by="Sample")
wgs.anno.her2e <- wgs.anno %>% filter(NCN.PAM50=="Her2")
View(wgs.anno %>% filter(NCN.PAM50=="Her2") %>% filter(toupper(`Final QC`)=="PASS"))

# cn data
cn.her2e <- names(segment.files)
setdiff(cn.her2e,wgs.anno.her2e$Sample)
setdiff(cn.her2e, her2e.samples)

# scanb anno
d <- loadRData("data/SCANB/1_clinical/processed/SCANB_clinData.RData")
View(d[d$Sample=="S003095",])

# file i sent to johan 16.11.
q <- loadRData("/Users/le7524ho/Desktop/her2samples_wgs_status.RData")
setdiff(q$rel4.Sample,her2e.samples)
setdiff(q$rel4.Sample,her2e.samples.df$Sample)

setdiff(q$rel4.Sample,wgs.anno.her2e$Sample)


sub.drivers <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_coding_and_drivers_3March23.xlsx", sheet = "AllCodingSubs"))
allut <- loadRData("data/SCANB/3_genomic/processed/driver_mutations_all.RData")
View(allut)
setdiff(unique(scanb.dmut$sample),her2e.samples)
setdiff(unique(scanb.dmut$sample),her2e.samples.df)

setdiff(unique(scanb.dmut$sample),cn.her2e)
setdiff(cn.her2e,unique(scanb.dmut$sample))
View()
# included her2e samples
cn.her2e
# file of all ERpos samples, with indication of ones that were included but are not required
all.samples.anno <- wgs.anno %>% 
  mutate(Missing = if_else(NCN.PAM50 =="Her2" & 
                             toupper(`Final QC`) == "PASS" & 
                             !(Sample %in% cn.her2e),1,0)) %>% 
  mutate(Incl_but_not_required = if_else(NCN.PAM50 !="Her2" & 
                                           Sample %in% cn.her2e,1,0))

# file of all samples that should be included, incl. column indicating the ones that are missing
her2e.samples.df <- wgs.anno %>% 
  filter(NCN.PAM50=="Her2") %>% 
  filter(toupper(`Final QC`)=="PASS") %>% 
  mutate(Missing = if_else(!(Sample %in% cn.her2e),1,0))

View(all.samples.anno)
View(her2e.samples.df)
save(all.samples.anno,file="data/SCANB/3_genomic/processed/not_required_samples.RData")
save(her2e.samples.df,file="data/SCANB/3_genomic/processed/her2e_samples.RData")

#######################


r <- loadRData("data/SCANB/2_transcriptomic/processed/mg_anno.RData")
r <- loadRData("data/METABRIC/2_transcriptomic/processed/mg_anno.RData")
head(r)

View(r[[2]])
r <- loadRData("data/SCANB/2_transcriptomic/processed/DE_results.RData")
head(r)
colnames(r)
rownames(r)
nrow(r)

# scnab % of DEGS compared to all
r <- loadRData("data/SCANB/2_transcriptomic/processed/DE_results.RData")
nrow(r[r$Her2.LumA.padj <= 0.05,])
(nrow(r[r$Her2.LumA.padj <= 0.05,]) / nrow(r))*100

nrow(r[r$Her2.LumB.padj <= 0.05,])
(nrow(r[r$Her2.LumB.padj <= 0.05,]) / nrow(r))*100

# Metabric % of DEGS compared to all
r <- loadRData("data/METABRIC/2_transcriptomic/processed/DE_results.RData")
nrow(r)
row.names(r)

nrow(r[r$Her2.LumA.padj <= 0.05,])
(nrow(r[r$Her2.LumA.padj <= 0.05,]) / nrow(r))*100

# double check method in de script - why so few, rerun de script
nrow(r[r$Her2.LumB.padj <= 0.05,])
(nrow(r[r$Her2.LumB.padj <= 0.05,]) / nrow(r))*100



####
r <- loadRData("data/SCANB/4_CN/raw/to_lennart/ascat.S000763_l_d_a_vs_B000930_d_a/S000763_l_d_a_ascat_output_gamma0.9_penalty100.RData")
head(r)
str(r)
r$ploidy
View(r)
View(r$segments)
nrow(r$segments)
View(r$nA)
