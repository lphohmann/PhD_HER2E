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


##

View(head(SCANB.ERpHER2n.data))


## her2p de results
x <- read_table("data/SCANB/2_transcriptomic/processed/DEG_analysis_results_DEGgenes_HER2E_vs_ERpHER2p.txt") %>% as.data.frame()
x$p.val.wilcox.HER2E_vs_ERpHER2pHER2E.fdr <- p.adjust(x$p.val.wilcox.HER2E_vs_ERpHER2pHER2E,
                                                      method = "fdr")
x$p.val.wilcox.HER2E_vs_ERpHER2pNonHER2E.fdr <- p.adjust(x$p.val.wilcox.HER2E_vs_ERpHER2pNonHER2E,
                                                      method = "fdr")

x %>% filter(p.val.wilcox.HER2E_vs_ERpHER2pHER2E.fdr <=0.05)
x %>% filter(p.val.wilcox.HER2E_vs_ERpHER2pNonHER2E.fdr <=0.05)

y <-x %>% filter(p.val.wilcox.HER2E_vs_ERpHER2pNonHER2E.fdr <=0.05) %>% pull(Gene.Name) %>% unique()
length(y)



##########################################################################
# create supplementary files
##########################################################################

# 1. DE results: HER2E SCANB, HER2E METBARIC, CORE SETS 1,2,3, PWE?
DE.res.scanb <- loadRData("./data/SCANB/2_transcriptomic/processed/DE_results.RData") %>% 
  dplyr::select(-c("Her2.LumA.de","Her2.LumB.de"))

DE.res.metabric <- loadRData("./data/METABRIC/2_transcriptomic/processed/DE_results.RData") %>% 
  dplyr::select(-c("Her2.LumA.de","Her2.LumB.de"))

# core sets
# define DEG set
DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene)
DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% pull(Gene) 
# core gex
HER2E.DEGs.core <- intersect(DEGs.scanb,DEGs.metabric)
# define DEG set
LumA.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene)
LumA.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumA.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene) 
# gex
LumA.DEGs.core <- intersect(LumA.DEGs.scanb,LumA.DEGs.metabric)
# define DEG set
LumB.DEGs.scanb <- DE.res.scanb %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene)
LumB.DEGs.metabric <- DE.res.metabric %>% 
  filter(Her2.LumB.padj <= 0.05) %>% 
  rownames_to_column("Gene") %>% 
  pull(Gene) 
# core top gex
LumB.DEGs.core <- intersect(LumB.DEGs.scanb,LumB.DEGs.metabric)


sq <- seq(max(length(HER2E.DEGs.core), length(LumA.DEGs.core), length(LumB.DEGs.core)))
coreset.df <- data.frame(HER2E.core.set = HER2E.DEGs.core[sq],
                         LumA.core.set = LumA.DEGs.core[sq],
                         LumB.core.set = LumB.DEGs.core[sq])

#pwe 
GSEA.HER2E.core.set <- openxlsx::read.xlsx("./output/supplementary_data/HER2n_pwenrichment_WikiPathway_2023_Human.xlsx",
                                           sheet="CoreDEGs")
GSEA.LumA.core.set <- openxlsx::read.xlsx("./output/supplementary_data/HER2n_pwenrichment_WikiPathway_2023_Human.xlsx",
                                           sheet="LumA_coreDEGs")
GSEA.LumB.core.set <- openxlsx::read.xlsx("./output/supplementary_data/HER2n_pwenrichment_WikiPathway_2023_Human.xlsx",
                                           sheet="LumB_coreDEGs")


#View(coreset.df)

# save results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb=wb, sheetName="DE_results_SCANB")
openxlsx::writeDataTable(wb=wb,sheet="DE_results_SCANB",x=DE.res.scanb %>% rownames_to_column(var="Gene"))

openxlsx::addWorksheet(wb=wb, sheetName="DE_results_METABRIC")
openxlsx::writeDataTable(wb=wb,sheet="DE_results_METABRIC",x=DE.res.metabric %>% rownames_to_column(var="Gene"))

openxlsx::addWorksheet(wb=wb, sheetName="Core_gene_sets")
openxlsx::writeDataTable(wb=wb,sheet="Core_gene_sets",keepNA = TRUE,x=coreset.df)

openxlsx::addWorksheet(wb=wb, sheetName="GSEA.HER2E.core.set")
openxlsx::writeDataTable(wb=wb,sheet="GSEA.HER2E.core.set",x=GSEA.HER2E.core.set)

openxlsx::addWorksheet(wb=wb, sheetName="GSEA.LumA.core.set")
openxlsx::writeDataTable(wb=wb,sheet="GSEA.LumA.core.set",x=GSEA.LumA.core.set)

openxlsx::addWorksheet(wb=wb, sheetName="GSEA.LumB.core.set")
openxlsx::writeDataTable(wb=wb,sheet="GSEA.LumB.core.set",x=GSEA.LumB.core.set)

openxlsx::saveWorkbook(wb=wb,file="./output/supplementary_data/Manuscript_files/Table_S2.xlsx",overwrite=TRUE)

# 2. core Set: HER2E core set, cytoscape cluster, PWE, TFBS

CS.cluster1 <- read.table("data/SCANB/2_transcriptomic/processed/DEGs_HER2E_Sub1_network_nCon4_corCut0.5_GENES.txt") %>% 
  pull(V1)
CS.cluster2 <-read.table("data/SCANB/2_transcriptomic/processed/DEGs_HER2E_Sub2_network_nCon4_corCut0.5_GENES.txt") %>% 
  pull(V1)
clusters.df <- as.data.frame(HER2E.DEGs.core) %>% 
  mutate(CS_clust1 = ifelse(HER2E.DEGs.core %in% CS.cluster1,1,0)) %>% 
  mutate(CS_clust2 = ifelse(HER2E.DEGs.core %in% CS.cluster2,1,0))

GSEA.cluster1 <- loadRData("./output/supplementary_data/HER2E_cluster1_PWE_WP2023.RData")
GSEA.cluster2 <- loadRData("./output/supplementary_data/HER2E_cluster2_PWE_WP2023.RData")
TFBS.transfak.cluster1 <- "data/SCANB/2_transcriptomic/processed/TF_Sub1_Transfac_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
TFBS.transfak.cluster2 <- "data/SCANB/2_transcriptomic/processed/TF_Sub2_Transfac_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
TFBS.transfak.all <- "data/SCANB/2_transcriptomic/processed/TF_Lenart_Transfac_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
TFBS.jaspar.all <- "data/SCANB/2_transcriptomic/processed/TF_Lenart_JASPAR_new_tfsearch_TFBS_significant_analysis_file_100000_p(1.0).txt"
TFBS.transfak.cluster1 <- read.table(TFBS.transfak.cluster1,header=TRUE,sep = "\t")
TFBS.transfak.cluster2 <- read.table(TFBS.transfak.cluster2,header=TRUE,sep = "\t")
TFBS.transfak.all <- read.table(TFBS.transfak.all,header=TRUE,sep = "\t")
TFBS.jaspar.all <- read.table(TFBS.jaspar.all,header=TRUE,sep = "\t")

HER2E_vs_ERpHER2p_DE_results <- read.table("./data/SCANB/2_transcriptomic/processed/DEG_analysis_results_DEGgenes_HER2E_vs_ERpHER2p.txt",header = TRUE)

# save results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb=wb, sheetName="HER2E.core.set.Clusters")
openxlsx::writeDataTable(wb=wb, sheet="HER2E.core.set.Clusters", x=clusters.df)

openxlsx::addWorksheet(wb=wb, sheetName="GSEA.cluster1")
openxlsx::writeDataTable(wb=wb,sheet="GSEA.cluster1",x=GSEA.cluster1)

openxlsx::addWorksheet(wb=wb, sheetName="GSEA.cluster2")
openxlsx::writeDataTable(wb=wb,sheet="GSEA.cluster2",x=GSEA.cluster2)

openxlsx::addWorksheet(wb=wb, sheetName="HER2E_vs_ERpHER2p_DE_results")
openxlsx::writeDataTable(wb=wb,sheet="HER2E_vs_ERpHER2p_DE_results",x=HER2E_vs_ERpHER2p_DE_results)

openxlsx::addWorksheet(wb=wb, sheetName="TFBS.jaspar.all")
openxlsx::writeDataTable(wb=wb,sheet="TFBS.jaspar.all",x=TFBS.jaspar.all)

openxlsx::addWorksheet(wb=wb, sheetName="TFBS.transfak.all")
openxlsx::writeDataTable(wb=wb,sheet="TFBS.transfak.all",x=TFBS.transfak.all)

openxlsx::addWorksheet(wb=wb, sheetName="TFBS.transfak.cluster1")
openxlsx::writeDataTable(wb=wb,sheet="TFBS.transfak.cluster1",x=TFBS.transfak.cluster1)

openxlsx::addWorksheet(wb=wb, sheetName="TFBS.transfak.cluster2")
openxlsx::writeDataTable(wb=wb,sheet="TFBS.transfak.cluster2",x=TFBS.transfak.cluster2)

openxlsx::saveWorkbook(wb=wb,file="./output/supplementary_data/Manuscript_files/Table_S3.xlsx",overwrite=TRUE)

# 3. CNA: signif. genes  One supp table dumping the significant CNA genes per comparison.
HER2E_vs_LumA_LumB_CNA_results <- loadRData("data/COMBINED/4_CN/processed/CNA_genelevel.RData")

#colnames(genes)
luma.genes.loss <- HER2E_vs_LumA_LumB_CNA_results %>% 
  filter(LumA.Loss.padj <= 0.05) %>% 
  pull(gene)
luma.genes.gain <- HER2E_vs_LumA_LumB_CNA_results %>% 
  filter(LumA.Gain.padj <= 0.05) %>% 
  pull(gene)
lumb.genes.loss <- HER2E_vs_LumA_LumB_CNA_results %>% 
  filter(LumB.Loss.padj <= 0.05) %>% 
  pull(gene)
lumb.genes.gain <- HER2E_vs_LumA_LumB_CNA_results %>% 
  filter(LumB.Gain.padj <= 0.05) %>% 
  pull(gene)

sq <- seq(max(length(luma.genes.loss),length(lumb.genes.loss)))

Genes_signif_LossFreqDiff <- data.frame(HER2E_vs_LumA=luma.genes.loss[sq],
                                       HER2E_vs_LumB=lumb.genes.loss[sq])
sq <- seq(max(length(luma.genes.gain),length(lumb.genes.gain)))
Genes_signif_GainFreqDiff <- data.frame(HER2E_vs_LumA=luma.genes.gain[sq],
                                       HER2E_vs_LumB=lumb.genes.gain[sq])

# save results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb=wb, sheetName="HER2E_vs_LumA_LumB_CNA_results")
openxlsx::writeDataTable(wb=wb, sheet="HER2E_vs_LumA_LumB_CNA_results", x=HER2E_vs_LumA_LumB_CNA_results)

openxlsx::addWorksheet(wb=wb, sheetName="Genes_signif_LossFreqDiff")
openxlsx::writeDataTable(wb=wb,sheet="Genes_signif_LossFreqDiff",x=Genes_signif_LossFreqDiff)

openxlsx::addWorksheet(wb=wb, sheetName="Genes_signif_GainFreqDiff")
openxlsx::writeDataTable(wb=wb,sheet="Genes_signif_GainFreqDiff",x=Genes_signif_GainFreqDiff)

openxlsx::saveWorkbook(wb=wb,file="./output/supplementary_data/Manuscript_files/Table_S4.xlsx",overwrite=TRUE)
#



# repository file prep
# prep all sheets with correct sampleIDs
Samples <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                               sheet="Samples") %>% mutate(Final.QC="PASS")
idkey <- Samples[1:2]
Samples <- Samples %>% dplyr::select(-c("Tumour","Missing","X9","NCN.PAM50")) %>% mutate(PAM50="HER2E")

HRDetect <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                               sheet="HRDetect") %>% 
  mutate(Call = ifelse(Probability >= 0.7,"HRD-high","HRD-low"))
HRDetect$Sample <- idkey$Sample[match(HRDetect$Sample,idkey$Tumour)]
HRDetect <- HRDetect %>% dplyr::select(-c("Sample_pair","Normal"))

Rearrangement_Sig <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                               sheet="RearrangementSig")
Rearrangement_Sig$Sample <- idkey$Sample[match(Rearrangement_Sig$Sample,idkey$Tumour)]
Rearrangement_Sig <- Rearrangement_Sig %>% dplyr::select(-c("Sample_pair","Normal"))

SBS_Sig <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                               sheet="SBSsig")
SBS_Sig$Sample <- idkey$Sample[match(SBS_Sig$Sample,idkey$Tumour)]
SBS_Sig <- SBS_Sig %>% dplyr::select(-c("Sample_pair","X3"))

Caveman_Drivers <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                               sheet="CavemanDrivers")
Caveman_Drivers$Sample <- idkey$Sample[match(Caveman_Drivers$Sample,idkey$Tumour)]
Caveman_Drivers <- Caveman_Drivers %>% dplyr::select(-c("Normal"))

Pindel_Drivers <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                               sheet="PindelDrivers")
Pindel_Drivers$Sample <- idkey$Sample[match(Pindel_Drivers$Sample,idkey$Tumour)]
Pindel_Drivers <- Pindel_Drivers %>% dplyr::select(-c("Normal"))

BRASS_Drivers <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                                         sheet="BRASS_dirvers")
BRASS_Drivers$sample <- idkey$Sample[match(BRASS_Drivers$Tumour,idkey$Tumour)]
BRASS_Drivers <- BRASS_Drivers %>% dplyr::select(-c("Tumour")) %>% dplyr::rename(Sample=sample) %>% mutate(QC="PASS")

AllCodingASMD_CLPM <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                               sheet="AllCodingASMD_CLPM")
AllCodingASMD_CLPM$Sample <- idkey$Sample[match(AllCodingASMD_CLPM$Sample,idkey$Tumour)]
AllCodingASMD_CLPM <- AllCodingASMD_CLPM %>% dplyr::select(-c("Normal"))

AllCodingPindel <- openxlsx::read.xlsx("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx",
                                      sheet="AllCodingPindel")
AllCodingPindel$Sample <- idkey$Sample[match(AllCodingPindel$Sample, idkey$Tumour)]
AllCodingPindel <- AllCodingPindel %>% dplyr::select(-c("Normal"))



# kataegis and TMB
MutationCounts <- loadRData("./data/SCANB/3_genomic/processed/mutation_counts.RData") %>% relocate(PAM50,.after=N_mut) %>% 
  relocate(sampleID,.before=caveman_count) %>% dplyr::rename(Sample=sampleID)

Kataegis <- loadRData("./data/SCANB/4_CN/processed/CN_kataegis.RData") %>% 
  filter(PAM50=="HER2E") %>% 
  relocate(PAM50,.after=Kat_binary) %>% 
  dplyr::rename(Sample=SampleID)

# get segment data with CN states
x <- loadRData("./data/SCANB/4_CN/processed/Segment_CN_states.RData")
# add sample name as column for each
x <- lapply(seq_along(x), y=x, n=names(x), function(i, y, n) {
  df <- y[[i]]
  sample <- n[[i]]
  df <- df %>% mutate(Sample = sample) %>% relocate(Sample,.before=chr)
  return(df)
})
# rbind together with sampleID column
CNA_segments <- do.call("rbind", x)

wb <- openxlsx::createWorkbook()
for(i in c("Samples","HRDetect","Rearrangement_Sig","SBS_Sig",
           "Caveman_Drivers","Pindel_Drivers","BRASS_Drivers",
           "AllCodingASMD_CLPM","AllCodingPindel","MutationCounts","Kataegis","CNA_segments")) {
  
  openxlsx::addWorksheet(wb=wb, sheetName=i)
  openxlsx::writeDataTable(wb=wb, sheet=i, x=eval(as.name(paste(i))))
  
}
openxlsx::saveWorkbook(wb=wb,file="./output/supplementary_data/Manuscript_files/Table_S1.xlsx",overwrite=TRUE)
