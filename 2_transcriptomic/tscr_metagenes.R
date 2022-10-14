# Script: Metagene analysis in SCAN-B and METABRIC

# TODO:

# empty environment
rm(list=ls())

# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

# indicate for which cohort the analysis is run 
cohort <- "SCANB" # SCANB or METABRIC

# set/create output directory for plots
output.path <- "output/plots/2_transcriptomic/"
dir.create(output.path)

# set/create output directory for processed data
data.path <- paste("data/",cohort,"/2_transcriptomic/processed/",sep="")
dir.create(data.path)

#packages
source("scripts/2_transcriptomic/src/tscr_functions.R")
library(ggplot2)
library(tidyverse)
#library(matrixStats)
library(readxl)
library(biomaRt)
library(gridExtra)
library(ggsignif)

#if (cohort=="SCANB") {
    
#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# load annotation data
clin.rel4 <- as.data.frame(
    read_excel("data/SCANB/1_clinical/raw/NPJ_release.xlsx"))
# load gex data
load("data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata")

anno <- clin.rel4 %>% 
    filter(Follow.up.cohort==TRUE) %>% 
    filter(NCN.PAM50 %in% c("LumA", "LumB", "Her2")) %>% 
    filter(ER=="Positive" & HER2=="Negative") %>% 
    dplyr::rename(sampleID = GEX.assay, PAM50 = NCN.PAM50)

gex.data <- as.data.frame(genematrix_noNeg[,colnames(genematrix_noNeg) %in% anno$sampleID])

#######################################################################
# 3. data processing # IDEA SCALE OVER WHOLE COHORT IN THE BEGINNING
#######################################################################

# log transformed FPKM data
gex.data <- as.data.frame(log2(gex.data + 1))

# sample annotation
sample.info <- anno[c("sampleID","PAM50")] %>% remove_rownames %>% column_to_rownames(var="sampleID") #%>% dplyr::rename(group=PAM50)
    
# metagene definitions
metagene.def <- as.data.frame(read_excel(
    "data/SCANB/2_transcriptomic/raw/metagene_definitions.XLSX")) %>% 
    dplyr::rename(entrezgene_id = `Entrez Gene ID`,
                  module = `Module Name`,
                  gene_symbol = `Gene symbol`)

# Convert ensembl to entrez IDs
# select ensembl ids
ensembl.ids <- as.data.frame(gex.data) %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    pull(ensembl_gene_id) 
ensembl.ids <- gsub("\\..*","",ensembl.ids) # remove characters after dot

# convert 
#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                         dataset = "hsapiens_gene_ensembl",
#                         host = "http://www.ensembl.org")
#res <- getBM(filters = "ensembl_gene_id",
#                 attributes = c("ensembl_gene_id","entrezgene_id"),
#                 values = ensembl.ids, 
#                 mart = mart)
#save(res,file = paste(data.path,"mart_res.RData",sep="") )
load(paste(data.path,"mart_res.RData",sep=""))
# select only the ids that are relevant for the metagenes
relevant.res <- res %>% 
    filter(entrezgene_id %in% metagene.def$entrezgene_id) #pull(ensembl_gene_id)

# no duplicates
#length(relevant.ensembl.ids)
#length(unique(relevant.ensembl.ids))

# genes in the metag def
m.genes <- metagene.def %>% pull(gene_symbol)

# merge
metagene.def <- merge(metagene.def, relevant.res, by = "entrezgene_id")

# check if genes were not in the scanb data
setdiff(m.genes,metagene.def$gene_symbol) #"IGHM" "TRAC"

# get the gex data in the correct format
gex.data <- gex.data %>% 
    rownames_to_column(var="ensembl_gene_id") %>% 
    mutate(ensembl_gene_id = gsub("\\..*","",ensembl_gene_id)) %>% 
    filter(ensembl_gene_id %in% metagene.def$ensembl_gene_id) %>% 
    column_to_rownames(var="ensembl_gene_id")

# scale the data
gex.data <- gex.data %>% select_if(~ !any(is.na(.))) # exclude column iwth NA
scaled.gex.data <- as.data.frame(t(apply(gex.data, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) %>% # for some rows there may be 0 variance so i have to handle these cases
    rownames_to_column(var="ensembl_gene_id") 

#######################################################################
# 4. calc. score for each metagene in each sample
#######################################################################

## remove later
source("scripts/2_transcriptomic/src/tscr_functions.R")

#test
basal.scores <- mgscore(metagene = "Basal",
                metagene.def = metagene.def,
                gex.data = scaled.gex.data) 

earlyresponse.scores <- mgscore(metagene = "Early_response",
                        metagene.def = metagene.def,
                        gex.data = scaled.gex.data)

m <- merge(basal.scores,earlyresponse.scores,by=0) %>% column_to_rownames(var = "Row.names")

head(m)
# basal
scaled_mg_scores <- mg_score("Basal",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Basal = mg)

er <- mg_score("Early_response",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Early_response = mg)
scaled_mg_scores <- merge(scaled_mg_scores,er,by=0) %>% column_to_rownames(var = "Row.names")

ir <- mg_score("IR",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(IR = mg)
scaled_mg_scores <- merge(scaled_mg_scores,ir,by=0) %>% column_to_rownames(var = "Row.names")

lp <- mg_score("Lipid",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Lipid = mg)
scaled_mg_scores <- merge(scaled_mg_scores,lp,by=0) %>% column_to_rownames(var = "Row.names")

mitc <- mg_score("Mitotic_checkpoint",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Mitotic_checkpoint = mg)
scaled_mg_scores <- merge(scaled_mg_scores,mitc,by=0) %>% column_to_rownames(var = "Row.names")

mitp <- mg_score("Mitotic_progression",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Mitotic_progression = mg)
scaled_mg_scores <- merge(scaled_mg_scores,mitp,by=0) %>% column_to_rownames(var = "Row.names")

sr <- mg_score("SR",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(SR = mg)
scaled_mg_scores <- merge(scaled_mg_scores,sr,by=0) %>% column_to_rownames(var = "Row.names")

str <- mg_score("Stroma",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Stroma = mg)
scaled_mg_scores <- merge(scaled_mg_scores,str,by=0) %>% column_to_rownames(var = "Row.names")

#View(scaled_mg_scores)

pvalues <- mg_ttest(scaled_mg_scores,mg_anno)

#View(pvalues)
pvalues <- pvalues %>% mutate_if(is.numeric, round, digits=5)

#######################################################################
# 5. compare scaled metagene scores between groups
#######################################################################

# add the pam50 annotations
mg_anno <- mg_anno %>% column_to_rownames(var="sampleID")
scaled_mg_scores <- merge(scaled_mg_scores,mg_anno,by=0) %>% column_to_rownames(var = "Row.names")

# save the file
save(scaled_mg_scores, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))
#load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))

# create pdf with all the boxplots
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/metagene_plots.pdf", sep =""), onefile= TRUE)
for(i in c(1:8)) {
    plot <- ggplot(scaled_mg_scores, aes(x=as.factor(PAM50),y=scaled_mg_scores[,i])) +
        geom_boxplot(fill="slateblue",alpha=0.2) +
        xlab("PAM50 subtype") +
        ylab("scaled metagene score") +
        ggtitle(colnames(scaled_mg_scores)[i]) +
        geom_text(data=as.data.frame(dplyr::count(x=mg_anno, PAM50)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -2,nudge_x = 0.3,size=5) +
        theme(axis.text.x = element_text(size = 15),
              axis.title.x = element_text(size = 20),
              axis.text.y = element_text(size = 15),
              axis.title.y = element_text(size = 20))
    print(plot)
}
dev.off()

#geom_text(data=as.data.frame(dplyr::count(sample_info, group)), aes(y = 0, label = paste("n=",n,sep = ""))))

#######################################


# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant
pvalues[1,1:2]
quickplot(scaled_mg_scores, 1, "****", "****", 3.5, c(-1.5,4))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[1],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[2,1:2]
quickplot(scaled_mg_scores, 2, "****", "*", 4, c(-2,4.6))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[2],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")
pvalues[3,1:2]
quickplot(scaled_mg_scores, 3, "****", "****", 5, c(-2,5.6))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[3],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")
pvalues[4,1:2]
quickplot(scaled_mg_scores, 4, "****", "ns", 3.5, c(-1.5,4))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[4],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")
pvalues[5,1:2]
quickplot(scaled_mg_scores, 5, "****", "ns", 4.5, c(-2,5))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[5],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")
pvalues[6,1:2]
quickplot(scaled_mg_scores, 6, "****", "ns", 4.5, c(-2,5))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[6],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")
pvalues[7,1:2]
quickplot(scaled_mg_scores, 7, "****", "****", 2.7, c(-2,3.2))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[7],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")
pvalues[8,1:2]
quickplot(scaled_mg_scores, 8, "ns", "****", 2.9, c(-3.5,3.5))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[8],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

#dev.off()



###############
###############
###############
###############
###############
###############
###############
###############



} else if (cohort=="Metabric") {
# load annotation data
load("Data/Metabric/Annotations/Merged_annotations.RData")
anno <- as.data.frame(anno) %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% filter(grepl('ERpHER2n', ClinGroup)) %>% dplyr::rename(sampleID=METABRIC_ID) %>% dplyr::select(sampleID,PAM50)
# load gex data
# load gex data
gex_data <- as.data.frame(read.table("Data/Metabric/data_mRNA_median_all_sample_Zscores.txt", sep="\t")) %>% row_to_names(row_number = 1) # samples = 1906 - genes = 24368

# remove rowns with na and duplicates
gex_data <- na.omit(gex_data)
gex_data <- gex_data[!duplicated(gex_data$Hugo_Symbol),]

# also save the gene annotation data
gene_anno <- gex_data[,1:2]
gex_data <- gex_data %>% remove_rownames() %>% column_to_rownames(var="Hugo_Symbol") %>% dplyr::select(-Entrez_Gene_Id)
gex_data <- gex_data[,colnames(gex_data) %in% anno$sampleID]

# remove samples from anno that are not in the gex data
anno <- anno %>% filter(sampleID %in% colnames(gex_data))
anno <- anno %>% remove_rownames %>% column_to_rownames(var="sampleID")

gex_data <- rownames_to_column(gex_data, "ensembl_gene_id") # name the column ensembl even though they arent so i can reuse the same functions as for the other cohorts


#######################################################################
# 3. Required data
#######################################################################

metag_def <- as.data.frame(read_excel("Data/metagene_definitions.XLSX")) %>% 
    dplyr::rename(entrezgene_id = `Entrez Gene ID`,
                  module = `Module Name`,
                  gene_symbol = `Gene symbol`)

relevant_ensembl_ids <- gex_data %>% filter(ensembl_gene_id %in% metag_def$gene_symbol) %>% pull(ensembl_gene_id)

#cheeky lil operator
'%!in%' <- function(x,y)!('%in%'(x,y))

#######################################################################
# 4. calc. score for each metagene in each sample
#######################################################################

# function to calc. the p value for each metagene
mg_ttest <- function(scaled_mg_scores,sample_info) {
    # sample ids 
    Her2_cols <- sample_info %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="Her2") %>% pull(sampleID)
    LumA_cols <- sample_info %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumA") %>% pull(sampleID)
    LumB_cols <- sample_info %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumB") %>% pull(sampleID) 
    
    #transpose
    tmg <- t(scaled_mg_scores)
    #initialize storing vector
    H2vsLA_pvalue <- rep(0,nrow(tmg))
    H2vsLB_pvalue <- rep(0,nrow(tmg))
    #loop
    for (i in 1:nrow(tmg)) {
        # vars
        hdata <- as.numeric(tmg[i,Her2_cols])
        adata <- as.numeric(tmg[i,LumA_cols])
        bdata <- as.numeric(tmg[i,LumB_cols])
        
        # for Her2 vs LumA
        # equal variance check
        if (var.test(unlist(hdata),unlist(adata), alternative = "two.sided")$p.value <= 0.05) {
            H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = FALSE)
        } else {
            H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = TRUE)
        }
        # save results
        H2vsLA_pvalue[i] <- H2vsLA_ttest_result$p.value
        
        # for Her2 vs LumB
        # equal variance check
        if (var.test(unlist(hdata),unlist(bdata), alternative = "two.sided")$p.value <= 0.05) {
            H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = FALSE)
        } else {
            H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = TRUE)
        }
        # save results
        H2vsLB_pvalue[i] <- H2vsLB_ttest_result$p.value }
    
    # add the results to the dataframe
    results <- as.data.frame(tmg) %>% add_column(H2vsLA_pvalue = H2vsLA_pvalue,H2vsLB_pvalue = H2vsLB_pvalue) %>% dplyr::select(H2vsLA_pvalue,H2vsLB_pvalue)
    return(results)
}

# function to calc. the score for each metagene for each sample
mg_score <- function(mg,method="scaled",metag_def, gex,cohort) {
    # extract the gex for each metagene
    #method <- "simple"
    mg_ids <- filter(metag_def,module == mg) %>% pull(gene_symbol)
    print(mg_ids)
    mg_gex <- as.data.frame(gex) %>% filter(ensembl_gene_id %in% mg_ids) %>% column_to_rownames(var = "ensembl_gene_id")
    mg_gex <- mg_gex %>% mutate(across(where(is.character), as.numeric))
    if (method == "simple") {
        result <- as.data.frame(apply(mg_gex, 2, median)) %>% dplyr::rename(mg = "apply(mg_gex, 2, median)")
    } else if (method == "scaled") {
        # first scale the rows
        mg_gex <- mg_gex %>% select_if(~ !any(is.na(.))) # exclude column iwth NA
        scaled_mg_gex <- as.data.frame(t(apply(mg_gex, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases
        result <- as.data.frame(apply(scaled_mg_gex, 2, median)) %>% dplyr::rename(mg = "apply(scaled_mg_gex, 2, median)") 
    }
    return(result)
}

#######################################################

# exclude samples that had NA for any of the genes in the metagenes
mg_gex_data <- as.data.frame(t(gex_data)) %>% row_to_names(1)
mg_gex_data <- as.data.frame(t(mg_gex_data[,colnames(mg_gex_data) %in% metag_def$gene_symbol]))
nasamples <- mg_gex_data %>% select_if(~any(is.na(.))) %>% colnames()
mg_gex_data <- mg_gex_data[,colnames(mg_gex_data) %!in% nasamples] %>% rownames_to_column(var="ensembl_gene_id")

mg_anno <- as.data.frame(t(anno)) 
mg_anno <- as.data.frame(t(mg_anno[,colnames(mg_anno) %!in% nasamples,]))
##

# basal
scaled_mg_scores <- mg_score("Basal",metag_def=metag_def, gex=mg_gex_data,cohort=cohort,method = "simple") %>% dplyr::rename(Basal = mg)

er <- mg_score("Early_response",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Early_response = mg)
scaled_mg_scores <- merge(scaled_mg_scores,er,by=0) %>% column_to_rownames(var = "Row.names")

ir <- mg_score("IR",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(IR = mg)
scaled_mg_scores <- merge(scaled_mg_scores,ir,by=0) %>% column_to_rownames(var = "Row.names")

lp <- mg_score("Lipid",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Lipid = mg)
scaled_mg_scores <- merge(scaled_mg_scores,lp,by=0) %>% column_to_rownames(var = "Row.names")

mitc <- mg_score("Mitotic_checkpoint",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Mitotic_checkpoint = mg)
scaled_mg_scores <- merge(scaled_mg_scores,mitc,by=0) %>% column_to_rownames(var = "Row.names")

mitp <- mg_score("Mitotic_progression",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Mitotic_progression = mg)
scaled_mg_scores <- merge(scaled_mg_scores,mitp,by=0) %>% column_to_rownames(var = "Row.names")

sr <- mg_score("SR",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(SR = mg)
scaled_mg_scores <- merge(scaled_mg_scores,sr,by=0) %>% column_to_rownames(var = "Row.names")

str <- mg_score("Stroma",metag_def=metag_def, gex=mg_gex_data,cohort=cohort) %>% dplyr::rename(Stroma = mg)
scaled_mg_scores <- merge(scaled_mg_scores,str,by=0) %>% column_to_rownames(var = "Row.names")

#View(scaled_mg_scores)
pvalues <- mg_ttest(scaled_mg_scores,mg_anno)
#View(head(scaled_mg_scores))
pvalues <- pvalues %>% mutate_if(is.numeric, round, digits=5)

#######################################################################
# 5. compare scaled metagene scores between groups
#######################################################################

# add the pam50 annotations
scaled_mg_scores <- merge(scaled_mg_scores,mg_anno,by=0) %>% column_to_rownames(var = "Row.names")

# save the file
save(scaled_mg_scores, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))
#load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))

# create pdf with all the boxplots
pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/metagene_plots.pdf", sep =""), onefile= TRUE)
for(i in c(1:8)) {
    plot <- ggplot(scaled_mg_scores, aes(x=as.factor(PAM50),y=scaled_mg_scores[,i])) +
        geom_boxplot(fill="slateblue",alpha=0.2) +
        xlab("PAM50 subtype") +
        ylab("scaled metagene score") +
        ggtitle(colnames(scaled_mg_scores)[i]) +
        geom_text(data=as.data.frame(dplyr::count(x=mg_anno, PAM50)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -2,nudge_x = 0.3,size=5) +
        theme(axis.text.x = element_text(size = 15),
              axis.title.x = element_text(size = 20),
              axis.text.y = element_text(size = 15),
              axis.title.y = element_text(size = 20))
    print(plot)
}
dev.off()


#geom_text(data=as.data.frame(dplyr::count(sample_info, group)), aes(y = 0, label = paste("n=",n,sep = ""))))

#######################################
# CONT HERE

# ill do it here without a loop to make the figures nice for publication
#pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/metagene_plots.pdf", sep =""), onefile= TRUE)

# 
quickplot <- function(scaled_mg_scores, plotnum, A_signs, B_signs, B_pos, ylim) {
    plot <- ggplot(scaled_mg_scores, aes(x=as.factor(mg_anno$PAM50),y=scaled_mg_scores[,plotnum],fill=as.factor(mg_anno$PAM50))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("PAM50 subtype") +
        ylab("scaled metagene score") +
        ylim(ylim) +
        ggtitle(colnames(scaled_mg_scores)[plotnum]) +
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations=B_signs, tip_length = 0.02, vjust=0.01, y_position = B_pos, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations=A_signs, tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) + 
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              legend.position = "none")
    print(plot)
}

# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant
pvalues[1,1:2]
quickplot(scaled_mg_scores, 1, "****", "ns", 4, c(-1.5,4.5))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[1],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[2,1:2]
quickplot(scaled_mg_scores, 2, "****", "*", 2.7, c(-2,3.2))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[2],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[3,1:2]
quickplot(scaled_mg_scores, 3, "**", "*", 3.5, c(-2,4))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[3],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[4,1:2]
quickplot(scaled_mg_scores, 4, "****", "*", 3.4, c(-1.5,3.9))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[4],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[5,1:2]
quickplot(scaled_mg_scores, 5, "****", "ns", 3.1, c(-2,3.6))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[5],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[6,1:2]
quickplot(scaled_mg_scores, 6, "****", "ns", 3.4, c(-2,3.9))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[6],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[7,1:2]
quickplot(scaled_mg_scores, 7, "****", "****", 2.5, c(-2,3))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[7],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")

pvalues[8,1:2]
quickplot(scaled_mg_scores, 8, "***", "ns", 2.8, c(-3.5,3.3))
ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/",colnames(scaled_mg_scores)[8],"_metagene_expression.pdf", sep =""),width = 300,height = 300,units = "mm")#dev.off()

    
    
    
}