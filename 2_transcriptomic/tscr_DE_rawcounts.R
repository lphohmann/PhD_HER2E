# Script: Differential gene expression based on raw counts in SCANB 
# use DESeq.
#BiocManager::install("DESeq2")
library(DESeq2)
#browseVignettes("DESeq2")

#read in a count matrix, which we will name cts, and the sample information table, which we will name coldata
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
#examine the count matrix and column data to see if they are consistent in terms of sample order.
head(cts,2)
coldata
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. also names have to be consistent
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# pre filtering
#minimal pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#Alternatively, a popular filter is to ensure at least X samples with a count of 10 or more, where X can be chosen as the sample size of the smallest group of samples
  
# unevaluated chunk...
keep <- rowSums(counts(dds) >= 10) >= X
dds <- dds[keep,]

# TELL WICH COMPARISONS SHOULD BE MODE USING THE CONTRAST ARGUMENT or order the levels of the $condition


# my part:
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")

#load count data
cts <- as.data.frame(read.csv("data/SCANB/2_transcriptomic/raw/gene_count_matrix-4.csv"))
View(head(cts))
# delete ensg and only keep hugo gene name
cts$gene_id <- sub('^.', '', gsub("^[^|]*|", "", cts$gene_id))
# df with samples as rownames and columns condition, type
coldata <- data.frame(condition=rep(NA,ncol(cts)-1), # this will be the pam50 subtypes
                      type=rep(NA,ncol(cts)-1),
                      row.names = tail(colnames(cts),-1))

# CONTINUE HERE AFTER KNOWING IF PAIRED-END OR SINGLE-READ

#coldata$condition <- factor(coldata$condition)
#coldata$type <- factor(coldata$type)

#examine the count matrix and column data to see if they are consistent in terms of sample order.
head(cts,2)
coldata

#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. also names have to be consistent
cts <- cts[, rownames(coldata)]
