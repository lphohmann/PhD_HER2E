# Script: prep. source file
# Author: Lennart Hohmann
# Date: 03.10.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_HER2E/")
#-------------------
# packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl,openxlsx)
#-------------------
# set/create output directories
output.path <- "./output/source_data/Excel_files/"
dir.create(output.path)
#-------------------
# input paths
infile.1 <- "./output/source_data/R_objects/Table_1_Metabric.RData"
infile.2 <- "./output/source_data/R_objects/Table_1_SCANB.RData"
infile.3 <- "./output/source_data/R_objects/Figure_2_SCANB.RData"
infile.4 <- "./output/source_data/R_objects/Figure_2_Metabric.RData"
infile.5 <- "./data/SCANB/2_transcriptomic/processed/mg_anno.RData"
infile.6 <- "./data/SCANB/2_transcriptomic/processed/mg_anno_HER2p.RData"
infile.7 <- "./output/source_data/R_objects/Figure_3_singlegex.RData"
infile.8 <- "./output/source_data/R_objects/Figure_3_umap.RData" #umap
infile.9 <- "./output/source_data/R_objects/Figure_3_PAM50corr.RData" # centroid correlations
# output paths
infile.10 <- "./output/source_data/R_objects/Figure_4_sigs.RData"
source.file <- "./output/source_data/source_file.xlsx" # excel file extension

################################################################################

loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
}

add_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)   # Add a new worksheet with the specified name
  writeData(wb, sheet_name, data)  # Write the provided data frame to the worksheet
}

################################################################################
# Create a new workbook
wb <- createWorkbook()

# table 1
t1.SCANB <- loadRData(infile.1) #mb
add_sheet(wb, "Table 1 -SCAN-B", t1.SCANB)
t1.METABRIC <- loadRData(infile.2) #scanb
add_sheet(wb, "Table 1 - METABRIC", t1.METABRIC)

# figure 2
f2.ACEG <- loadRData(infile.3)
add_sheet(wb, "Figure 2 - Panels A,C,E,G", f2.ACEG)
f2.BDFH <- loadRData(infile.4)
add_sheet(wb, "Figure 2 - Panels B,D,F,H", f2.BDFH)

# figure 3
f3.ABC <- loadRData(infile.5)[[1]]
f3.ABC <- f3.ABC[c("sampleID","PAM50","SR","IR","Mitotic_progression")]
add_sheet(wb, "Figure 3 - Panels A-C", f3.ABC)
f3.D <-  data.frame(Cell.type = c("B.cells","CD8.T.cells","CD4.T.cells","NK.cells","Monocytes","Neutrophils"),ir.ciber.cor=c(0.32475404,0.51323221,0.02367644,0.16116747,0.52607297,-0.02584477))
add_sheet(wb, "Figure 3 - Panel D", f3.D)

f3.E <- read_excel("./output/source_data/Excel_files/Fig3_sourceData.xlsx")
add_sheet(wb, "Figure 3 - Panel E", f3.E)
f3.FGHI <- loadRData(infile.7)
add_sheet(wb, "Figure 3 - Panels F-I", f3.FGHI)
f3.J <- loadRData(infile.8)
add_sheet(wb, "Figure 3 - Panel J", f3.J)
f3.KL <- loadRData(infile.9)
add_sheet(wb, "Figure 3 - Panels K,L", f3.KL)

# figure 4
f4.ABCDE <- loadRData(infile.10) # maybe get basis sample ids as well
add_sheet(wb, "Figure 4 - Panels A-E", f4.ABCDE)
f4.F <- loadRData("./output/source_data/R_objects/Figure_4_TMB.RData")
add_sheet(wb, "Figure 4 - Panel F", f4.F)
f4.G <- loadRData("./data/SCANB/3_genomic/processed/driver_mutations_all.RData")
add_sheet(wb, "Figure 4 - Panel G", f4.G)
f4.HI <- loadRData("./output/source_data/R_objects/Figure_4_mutfreqs.RData")
add_sheet(wb, "Figure 4 - Panels H,I", f4.HI)
f4.J <- loadRData("./output/source_data/R_objects/Figure_4_genomealtered.RData")
add_sheet(wb, "Figure 4 - Panel J", f4.J)
f4.K <- loadRData("./output/source_data/R_objects/Figure_4_genediff.RData")
add_sheet(wb, "Figure 4 - Panel K", f4.K)
f4.L <- read_excel("./output/supplementary_data/Manuscript_files/Table_S7.xlsx", sheet = "HER2E_vs_LumA_LumB_CNA_results")
add_sheet(wb, "Figure 4 - Panel L", f4.L)

# figure 5
f5.ABC <- loadRData(infile.6)[[1]]
f5.ABC <- f5.ABC[c("sampleID","Group","SR","IR","Mitotic_progression")]
add_sheet(wb, "Figure 5 - Panels A-C", f5.ABC)
f5.DEF <- loadRData("./output/source_data/R_objects/Figure_5_singlegexHER2p.RData")
add_sheet(wb, "Figure 5 - Panels D-F", f5.DEF)
f5.GH <- loadRData("./output/source_data/R_objects/Figure_5_mutfreqsHER2p.RData")
add_sheet(wb, "Figure 5 - Panels G,H", f5.GH)
f5.I <- loadRData("./output/source_data/R_objects/Figure_5_umapHER2p.RData")
add_sheet(wb, "Figure 5 - Panel I", f5.I)
f5.J <- loadRData("./output/source_data/R_objects/Figure_5_wfHER2p.RData")
add_sheet(wb, "Figure 5 - Panel J", f5.J)

# figure 6
existing_wb <- loadWorkbook("./output/source_data/Excel_files/Fig6_sourceData.xlsx")

# Copy sheets from the existing workbook to the combined workbook
for (sheet_name in getSheetNames("./output/source_data/Excel_files/Fig6_sourceData.xlsx")) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, read.xlsx("./output/source_data/Excel_files/Fig6_sourceData.xlsx", sheet = sheet_name))
}

# Save the workbook to an Excel file
saveWorkbook(wb, source.file, overwrite = TRUE)
