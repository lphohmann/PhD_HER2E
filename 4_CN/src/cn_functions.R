# Function definitions for cn analyses

################################################################################
# function to process data 
################################################################################

# function
processing <- function(cn.data,chr.lengths) {
    
    # rename X chromosome to 23
    if (!is.numeric(cn.data$Chr)) {
        cn.data$Chr[cn.data$Chr == "X"] <- "23" 
        cn.data$Chr <- as.numeric(cn.data$Chr)
    }
    
    # add chr lengths to df
    cn.data <- cn.data %>% group_by(Chr) %>% 
        mutate(chr_genomelengths = 
                   case_when(Chr==1 ~ chr.lengths$genome[1],
                             Chr==2 ~ chr.lengths$genome[2],
                             Chr==3 ~ chr.lengths$genome[3],
                             Chr==4 ~ chr.lengths$genome[4],
                             Chr==5 ~ chr.lengths$genome[5],
                             Chr==6 ~ chr.lengths$genome[6],
                             Chr==7 ~ chr.lengths$genome[7],
                             Chr==8 ~ chr.lengths$genome[8],
                             Chr==9 ~ chr.lengths$genome[9],
                             Chr==10 ~ chr.lengths$genome[10],
                             Chr==11 ~ chr.lengths$genome[11],
                             Chr==12 ~ chr.lengths$genome[12],
                             Chr==13 ~ chr.lengths$genome[13],
                             Chr==14 ~ chr.lengths$genome[14],
                             Chr==15 ~ chr.lengths$genome[15],
                             Chr==16 ~ chr.lengths$genome[16],
                             Chr==17 ~ chr.lengths$genome[17],
                             Chr==18 ~ chr.lengths$genome[18],
                             Chr==19 ~ chr.lengths$genome[19],
                             Chr==20 ~ chr.lengths$genome[20],
                             Chr==21 ~ chr.lengths$genome[21],
                             Chr==22 ~ chr.lengths$genome[22],
                             Chr==23 ~ chr.lengths$genome[23])) %>% 
        relocate(chr_genomelengths, .after=Position) %>% ungroup()
    
    # add up position + chr length to get genome position 
    cn.data$genome_pos <- cn.data$Position + cn.data$chr_genomelengths 
    cn.data <- cn.data %>% relocate(genome_pos, .after=chr_genomelengths)
    cn.data$chr_genomelengths <- NULL
    cn.data$genome <- NULL
    
    return(cn.data)
}

################################################################################
# function to process data post liftover
################################################################################

# function
postliftprocess <- function(hg38.pos.df, fData, CN_Gain, CN_Loss) {
    fData <- as.data.frame(merge(hg38.pos.df, fData, by = "reporterId")) %>% dplyr::rename(position = start) %>% dplyr::select(-c("chromosome","centerPosition"))
    gaindf <- as.data.frame(CN_Gain) %>% rownames_to_column(var="reporterId") %>% dplyr::rename(Gain=2)
    lossdf <- as.data.frame(CN_Loss) %>% rownames_to_column(var="reporterId") %>% dplyr::rename(Loss=2)
    metadf <- as.data.frame(fData)
    combdf <- as.data.frame(merge(merge(gaindf,lossdf,by="reporterId"),metadf, by="reporterId"))
    return(combdf)
}
