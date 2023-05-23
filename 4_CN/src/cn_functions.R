# Function definitions for cn analyses

################################################################################
# function to process data 
################################################################################

# loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
}

################################################################################
# function to add genome position of probes
################################################################################

add_genomepos <- function(cn.data, chr.lengths) {
  res <- cn.data %>% 
    # add new chr genome position column
    mutate(genome = 0) %>% 
    # update the genome col to fill in the actual chr positions
    rows_update(chr.lengths[c("Chr","genome")]) %>% 
    # add a column with the genome position of each probe
    mutate(Genome_pos = Position + genome) %>% 
    relocate(c(genome,Genome_pos), .after=Position) %>% 
    dplyr::select(-c(genome))
  
  return(res)
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
