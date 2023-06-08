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