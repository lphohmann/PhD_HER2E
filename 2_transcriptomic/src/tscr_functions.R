# Function definitions for transcriptomic data analyses

################################################################################
# functions
################################################################################

# define "not in" operator 
'%!in%' <- function(x,y)!('%in%'(x,y))

# function to calc. the score for each metagene for each sample
mg_score <- function(mg,method="scaled",metag_def,gex,cohort) {
    
    # extract the gex for each metagene
    mg_ids <- filter(metag_def, module == mg) %>% pull(ensembl_gene_id)
    
    mg_gex <- gex %>% 
        filter(ensembl_gene_id %in% mg_ids) %>% 
        column_to_rownames(var = "ensembl_gene_id") 
    
    # calc. the score for each sample
    result <- as.data.frame(apply(scaled_mg_gex, 2, median)) %>% dplyr::rename(mg = "apply(scaled_mg_gex, 2, median)") 
    
    return(result)
}

# function to calc. the p value for each metagene
mg_ttest <- function(scaled_mg_scores,sample_info) {
    # sample ids 
    Her2_cols <- sample_info %>% filter(PAM50=="Her2") %>% pull(sampleID)
    LumA_cols <- sample_info %>% filter(PAM50=="LumA") %>% pull(sampleID)
    LumB_cols <- sample_info %>% filter(PAM50=="LumB") %>% pull(sampleID) 
    #transpose
    tmg <- t(scaled_mg_scores)
    #initialize storing vector
    H2vsLA_pvalue <- rep(0,nrow(tmg))
    H2vsLB_pvalue <- rep(0,nrow(tmg))
    #loop
    for (i in 1:nrow(tmg)) {
        # vars
        print(colnames(tmg))
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
        #print(rownames(tmg)[i])
        #print(H2vsLA_ttest_result)
        H2vsLA_pvalue[i] <- H2vsLA_ttest_result$p.value
        
        # for Her2 vs LumB
        # equal variance check
        if (var.test(unlist(hdata),unlist(bdata), alternative = "two.sided")$p.value <= 0.05) {
            H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = FALSE)
        } else {
            H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = TRUE)
        }
        # save results
        #print(rownames(tmg)[i])
        #print(H2vsLB_ttest_result)
        H2vsLB_pvalue[i] <- H2vsLB_ttest_result$p.value }
    
    # add the results to the dataframe
    results <- as.data.frame(tmg) %>% add_column(H2vsLA_pvalue = H2vsLA_pvalue,H2vsLB_pvalue = H2vsLB_pvalue) %>% dplyr::select(H2vsLA_pvalue,H2vsLB_pvalue)
    return(results)
}


# quickplot function
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