# Function definitions for transcriptomic data analyses
#TODO just write one quickplot function and one ttest function with enough parameters to adapt them to all situations

################################################################################
# functions
################################################################################

# define "not in" operator 
'%!in%' <- function(x,y)!('%in%'(x,y))

### functions for metagenes script (LUMA,LUMB) ###

# function to calc. the score for each metagene for each sample
mgscore <- function(metagene,metagene.def,gex.data) {
    
    # extract the gex for each metagene
    mg.ids <- filter(metagene.def, module == metagene) %>% 
        pull(entrezgene_id)
    
    # filter out the corresponding gex data
    mg.gex.data <- gex.data %>% 
        rownames_to_column(var = "entrezgene_id") %>% 
        filter(entrezgene_id %in% mg.ids) %>% 
        column_to_rownames(var = "entrezgene_id") 
    
    # calc. the score for each sample
    result <- as.data.frame(apply(mg.gex.data, 2, mean)) %>% 
        dplyr::rename(!!metagene := "apply(mg.gex.data, 2, mean)")
    
    return(result)
}

# function to calc. the p value for each metagene
mgtest <- function(metagene.scores,anno) {
    
    # sample ids 
    Her2.cols <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
    LumA.cols <- anno %>% filter(PAM50=="LumA") %>% pull(sampleID)
    LumB.cols <- anno %>% filter(PAM50=="LumB") %>% pull(sampleID)
    
    # transpose
    tmg <- t(metagene.scores)
    
    #initialize storing vector
    H2vsLA.pval <- rep(0,nrow(tmg))
    H2vsLB.pval <- rep(0,nrow(tmg))
    
    # calc. pvalues
    for (i in 1:nrow(tmg)) {
        
        # group data
        her2.data <- as.numeric(tmg[i,Her2.cols])
        luma.data <- as.numeric(tmg[i,LumA.cols])
        lumb.data <- as.numeric(tmg[i,LumB.cols])
        
        # for Her2 vs LumA
        # equal variance check
        if (var.test(unlist(her2.data),unlist(luma.data), alternative = "two.sided")$p.value <= 0.05) {
            H2vsLA.res <- t.test(her2.data,luma.data, var.equal = FALSE)
        } else {
            H2vsLA.res <- t.test(her2.data,luma.data, var.equal = TRUE)
        }
        # save results
        H2vsLA.pval[i] <- H2vsLA.res$p.value
        
        # for Her2 vs LumB
        # equal variance check
        if (var.test(unlist(her2.data),unlist(lumb.data), alternative = "two.sided")$p.value <= 0.05) {
            H2vsLB.res <- t.test(her2.data,lumb.data, var.equal = FALSE)
        } else {
            H2vsLB.res <- t.test(her2.data,lumb.data, var.equal = TRUE)
        }
        # save results
        H2vsLB.pval[i] <- H2vsLB.res$p.value 
    
    }
    
    # add the results to the dataframe
    results <- as.data.frame(tmg) %>% 
        add_column(H2vsLA_pval = H2vsLA.pval,H2vsLB_pval = H2vsLB.pval) %>% 
        dplyr::select(H2vsLA_pval,H2vsLB_pval)
    
    return(results)
}


# quickplot function
quickplot <- function(mg.anno, metagene, lumb.sig, lumb.pos, luma.sig, luma.pos, ylim) {
    plot <- ggplot(mg.anno, aes(x=as.factor(PAM50),y=.data[[metagene]],fill=as.factor(PAM50))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("PAM50 subtype") +
        ylab("Metagene score") +
        ylim(ylim) +
        ggtitle(paste(metagene," metagene scores in PAM50 subtypes (ERpHER2n)",sep="")) +
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations=luma.sig, tip_length = 0.02, vjust=0.01, y_position = luma.pos, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations=lumb.sig, tip_length = 0.02, vjust=0.01, y_position = lumb.pos, size = 2, textsize = 15) + 
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              plot.title = element_text(size=30),
              legend.position = "none") +
        scale_fill_manual(values=c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb")) +
        scale_x_discrete(limits = c("LumA","LumB","Her2"), labels = c("LUMA","LUMB","HER2E")) # check that these are in the correct order
    return(plot)
}


### functions for metagenes script (HER2pHER2E,HER2p) ###

# quickplot function for her2p samples
her2p_quickplot <- function(mg.anno, metagene, pher2e.sig, pher2e.pos, nonher2e.sig, nonher2e.pos, ylim) { 
    plot <- ggplot(mg.anno, aes(x=as.factor(Group),y=.data[[metagene]],fill=as.factor(Group))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("Subtype") +
        ylab("Metagene score") +
        ggtitle(paste(
            metagene," metagene scores in HER2 subtypes (ERp)",sep="")) +
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              plot.title = element_text(size=30),
              legend.position = "none") +
        scale_fill_manual(values=c(HER2p_nonHER2E = "#d8b365", HER2p_HER2E = "#5ab4ac", HER2n_HER2E ="#d334eb")) +
        scale_x_discrete(limits = levels(mg.anno$Group)) +
        ylim(ylim) +
        geom_signif(comparisons=list(c("HER2n_HER2E", "HER2p_nonHER2E")), annotations=nonher2e.sig, tip_length = 0.02, vjust=0.01, y_position = nonher2e.pos, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("HER2n_HER2E", "HER2p_HER2E")), annotations=pher2e.sig, tip_length = 0.02, vjust=0.01, y_position = pher2e.pos, size = 2, textsize = 15) 
        
    return(plot)
}

# function to calc. the p value for each metagene in her2p,her2pher2e,her2e
her2p_mgtest <- function(metagene.scores,anno) {
    
    # sample ids 
    HER2nHER2E.cols <- anno %>% filter(Group=="HER2n_HER2E") %>% pull(sampleID)
    HER2pHER2E.cols <- anno %>% filter(Group=="HER2p_HER2E") %>% pull(sampleID)
    HER2p.cols <- anno %>% filter(Group=="HER2p_nonHER2E") %>% pull(sampleID)
    
    # transpose
    tmg <- t(metagene.scores)
    
    #initialize storing vector
    HER2nHER2E.vs.HER2pHER2E.pval <- rep(0,nrow(tmg))
    HER2nHER2E.vs.HER2p.pval <- rep(0,nrow(tmg))
    
    # calc. pvalues
    for (i in 1:nrow(tmg)) {
        
        # group data
        HER2nHER2E.data <- as.numeric(tmg[i,HER2nHER2E.cols])
        HER2pHER2E.data <- as.numeric(tmg[i,HER2pHER2E.cols])
        HER2p.data <- as.numeric(tmg[i,HER2p.cols])
        
        # for HER2nHER2E vs HER2pHER2E
        # equal variance check
        if (var.test(unlist(HER2nHER2E.data),unlist(HER2pHER2E.data), alternative = "two.sided")$p.value <= 0.05) {
            HER2nHER2E.vs.HER2pHER2E.res <- t.test(HER2nHER2E.data,HER2pHER2E.data, var.equal = FALSE)
        } else {
            HER2nHER2E.vs.HER2pHER2E.res <- t.test(HER2nHER2E.data,HER2pHER2E.data, var.equal = TRUE)
        }
        # save results
        HER2nHER2E.vs.HER2pHER2E.pval[i] <- HER2nHER2E.vs.HER2pHER2E.res$p.value
        
        # for HER2nHER2E vs HER2p
        # equal variance check
        if (var.test(unlist(HER2nHER2E.data),unlist(HER2p.data), alternative = "two.sided")$p.value <= 0.05) {
            HER2nHER2E.vs.HER2p.res <- t.test(HER2nHER2E.data,HER2p.data, var.equal = FALSE)
        } else {
            HER2nHER2E.vs.HER2p.res <- t.test(HER2nHER2E.data,HER2p.data, var.equal = TRUE)
        }
        # save results
        HER2nHER2E.vs.HER2p.pval[i] <- HER2nHER2E.vs.HER2p.res$p.value
        
    }
    
    # add the results to the dataframe
    results <- as.data.frame(tmg) %>% 
        add_column(HER2p_HER2E_pval = HER2nHER2E.vs.HER2pHER2E.pval,HER2p_nonHER2E_pval = HER2nHER2E.vs.HER2p.pval) %>% 
        dplyr::select(HER2p_HER2E_pval,HER2p_nonHER2E_pval)
    
    return(results)
}


### functions for singlegex script ###
# delete?
## do this later for the extracted gex data

# get gex, sampleid, pam50
get_gex <- function(id,gex.data,anno) {
    # extract the gex for each gene
    single.gex <- gex.data %>% 
        select_if(~ !any(is.na(.))) %>% 
        rownames_to_column(var = "gene_id") %>% 
        filter(gene_id == id) %>% 
        column_to_rownames(var = "gene_id") %>% 
        t(.) %>% as.data.frame() %>% 
        rownames_to_column(var="sampleID")
    # combine with anno
    single.gex <- merge(single.gex, anno[,c("PAM50","sampleID")],by="sampleID") %>% relocate(PAM50,.after=sampleID)
}

# summary stats.
get_stats <- function(data,grouping_var,stat_var) {
    res <- data %>% 
        group_by(!!rlang::ensym(grouping_var)) %>%
        get_summary_stats(!!rlang::ensym(stat_var), type = "mean_sd")
    return(res)
}

# ttest function with equal var check
var_ttest <- function(dat1,dat2) {
    # equal var check
    if (var.test(dat1,dat2)$p.value <= 0.05) {
        res <- t.test(dat1,dat2, var.equal = FALSE)
    } else {
        res <- t.test(dat1,dat2, var.equal = TRUE)
    }
    return(res)
}

# function for pairwise ttests (g1vsg2,g1vsg3)
pair_ttest <- function(data,anno=NULL,group.var,test.var,g1,g2,g3) {
    
    # if anno data was provided it means that the data does not contain the sample pam50 annotations and is assumed to have the structure(rownames=sampleID; columns=data). The anno object does have to have the structure of sampleID as column, PAM50 (grouping var) as column
    if(!is.null(anno)) { 
        data <- merge(data %>% rownames_to_column(var="sampleID"),anno[,c("sampleID",group.var)],by="sampleID")
    }
    
    # group data
    dat1 <- data %>% 
        filter(!!rlang::ensym(group.var)==g1) %>% 
        pull(!!rlang::ensym(test.var))
    
    dat2 <- data %>% 
        filter(!!rlang::ensym(group.var)==g2) %>% 
        pull(!!rlang::ensym(test.var))
    dat3 <- data %>% 
        filter(!!rlang::ensym(group.var)==g3) %>% 
        pull(!!rlang::ensym(test.var))
    
    # test and output formatting
    var_pair <- c(paste(g1,".",g2,sep=""),paste(g1,".",g3,sep=""))
    pval <- c(var_ttest(dat1,dat2)$p.value, var_ttest(dat1,dat3)$p.value)
    signif <- c(symnum(pval[1], corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                       symbols = c("****", "***", "**", "*", "ns")),
                symnum(pval[2], corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                       symbols = c("****", "***", "**", "*", "ns")))
    
    return(data.frame(var_pair,pval,signif))
}


# universal boxplot function
uni_quickplot <- function(data, group.var, test.var, pair_ttest.res, lumb.pos, luma.pos, ylim, xlab="PAM50 subtype",
                           ylab, title) {
    plot <- ggplot(data, aes(x=as.factor(.data[[group.var]]),y=.data[[test.var]],fill=as.factor(.data[[group.var]]))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab(xlab) +
        ylab(ylab) +
        ylim(ylim) +
        ggtitle(title) +
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations=pair_ttest.res[1,3], tip_length = 0.02, vjust=0.01, y_position = luma.pos, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations=pair_ttest.res[2,3], tip_length = 0.02, vjust=0.01, y_position = lumb.pos, size = 2, textsize = 15) + 
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              plot.title = element_text(size=30),
              legend.position = "none") +
        scale_fill_manual(values=c(LumA = "#2176d5", LumB = "#34c6eb", Her2 ="#d334eb")) +
        scale_x_discrete(limits = c("LumA","LumB","Her2"), labels = c("LUMA","LUMB","HER2E")) # check that these are in the correct order
    return(plot)
}
