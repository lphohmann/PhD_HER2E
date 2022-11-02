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

# DELETE THIS ONE AFTER CHANGING THE SCRIPT
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

### functions for singlegex script ###

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


# boxplot function for pam50 groups
uni_quickplot <- function(data, group.var, test.var, lumb.pos, lumb.sign, luma.pos, luma.sign, ylim, xlab="PAM50 subtype",
                           ylab, title) {
    plot <- ggplot(data, aes(x=as.factor(.data[[group.var]]),y=.data[[test.var]],fill=as.factor(.data[[group.var]]))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab(xlab) +
        ylab(ylab) +
        ylim(ylim) +
        ggtitle(title) +
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations=luma.sign, tip_length = 0.02, vjust=0.01, y_position = luma.pos, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations=lumb.sign, tip_length = 0.02, vjust=0.01, y_position = lumb.pos, size = 2, textsize = 15) + 
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

# boxplot real unviersal now
hp_uni_quickplot <- function(data, group.var, test.var, g1, g2, g3, g1.col = "#d334eb", g3.pos, g3.sign, g3.col = "#34c6eb", g2.pos, g2.sign, g2.col = "#2176d5", ylim, xlab="PAM50 subtype", ylab, title) {
    
    g1 <- ensym(g1)
    g2 <- ensym(g2)
    g3 <- ensym(g3)

    plot <- ggplot(data, aes(x=as.factor(.data[[group.var]]),y=.data[[test.var]],fill=as.factor(.data[[group.var]]))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab(xlab) +
        ylab(ylab) +
        ylim(ylim) +
        ggtitle(title) +
        geom_signif(comparisons=list(c(g1, g2)), annotations = g2.sign, tip_length = 0.02, vjust=0.01, y_position = g2.pos, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c(g1, g3)), annotations = g3.sign, tip_length = 0.02, vjust=0.01, y_position = g3.pos, size = 2, textsize = 15) + 
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              plot.title = element_text(size=30),
              legend.position = "none") +
        scale_fill_manual(values=c(!!rlang::ensym(g2) = g2.col, !!rlang::ensym(g3) = g3.col, !!rlang::ensym(g1) = g1.col)) +
        scale_x_discrete(limits = c(g2,g3,g1), labels = c(toupper(g2),toupper(g3),toupper(g1))) # check that these are in the correct order
    #scale_fill_manual(values=c(!!rlang::ensym(g2) = g2.col, !!rlang::ensym(g3) = g3.col, !!rlang::ensym(g1) = g1.col)) +
    #scale_x_discrete(limits = c(g2,g3,g1), labels = c(toupper(g2),toupper(g3),toupper(g1))) # check that these are in the correct order
    return(plot)
}
