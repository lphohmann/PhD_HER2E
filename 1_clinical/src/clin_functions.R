# Function definitions for clinical data analyses

################################################################################
# functions
################################################################################

# function that creates KM plot for specified OM
KMplot <- function(OM,OMbin,OMstring,group.cohort.version,sdata) {
    
    # surv object
    data.surv <- Surv(OM, OMbin) 
    
    # fit
    fit <- survminer::surv_fit(data.surv~PAM50, data=sdata, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
        fit,
        censor.size = 6,
        censor.shape = "|",
        size = 3,
        risk.table = FALSE,       
        pval = TRUE,
        pval.size = 6,
        pval.coord = c(0,0.1),
        conf.int = FALSE,         
        xlim = c(0,max(OM[is.finite(OM)])),         
        xlab = paste(OMstring," (days)", sep = ""),
        ylab = paste(OMstring," event probability", sep = ""),
        ylim = c(0,1),
        palette = c("#d334eb", "#2176d5", "#34c6eb"), 
        legend = c(0.85,0.90),
        ggtheme = theme(legend.title = element_text(size=20), #20
                        legend.key.size = unit(0.5,"cm"), 
                        legend.text = element_text(size = 20), #20
                        axis.text.x = element_text(size = 20), #20
                        axis.title.x = element_text(size = 25), #25
                        axis.text.y = element_text(size = 20), #20
                        axis.title.y = element_text(size = 25),
                        plot.title = element_text(size=22)),
        title= paste(OMstring, ": ",group.cohort.version, sep = ""),
        legend.title = "Subtypes",
        legend.labs = c(paste("HER2E"," (",(table(sdata[!is.na(OM),]$PAM50)[1]),")",sep = ""),
                        paste("LUMA"," (",(table(sdata[!is.na(OM),]$PAM50)[2]),")",sep = ""),
                        paste("LUMB"," (",(table(sdata[!is.na(OM),]$PAM50)[3]),")",sep = "")),
        break.x.by = 500, # break X axis in time intervals of x (what is nicest here? maybe 365)
        break.y.by = 0.1)
    
    print(plot)
    
}
