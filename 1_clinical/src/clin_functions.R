# Function definitions for clinical data analyses

################################################################################
# functions
################################################################################

# loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
}

# function that creates KM plot for specified OM
KMplot <- function(OM,OMbin,OMstring,group.cohort.version,sdata) {
    
    # surv object
    data.surv <- Surv(OM, OMbin) 
    
    # fit
    fit <- survminer::surv_fit(data.surv~PAM50, data=sdata, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
        fit,
        censor.size = 8,
        censor.shape = "|",
        size = 5,
        risk.table = FALSE,       
        pval = TRUE,
        pval.size = 8,
        pval.coord = c(0,0.1),
        conf.int = FALSE,         
        xlim = c(0,max(OM[is.finite(OM)])),         
        xlab = paste(OMstring," (years)", sep = ""),
        ylab = paste(OMstring," event probability", sep = ""), # ggf just label as "event probability"
        ylim = c(0,1),
        palette = c("#2176d5", "#d334eb", "#34c6eb"), 
        legend = c(0.9,0.96),
        ggtheme = theme(legend.title = element_text(size=25), #20
                        legend.key.size = unit(0.5,"cm"), 
                        legend.text = element_text(size = 25), #20
                        axis.text.x = element_text(size = 25), #20
                        axis.title.x = element_text(size = 30), #25
                        axis.text.y = element_text(size = 25), #20
                        axis.title.y = element_text(size = 30),
                        plot.title = element_text(size=30)),
        title= paste(OMstring, ": ",group.cohort.version, sep = ""),
        legend.title = "Subtypes",
        legend.labs = c(paste("LUMA"," (",table(sdata[!is.na(OM),]$PAM50)[1],")",sep = ""),
                        paste("HER2E"," (",table(sdata[!is.na(OM),]$PAM50)[2],")",sep = ""),
                        paste("LUMB"," (",table(sdata[!is.na(OM),]$PAM50)[3],")",sep = "")),
        break.x.by = 1, # break X axis in time intervals of x (what is nicest here? maybe 365)
        break.y.by = 0.1)
    
    return(plot)
    
}

# function that created univariate Cox proportional hazards model forest plot
unicox <- function(data,surv,title) {
    
    # Model construction
    main.pam50 <- coxph(surv~PAM50, data=data)
    
    # Result
    print(summary(main.pam50))

    # forest 
    ggforest(main.pam50,fontsize = 2,data=data,main=title)
}

# function that created multivariate Cox proportional hazards model forest plot
mvcox <- function(data,surv,title) {
    
    # Model construction 
    # parameters to incl: PAM50, Age, Grade, TumSize
    main.all <- coxph(surv~PAM50+Age+LN+TumSize+Grade, data=data) 
    
    # Result
    print(summary(main.all))
    
    # Plot forest 
    ggforest(main.all,fontsize = 2,cpositions = c(0.01,0.13,0.35),data=data,main=title)
}

# function that creates KM plot for HER2p for specified OM
HER2p_KMplot <- function(OM,OMbin,OMstring,group.cohort.version,sdata) {
    
    # surv object
    data.surv <- Surv(OM, OMbin) 
    
    # fit
    fit <- survminer::surv_fit(data.surv~PAM50, data=sdata, conf.type="log-log") # weird bug: survival::survfit() cant be passed data in function call ?! so i use survminer::surv_fit()
    #survdiff(data.surv ~ PAM50, data = sdata) 
    
    plot <- ggsurvplot(
        fit,
        censor.size = 8,
        censor.shape = "|",
        size = 5,
        risk.table = FALSE,       
        pval = TRUE,
        pval.size = 8,
        pval.coord = c(0,0.1),
        conf.int = FALSE,         
        xlim = c(0,max(OM[is.finite(OM)])),         
        xlab = paste(OMstring," (days)", sep = ""),
        ylab = paste(OMstring," event probability", sep = ""), # ggf just label as "event probability"
        ylim = c(0,1),
        palette = c("#d334eb", "#2176d5"), 
        legend = c(0.9,0.96),
        ggtheme = theme(legend.title = element_text(size=25), #20
                        legend.key.size = unit(0.5,"cm"), 
                        legend.text = element_text(size = 25), #20
                        axis.text.x = element_text(size = 25), #20
                        axis.title.x = element_text(size = 30), #25
                        axis.text.y = element_text(size = 25), #20
                        axis.title.y = element_text(size = 30),
                        plot.title = element_text(size=30)),
        title= paste(OMstring, ": ",group.cohort.version, sep = ""),
        legend.title = "Subtypes",
        legend.labs = c(paste("HER2E"," (",table(sdata[!is.na(OM),]$PAM50)[1],")",sep = ""),
                        paste("nonHER2E"," (",table(sdata[!is.na(OM),]$PAM50)[2],")",sep = "")),
        break.x.by = 500, # break X axis in time intervals of x (what is nicest here? maybe 365)
        break.y.by = 0.1)
    
    print(plot)
    }