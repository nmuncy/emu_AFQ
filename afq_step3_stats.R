

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")
library("ez")
library("dplyr")
library("lme4")



### Set Up
# Orienting paths - set globally
dataDir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
privateDir <- "/Users/nmuncy/Projects/emu_data/emu_private/"

plotDir_gam <- paste0(dataDir, "plots_gam/")
statsDir_gam <- paste0(dataDir, "stats_gam/")
plotDir_lm <- paste0(dataDir, "plots_lm/")
statsDir_lm <- paste0(dataDir, "stats_lm/")
tableDir <- paste0(dataDir, "tables/")


# set lists
groupType <- 1:3
tractList <- c("UNC_L", "UNC_R", "FA")


# Functions - general
func_df_master <- function(gType){
  
  ### --- Notes:
  #
  # This function will make a master dataframe that
  # contains node FA values, group membership, sex, PDS
  # and memory measures (LG/DI)
  #
  # Writes dataDir/Master_dataframe.csv
  
  # Get data
  df_afq <- read.delim(
    paste0(dataDir, "tract_profiles.csv"), sep = ",", header = T
  )
  colnames(df_afq) <- c("Counter", colnames(df_afq)[-1])
  
  df_full <- read.delim(
    paste0(privateDir, "emuR01_full_latest.csv"), sep = ",", header=T
  )
  
  df_pds <- read.delim(
    paste0(privateDir, "emuR01_pds_latest.csv"), sep = ",", header=T
  )
  
  df_adis <- read.delim(
    paste0(privateDir, "emuR01_adis.csv"), sep = ",", header=T
  )
  
  # make lists
  subjList <- unique(df_afq$subjectID)
  tractList <- unique(df_afq$tractID)
  nodeList <- unique(df_afq$nodeID)

  # add group (Adis), pars6, pds
  #   add LG/DI, sex, age
  df_afq$Sex <- df_afq$Age <- NA
  df_afq$Group <- df_afq$Pars6 <- df_afq$PDS <- NA
  df_afq$NegLDI <- df_afq$NeuLDI <- df_afq$NegLGI <- df_afq$NeuLGI <- NA
  
  for(subj in subjList){
    
    # determine subj indices
    ind_afq <- which(df_afq$subjectID == subj)
    ind_full <- which(df_full$src_subject_id == subj)
    ind_pds <- which(df_pds$emu_study_id == subj)
    ind_adis <- which(df_adis$Participant.ID == subj)
    
    # get pars
    h_anx <- df_full[ind_full,]$pars_6
    
    # Determine group in one of three ways:
    #   a) 0 = con, 1 = anx
    #       skip subj when not anx/phobia/control
    #   b) 0 = con, 1 = gad, 2 = social/separation
    #       1 = GAD in dx.1, or dx GAD but SAD not dx.1
    #   c) 0 = low, 1 = med, 2 = high pars
    #       low=0-3, med=4-12, high>12
    if(gType == 1){
      
      h_search <- c("Anxiety", "Phobia")
      
      if(sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis,])) != 0){
        h_group <- 1
      }else if(length(grep("None", df_adis[ind_adis,])) != 0){
        h_group <- 0
      }else{
        next
      }
      
    }else if(gType == 2){
      
      h_search <- c("Separation", "Social")
      
      if(
        grepl("Gen", df_adis[ind_adis,]$Diagnosis.1) == T |
        (sum(grep("Gen", df_adis[ind_adis,])) != 0 &
         sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis,])) == 0)
      ){
        h_group <- 1
      }else if(
        sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis,])) != 0
      ){
        h_group <- 2
      }else if(length(grep("None", df_adis[ind_adis,])) != 0){
        h_group <- 0
      }else{
        next
      }
      
    }else if(gType == 3){
      if(h_anx <= 3){
        h_group <- 0
      }else if(h_anx > 3 & h_anx < 13){
        h_group <- 1
      }else if(h_anx > 12){
        h_group <- 2
      }
    }
    
    # get pds
    h_pds <- df_pds[ind_pds,]$pds_shirtcliff
    
    # get age, sex
    h_age <- df_full[ind_full,]$pinf_age
    h_sex <- substr(df_full[ind_full,]$sex, 1, 1)
    if(h_sex == "f"){h_sexF <- 0}else if(h_sex == "m"){h_sexF <- 1}
    
    # get Beh counts
    neg_num_hit <- df_full[ind_full,]$negtarght_cnt_1WK
    neg_num_miss <- df_full[ind_full,]$negtargms_cnt_1WK
    neg_num_Lcr <- df_full[ind_full,]$neglurecr_cnt_1WK
    neg_num_Lfa <- df_full[ind_full,]$neglurefa_cnt_1WK
    neg_num_Fcr <- df_full[ind_full,]$negfoilcr_cnt_1WK
    neg_num_Ffa <- df_full[ind_full,]$negfoilfa_cnt_1WK
    
    neu_num_hit <- df_full[ind_full,]$neutarght_cnt_1WK
    neu_num_miss <- df_full[ind_full,]$neutargms_cnt_1WK
    neu_num_Lcr <- df_full[ind_full,]$neulurecr_cnt_1WK
    neu_num_Lfa <- df_full[ind_full,]$neulurefa_cnt_1WK
    neu_num_Fcr <- df_full[ind_full,]$neufoilcr_cnt_1WK
    neu_num_Ffa <- df_full[ind_full,]$neufoilfa_cnt_1WK
    
    # adjust 0 counts
    for(check in c("neg_num_hit",
                   "neg_num_miss",
                   "neg_num_Lcr",
                   "neg_num_Lfa",
                   "neg_num_Fcr",
                   "neg_num_Ffa",
                   "neu_num_hit",
                   "neu_num_miss",
                   "neu_num_Lcr",
                   "neu_num_Lfa",
                   "neu_num_Fcr",
                   "neu_num_Ffa")){
      check_val = get(check)
      if(check_val == 0){
        assign(check, 0.001)
      }
    }
    
    # Calculate Neg/Neu LD/GI
    #   LDI = p(N|L) - p(N|T)
    #   LGI = p(O|L) - p(O|F)
    neg_LDI <- round(
      (neg_num_Lcr / (neg_num_Lcr + neg_num_Lfa)) - 
        (neg_num_miss / (neg_num_miss + neg_num_hit)), 
      2)
    
    neu_LDI <- round(
      (neu_num_Lcr / (neu_num_Lcr + neu_num_Lfa)) - 
        (neu_num_miss / (neu_num_miss + neu_num_hit)), 
      2)
    
    neg_LGI <- round(
      (neg_num_Lfa / (neg_num_Lcr + neg_num_Lfa)) - 
        (neg_num_Ffa / (neg_num_Ffa + neg_num_Fcr)), 
      2)
    
    neu_LGI <- round(
      (neu_num_Lfa / (neu_num_Lcr + neu_num_Lfa)) - 
        (neu_num_Ffa / (neu_num_Ffa + neu_num_Fcr)), 
      2)
    
    # fill
    df_afq[ind_afq,]$Group <- h_group
    df_afq[ind_afq,]$Pars6 <- h_anx
    df_afq[ind_afq,]$PDS <- h_pds
    df_afq[ind_afq,]$Age <- h_age
    df_afq[ind_afq,]$Sex <- h_sexF
    
    df_afq[ind_afq,]$NegLDI <- neg_LDI
    df_afq[ind_afq,]$NeuLDI <- neu_LDI
    df_afq[ind_afq,]$NegLGI <- neg_LGI
    df_afq[ind_afq,]$NeuLGI <- neu_LGI
  }
  
  # clean NA (from Group skip), write csv
  df_out <- df_afq[complete.cases(df_afq$Group),]
  outFile <- paste0(dataDir, "Master_dataframe_G", gType,".csv")
  write.csv(df_out, file=outFile, quote=F, row.names = F)
  # return(df_out)
}

func_stat_mem <- function(df_afq){
  
  ### --- Notes:
  #
  # This function will conduct a
  # WBRM ANOVA on memory measures.
  
  # get data, make lists
  # df_afq <- read.csv(paste0(dataDir, "Master_dataframe.csv"))
  df_afq$Group <- factor(df_afq$Group)
  df_afq$Sex <- factor(df_afq$Sex)
  
  subjList <- unique(df_afq$subjectID)
  tractList <- unique(df_afq$tractID)
  nodeList <- unique(df_afq$nodeID)
  
  # set up LDI df
  df_LDI <-  as.data.frame(matrix(NA,nrow=2*length(subjList), ncol=4))
  colnames(df_LDI) <- c("Subj", "Group", "Meas", "Value")
  
  df_LDI$Subj <- rep(subjList, 2)
  df_LDI$Meas <- c(
    rep("NegLDI", length(subjList)), 
    rep("NeuLDI", length(subjList))
  )
  
  df_LDI$Group <- rep(
    df_afq[which(
      df_afq$nodeID == nodeList[1] & df_afq$tractID == tractList[1]
    ),]$Group, 
    2)
  
  df_LDI$Value <- c(
    df_afq[which(
      df_afq$nodeID == nodeList[1] & df_afq$tractID == tractList[1]
    ),]$NegLDI,
    df_afq[which(
      df_afq$nodeID == nodeList[1] & df_afq$tractID == tractList[1]
    ),]$NeuLDI
  )
  
  # plot
  boxplot(Value ~ Group*Meas, data=df_LDI)
  hist(df_LDI[which(df_LDI$Meas == "NegLDI"),]$Value)
  hist(df_LDI[which(df_LDI$Meas == "NeuLDI"),]$Value)
  
  # stats
  stats_LDI <- ezANOVA(df_LDI,
                       dv = Value,
                       wid = Subj,
                       between = Group,
                       within = Meas,
                       type = "III"
  )
  
  # write
  capture.output(stats_LDI, file = paste0(statsDir_gam, "Stats_AN_LDI.txt"))
  
  
  # repeat for LGI
  df_LGI <-  as.data.frame(matrix(NA,nrow=2*length(subjList), ncol=4))
  colnames(df_LGI) <- c("Subj", "Group", "Meas", "Value")
  
  df_LGI$Subj <- rep(subjList, 2)
  df_LGI$Meas <- c(
    rep("NegLGI", length(subjList)), 
    rep("NeuLGI", length(subjList))
  )
  
  df_LGI$Group <- rep(
    df_afq[which(
      df_afq$nodeID == nodeList[1] & df_afq$tractID == tractList[1]
    ),]$Group, 
    2)
  
  df_LGI$Value <- c(
    df_afq[which(
      df_afq$nodeID == nodeList[1] & df_afq$tractID == tractList[1]
    ),]$NegLGI,
    df_afq[which(
      df_afq$nodeID == nodeList[1] & df_afq$tractID == tractList[1]
    ),]$NeuLGI
  )
  
  boxplot(Value ~ Group*Meas, data=df_LGI)
  hist(df_LGI[which(df_LGI$Meas == "NegLGI"),]$Value)
  hist(df_LGI[which(df_LGI$Meas == "NeuLGI"),]$Value)
  
  stats_LGI <- ezANOVA(df_LGI,
                       dv = Value,
                       wid = Subj,
                       between = Group,
                       within = Meas,
                       type = "III"
  )
  
  capture.output(stats_LGI, file = paste0(statsDir_gam, "Stats_AN_LGI.txt"))
}

func_switch_g1 <- function(value){
  
  ### --- Notes:
  #
  # Switch for determining group, coloring
  #   for grouping set 1
  
  x_col <- switch(
    value,
    "0" = "blue",
    "1" = "darkred",
  )
  
  x_label <- switch(
    value,
    "0" = "Con",
    "1" = "Anx",
  )
  return(list(x_col, x_label))
}

func_switch_g2 <- function(value){
  
  ### --- Notes:
  #
  # Switch for determining group, coloring
  #   for grouping set 2
  
  x_col <- switch(
    value,
    "0" = "blue",
    "1" = "darkred",
    "2" = "black"
  )
  
  x_label <- switch(
    value,
    "0" = "Con",
    "1" = "GAD",
    "2" = "SAD"
  )
  return(list(x_col, x_label))
}

func_switch_g3 <- function(value){
  
  ### --- Notes:
  #
  # Switch for determining group, coloring
  #   for grouping set 3
  #
  # TODO maybe combine into func_switch_g3
  
  x_col <- switch(
    value,
    "0" = "blue",
    "1" = "darkred",
    "2" = "black"
  )
  
  x_label <- switch(
    value,
    "0" = "Low",
    "1" = "Med",
    "2" = "High"
  )
  return(list(x_col, x_label))
}

func_switch_name <- function(tract){
  
  ### --- Notes:
  #
  # Switch for determining name
  #   of tracts
  
  x_tract <- switch(
    tract,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "FA" = "A. Forceps",
  )
  return(x_tract)
}


# Functions - GAM
func_plot_gam <- function(model, tract, gType, df_tract){
  
  ### --- Notes:
  #
  # Will plot the GAM model of
  # dti data
  
  # plot
  df_pred <- predict.bam(
    model,
    exclude_terms = c("PDS", "Sex", "subjectID"),
    values=list(PDS = NULL, Sex = NULL),
    se.fit=T,
    type="response")
  
  df_pred <- data.frame(Group=df_tract$Group,
                        Sex=df_tract$Sex,
                        subjectID=df_tract$subjectID,
                        PDS=df_tract$PDS,
                        nodeID=df_tract$nodeID,
                        fit=df_pred$fit,
                        se.fit=df_pred$se.fit)
  
  h_tract <- func_switch_name(tract)
  h_title = paste0("GAM Fit of ", h_tract," FA Values")
  
  if(gType == 1){
    
    h_cols = c(func_switch_g1("0")[[1]][1], func_switch_g1("1")[[1]][1])
    names(h_cols) <- c("0", "1")
    h_breaks <- c("0", "1")
    h_labels <- c(func_switch_g1("0")[[2]][1], func_switch_g1("1")[[2]][1])
    
  }else if(gType == 2){

    h_cols = c(
      func_switch_g2("0")[[1]][1], 
      func_switch_g2("1")[[1]][1], 
      func_switch_g2("2")[[1]][1]
    )
    names(h_cols) <- c("0", "1", "2")
    h_breaks <- c("0", "1", "2")
    h_labels <- c(
      func_switch_g2("0")[[2]][1], 
      func_switch_g2("1")[[2]][1], 
      func_switch_g2("2")[[2]][1]
    )
    
  }else if(gType == 3){
    
    h_cols = c(
      func_switch_g3("0")[[1]][1], 
      func_switch_g3("1")[[1]][1], 
      func_switch_g3("2")[[1]][1]
    )
    names(h_cols) <- c("0", "1", "2")
    h_breaks <- c("0", "1", "2")
    h_labels <- c(
      func_switch_g3("0")[[2]][1], 
      func_switch_g3("1")[[2]][1], 
      func_switch_g3("2")[[2]][1]
    )
  }

  p <- ggplot(data = df_pred) +
    geom_smooth(mapping = aes(x=nodeID, y=fit, color=Group)) +
    ggtitle(h_title) +
    ylab("Fit FA") +
    xlab("Node ID")
  
  p + scale_color_manual(
    values = h_cols,
    breaks = h_breaks, 
    labels = h_labels
  )
  
  ggsave(paste0(plotDir_gam, "Plot_GAM_", tract, "_", "G", gType, ".png"))
}

func_stat_gam <- function(tract, df_tract, gType){
  
  ### --- Notes:
  #
  # This function has 3 main steps
  #
  # 1) Determine best GAM family.
  #   This is based off the distribution of the
  #   data. It was determined that a gamma
  #   fit better than beta for all tracts.
  #   It was also determined that k=40 worked
  #   for all tracts.
  #
  # 2) Model with covariates, compare models.
  #   All tracts had better model fit when PDS
  #   was added. The covariate model is then
  #   plotted.
  #
  # 3) A-B differences between splines are
  #   then plotted. The node which showed
  #   maximum difference is then used
  #   to make df_max.
  #
  # Note: this work can only be a function since
  #   it happened that all tracts had the same
  #   distribution, family, K, and best model.
  #   Individual 
  
  
  ### Model tract, no covariates
  
  # # plot mean data
  # ggplot(data = df_tract) +
  #   geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))
  # 
  # ggplot(data = df_tract) +
  #   geom_point(mapping = aes(x=nodeID, y=dti_fa, color=Group),size=0.3) +
  #   geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

  # # determine distribution
  # descdist(df_tract$dti_fa, discrete=F) # Could be beta or gamma
  # 
  # fit.beta <- fitdist(df_tract$dti_fa, "beta")
  # plot(fit.beta)
  # fit.beta$aic
  # 
  # fit.gamma <- fitdist(df_tract$dti_fa, "gamma")
  # plot(fit.gamma)
  # fit.gamma$aic
  
  #  determine k, compare families
  fit_gamma <- bam(dti_fa ~ Group +
                     Sex +
                     s(nodeID, by=Group, k=40) +
                     s(subjectID, bs="re"),
                   data = df_tract,
                   family = Gamma(link = "logit"),
                   method = "REML")
  
  # gam.check(fit_gamma, rep = 500)
  # 
  # fit_beta <- bam(dti_fa ~ Group +
  #                   Sex +
  #                   s(nodeID, by=Group, k=40) +
  #                   s(subjectID, bs="re"),
  #                 data = df_tract,
  #                 family = betar(link = "logit"),
  #                 method = "REML")
  # 
  # gam.check(fit_beta, rep = 500)
  # 
  # infoMessages('on')
  # compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  
  # get stats
  # summary(fit_gamma)
  capture.output(
    summary(fit_gamma), 
    file = paste0(
      statsDir_gam, "Stats_GAM-gamma_", tract, "_", "G", gType, ".txt"
    )
  )
  
  
  ### Model tract with covariates
  fit_cov_pds <- bam(dti_fa ~ Group +
                       Sex +
                       s(nodeID, by=Group, k=40) +
                       s(PDS, by=Sex) +
                       s(subjectID, bs="re"),
                     data = df_tract,
                     family = Gamma(link = "logit"),
                     method = "REML")
  
  # gam.check(fit_cov_pds, rep = 500)
  
  # Test cov model against gamma
  capture.output(
    compareML(fit_gamma, fit_cov_pds), 
    file = paste0(
      statsDir_gam, "Stats_GAM-comp_", tract, "_", "G", gType, ".txt"
    )
  )
  capture.output(
    summary(fit_cov_pds), 
    file = paste0(
      statsDir_gam, "Stats_GAM-cov_", tract, "_", "G", gType, ".txt"
    )
  )
  
  return(fit_cov_pds)
}

func_plot_diff1 <- function(model, tract, gType){
  
  ### --- Notes:
  #
  # This will draw plots and write tables of sig
  #   node differences for GAM when group style=1
  
  png(filename = paste0(
    plotDir_gam, "Plot_Diff_", tract, "_", "G", gType, ".png"), 
    width = 600, height = 600
  )
  
  gA <- func_switch_g1("0")[[2]][1]
  gB <- func_switch_g1("1")[[2]][1]
  
  par(mar=c(5,5,4,2))
  capture.output(plot_diff(model,
                           view="nodeID",
                           comp=list(Group=c("0", "1")),
                           rm.ranef = T,
                           main = paste0("Difference Scores, ", gA, "-", gB),
                           ylab = "Est. FA difference",
                           xlab = "Tract Node",
                           cex.lab = 2,
                           cex.axis = 2,
                           cex.main = 2.5,
                           cex.sub = 1.5),
                 file = paste0(tableDir, 
                               "Table_Diff_", 
                               tract, "_", "G",
                               gType, "_01.txt"
                 )
  )
  par(mar=c(5,4,4,2))
  dev.off()
}

func_plot_diff2 <- function(model, tract, gType){
  
  ### --- Notes:
  #
  # This will draw plots and write tables of sig
  #   node differences for GAM when group style=2
  
  png(filename = paste0(
    plotDir_gam, "Plot_Diff_", tract, "_", "G", gType, ".png"
  ), width = 1800, height = 600
  )
  
  if(gType == 2){
    gA <- func_switch_g2("0")[[2]][1]
    gB <- func_switch_g2("1")[[2]][1]
    gC <- func_switch_g2("2")[[2]][1]
  }else if(gType == 3){
    gA <- func_switch_g3("0")[[2]][1]
    gB <- func_switch_g3("1")[[2]][1]
    gC <- func_switch_g3("2")[[2]][1]
  }
  
  par(mfrow=c(1,3), mar=c(5,5,4,2))
  capture.output(plot_diff(model,
                           view="nodeID",
                           comp=list(Group=c("0", "1")),
                           rm.ranef = T,
                           main = paste0("Difference Scores, ", gA, "-", gB),
                           ylab = "Est. FA difference",
                           xlab = "Tract Node",
                           cex.lab = 2,
                           cex.axis = 2,
                           cex.main = 2.5,
                           cex.sub = 1.5),
                 file = paste0(tableDir, 
                               "Table_Diff_", 
                               tract, "_", "G", 
                               gType, "_01.txt"
                 )
  )
  
  par(mar=c(5,3,4,2))
  capture.output(plot_diff(model,
                           view="nodeID",
                           comp=list(Group=c("0", "2")),
                           rm.ranef = T,
                           main = paste0("Difference Scores, ", gA, "-", gC),
                           ylab = "",
                           xlab = "Tract Node",
                           cex.lab = 2,
                           cex.axis = 2,
                           cex.main = 2.5,
                           cex.sub = 2),
                 file = paste0(tableDir, 
                               "Table_Diff_", 
                               tract, "_", "G", 
                               gType, "_02.txt"
                 )
  )
  
  capture.output(plot_diff(model,
                           view="nodeID",
                           comp=list(Group=c("1", "2")),
                           rm.ranef = T,
                           main = paste0("Difference Scores, ", gB, "-", gC),
                           ylab = "",
                           xlab = "Tract Node",
                           cex.lab = 2,
                           cex.axis = 2,
                           cex.main = 2.5,
                           cex.sub = 2),
                 file = paste0(tableDir, 
                               "Table_Diff_", 
                               tract, "_", "G", 
                               gType, "_12.txt"
                 )
  )
  par(mfrow=c(1,1), mar=c(5,4,4,2))
  dev.off()
}

func_mkdf_diff1 <- function(model){
  
  ### --- Notes:
  #
  # Returns a dataframe of difference
  #   scores for each node.
  
  p01 <- plot_diff(model,
                   view="nodeID",
                   comp=list(Group=c("0", "1")),
                   rm.ranef = T,
                   plot = F)
  
  colnames(p01) <- c(colnames(p01[,1:4]), "Comp")
  p01$Comp <- "01"
  return(p01)
}

func_mkdf_diff2 <- function(model){
  
  ### --- Notes:
  #
  # Returns a dataframe of difference
  #   scores for each node.
  
  p01 <- plot_diff(model,
                   view="nodeID",
                   comp=list(Group=c("0", "1")),
                   rm.ranef = T,
                   plot = F)
  colnames(p01) <- c(colnames(p01[,1:4]), "Comp")
  p01$Comp <- "01"
  
  p02 <- plot_diff(model,
                   view="nodeID",
                   comp=list(Group=c("0", "2")),
                   rm.ranef = T,
                   plot = F)
  colnames(p02) <- c(colnames(p02[,1:4]), "Comp")
  p02$Comp <- "02"
  
  p12 <- plot_diff(model,
                   view="nodeID",
                   comp=list(Group=c("1", "2")),
                   rm.ranef = T, 
                   plot = F)
  colnames(p12) <- c(colnames(p12[,1:4]), "Comp")
  p12$Comp <- "12"
  
  df_out <- as.data.frame(matrix(NA, nrow=3*dim(p01)[1], ncol=dim(p01)[2]))
  colnames(df_out) <- colnames(p01)
  df_out <- rbind(p01, p02, p12)
  return(df_out)
}

func_stat_diff <- function(model, tract, gType){
  
  ### --- Notes:
  #
  # 1) Makes plots and tables
  #     Tables are not currently used,
  #     func_plot_diff and func_mkdf_diff could
  #     be combined.
  #
  # 2) Make dataframes of difference estimates
  #
  # 3) Determine nodes that differ in difference
  #     estimation. Repeats what is done in tables.
  #     Also gets max node.
  #
  # Returns list of nodes, max node
  
  # make plots and tables
  if(gType == 1){
    func_plot_diff1(model, tract, gType)
  }else if(gType == 2 | gType == 3){
    func_plot_diff2(model, tract, gType)
  }
  
  # # use tables to get sig nodes, save in df_sigNodes
  # df_sigNodes <- as.data.frame(matrix(NA, nrow=1, ncol=4))
  # colnames(df_sigNodes) <- c("Comparison", "Section", "Start", "End")
  # 
  # if(gType == 1){
  #   compList <- "01"
  # }else if(gType == 2 | gType == 3){
  #   compList <- c("01", "02", "12")
  # }
  # 
  # # get table made by func_plot_diff for e/comparison
  # for(comp in compList){
  #   
  #   # write, use bash
  #   h_cmd = paste0(
  #     "tail -n +10 ", 
  #     tableDir, "Table_Diff_", tract, "_", "G", gType, "_", comp, 
  #     ".txt | sed 's/-/,/g'"
  #   )
  #   h_lines <- system(h_cmd, intern = T)
  #   
  #   # make df
  #   h_df <- read.table(
  #     text=paste(h_lines, collapse = "\n"), 
  #     header = F, 
  #     stringsAsFactors = F, 
  #     sep = ","
  #   )
  #   
  #   # write to df_sigNodes
  #   for(i in 1:dim(h_df)[1]){
  #     df_sigNodes <- rbind(df_sigNodes, c(comp, i, h_df[i,1], h_df[i,2]))
  #   }
  # }
  # df_sigNodes <- na.omit(df_sigNodes)
  
  
  ## Step 2 - get plot_diff dataframes
  if(gType == 1){
    df_estDiff <- func_mkdf_diff1(model)
  }else if(gType == 2 | gType == 3){
    df_estDiff <- func_mkdf_diff2(model)
  }
  
  
  ## Step 3 - for grouping 2, determine
  #   nodes where all groups differ from e/o
  nodeList <- unique(df_estDiff$nodeID)
  diffList <- vector()
  for(node in nodeList){
    ind_node <- which(df_estDiff$nodeID == node)
    if(gType == 1){
      if((abs(df_estDiff[ind_node[1],]$est) - df_estDiff[ind_node[1],]$CI) > 0){
        diffList <- c(diffList, node)
      }
    }else if(gType == 2 | gType == 3){
      if(
        (abs(df_estDiff[ind_node[1],]$est) - df_estDiff[ind_node[1],]$CI) > 0 &
        (abs(df_estDiff[ind_node[2],]$est) - df_estDiff[ind_node[2],]$CI) > 0 &
        (abs(df_estDiff[ind_node[3],]$est) - df_estDiff[ind_node[3],]$CI) > 0
      ){
        diffList <- c(diffList, node)
      }
    }
  }
  
  if(length(diffList) == 0){
    return(NA)
    stop
  }
  
  # find max diff
  if(gType == 1){
    
    h_df <- subset(df_estDiff, nodeID %in% diffList)
    ind_max <- which(abs(h_df$est) == max(abs(h_df$est)))
    node_max <- h_df[ind_max,]$nodeID
    
  }else if(gType == 2 | gType == 3){
    f_sum <- function(x){
      sum(abs(df_estDiff[which(df_estDiff$nodeID == x),]$est))
    }
    h_sum <- sapply(diffList, f_sum)
    names(h_sum) <- diffList
    h_max <- which(h_sum == max(h_sum))
    node_max <- as.numeric(names(h_sum[h_max]))
  }
  
  return(list(diffList, node_max))
}  


# Functions - LM
func_mkdf_lm <- function(df_tract, node_list, gType, avg_max){
  
  ### --- Notes:
  #
  # Returns a dataframe for linear models
  #   with either average FA or max FA 
  #   difference
  
  subjList <- unique(df_tract$subjectID)
  df_out <- as.data.frame(matrix(NA, nrow=length(subjList), ncol=8))
  colnames(df_out) <- c("Subj", "FAvalue", "Pars6", "Group", 
                        "NeuLGI", "NeuLDI","NegLGI", "NegLDI")
  df_out$Subj <- subjList
  
  for(subj in subjList){
    
    ind_out <- which(df_out$Subj == subj)
    df_subj <- df_tract[which(df_tract$subjectID == subj),]
    
    df_out[ind_out,]$Pars6 <- df_subj[1,]$Pars6
    df_out[ind_out,]$Group <- as.numeric(df_subj[1,]$Group) - 1
    df_out[ind_out,]$NeuLGI <- df_subj[1,]$NeuLGI
    df_out[ind_out,]$NeuLDI <- df_subj[1,]$NeuLDI
    df_out[ind_out,]$NegLGI <- df_subj[1,]$NegLGI
    df_out[ind_out,]$NegLDI <- df_subj[1,]$NegLDI
    
    if(avg_max == "Avg"){
      df_node <- subset(df_subj, nodeID %in% node_list)
      df_out[ind_out,]$FAvalue <- round(mean(df_node$dti_fa), 4)
    }else if(avg_max == "Max"){
      ind_max <- which(df_subj$nodeID == node_list)
      df_out[ind_out,]$FAvalue <- round(df_subj[ind_max,]$dti_fa, 4)
    }
  }
  return(df_out)
}

func_plot_lm1 <- function(df_plot, avg_max, mem){
  
  ### --- Notes:
  #
  # Plots A and B
  
  # plot if group diff
  h_labels <- c(func_switch_g1("0")[[2]][1], func_switch_g1("1")[[2]][1])
  h_tract <- func_switch_name(tract)
  
  plot1.x <- df_plot[which(df_plot$Group == 0),]$FAvalue
  plot1.y <- df_plot[which(df_plot$Group == 0),]$MemScore
  
  plot2.x <- df_plot[which(df_plot$Group == 1),]$FAvalue
  plot2.y <- df_plot[which(df_plot$Group == 1),]$MemScore
  
  h_title <- paste(h_tract, "Spline Differences Predicting Memory Performance")
  x_title <- ifelse(avg_max == "Avg", "Mean FA", "Max FA")
  
  png(filename = paste0(
      plotDir_lm, "Plot_LM-", avg_max, "_", tract, "_G", gType, "_", mem, ".png"
    ),
    width = 8,
    height = 4,
    units = "in",
    res = 600
  )
  
  par(mfrow=c(1,2), 
      oma=c(0,0,2,0), 
      family="Times New Roman"
  )
  
  plot(plot1.x, plot1.y, 
       xlab = x_title, 
       ylab = mem, 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)),
       main = h_labels[1])
  abline(lm(plot1.y ~ plot1.x))
  
  plot(plot2.x, plot2.y, 
       xlab = x_title, 
       ylab = "",
       main = h_labels[2], 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)))
  abline(lm(plot2.y ~ plot2.x))
  
  mtext(h_title, outer = T, cex = 1.5)
  dev.off()
  par(mfrow=c(1,1))
}

func_plot_lm2 <- function(df_plot, avg_max, mem, gType){
  
  ### --- Notes:
  #
  # Plots A, B, and C
  
  # plot if group diff
  if(gType == 2){
    h_labels <- c(
      func_switch_g2("0")[[2]][1], 
      func_switch_g2("1")[[2]][1], 
      func_switch_g2("2")[[2]][1]
    )
  }else if(gType == 3){
    h_labels <- c(
      func_switch_g3("0")[[2]][1], 
      func_switch_g3("1")[[2]][1], 
      func_switch_g3("2")[[2]][1]
    )
  }
  
  h_tract <- func_switch_name(tract)
  
  plot1.x <- df_plot[which(df_plot$Group == 0),]$FAvalue
  plot1.y <- df_plot[which(df_plot$Group == 0),]$MemScore
  
  plot2.x <- df_plot[which(df_plot$Group == 1),]$FAvalue
  plot2.y <- df_plot[which(df_plot$Group == 1),]$MemScore
  
  plot3.x <- df_plot[which(df_plot$Group == 2),]$FAvalue
  plot3.y <- df_plot[which(df_plot$Group == 2),]$MemScore
  
  h_title <- paste(h_tract, "Spline Differences Predicting Memory Performance")
  x_title <- ifelse(avg_max == "Avg", "Mean FA", "Max FA")
  
  png(filename = paste0(
      plotDir_lm, "Plot_LM-", avg_max, "_", tract, "_G", gType, "_", mem, ".png"
    ),
    width = 8,
    height = 4,
    units = "in",
    res = 600
  )
  
  par(mfrow=c(1,3), 
      oma=c(0,0,2,0), 
      family="Times New Roman"
  )
  
  plot(plot1.x, plot1.y, 
       xlab = x_title, 
       ylab = mem, 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)),
       main = h_labels[1])
  abline(lm(plot1.y ~ plot1.x))
  
  plot(plot2.x, plot2.y, 
       xlab = x_title, 
       ylab = "",
       main = h_labels[2], 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)))
  abline(lm(plot2.y ~ plot2.x))
  
  plot(plot3.x, plot3.y, 
       xlab = x_title, 
       ylab = "",
       main = h_labels[3], 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)))
  abline(lm(plot3.y ~ plot3.x))
  
  mtext(h_title, outer = T, cex = 1.5)
  dev.off()
  par(mfrow=c(1,1))
}

func_stat_lm <- function(df_lm, tract, gType, avg_max){
  
  ### --- Notes:
  #
  # Conduct linear model for list of mem scores
  #   then make plots.
  
  memList <- c("NeuLGI", "NeuLDI","NegLGI", "NegLDI")
  
  for(mem in memList){
    
    df_mem <- as.data.frame(matrix(NA, nrow = dim(df_lm)[1], ncol=4))
    colnames(df_mem) <- c("Subj", "FAvalue", "Group", "MemScore")
    
    df_mem$Subj <- df_lm$Subj
    df_mem$FAvalue <- df_lm$FAvalue
    df_mem$Group <- df_lm$Group
    
    ind_mem <- grep(mem, colnames(df_lm))
    df_mem$MemScore <- df_lm[,ind_mem]
    
    fit.int <- lm(MemScore ~ FAvalue*Group, data = df_mem)
    capture.output(
      summary(fit.int),
      file = paste0(statsDir_lm,
                    "Stats_LM-", avg_max, "_", tract, "_G", gType, "_", mem, ".txt"
      )
    )
    capture.output(
      anova(fit.int),
      file = paste0(statsDir_lm,
                    "Stats_AN-", avg_max, "_", tract, "_G", gType, "_", mem, ".txt"
      )
    )
    
    if(gType == 1){
      func_plot_lm1(df_mem, avg_max, mem)
    }else if(gType == 2 | gType == 3){
      func_plot_lm2(df_mem, avg_max, mem, gType)
    }
  }
}



### --- Work:
#
# Two analyses (grouping types):
#   1) Con vs Anx
#   2) Con vs GAD vs SAD

for(gType in groupType){

  # make/get data, assign factor
  dataFile <- paste0(dataDir, "Master_dataframe_G", gType,".csv")
  
  if( ! file.exists(dataFile)){
    func_df_master(gType)
  }
  
  df_afq <- read.csv(dataFile)
  df_afq$Group <- factor(df_afq$Group)
  df_afq$Sex <- factor(df_afq$Sex)
  
  # # Check Memory behavior
  # func_stat_mem()
  
  for(tract in tractList){
    
    # subset df_afq for tract
    df_tract <- df_afq[which(df_afq$tractID == tract), ]
    df_tract$dti_fa <- round(df_tract$dti_fa, 3)
    
    # run gam, plot
    gamFile <- paste0(dataDir, "G", gType, "_gam_", tract, ".Rda")
    
    if( ! file.exists(gamFile)){
      h_gam <- func_stat_gam(tract, df_tract, gType)
      saveRDS(h_gam, file=gamFile)
      rm(h_gam)
    }
    
    gam_model <- readRDS(gamFile)
    func_plot_gam(gam_model, tract, gType, df_tract)
    
    # determine group differences
    #   Note: for group 2, 3 nodeList contains
    #     list where all groups differed from e/o
    nodeList <- func_stat_diff(gam_model, tract, gType)
    
    # deal w/no differences
    if(is.na(nodeList)){
      next
    }
    avg_nList <- nodeList[[1]]
    max_nList <- nodeList[[2]]
    
    # avg lm
    df_avg <- func_mkdf_lm(df_tract, avg_nList, gType, "Avg")
    func_stat_lm(df_avg, tract, gType, "Avg")
    
    # max lm
    df_max <- func_mkdf_lm(df_tract, max_nList, gType, "Max")
    func_stat_lm(df_max, tract, gType, "Max")
  }
}


# For Demographics
gTYpe <- 1
tract <- "UNC_L"

dataFile <- paste0(dataDir, "Master_dataframe_G", gType,".csv")
df_afq <- read.csv(dataFile)
df_subset <- df_afq[which(df_afq$nodeID == 0 & df_afq$tractID == tract),]

num_subj <- dim(df_subset)[1]
num_female <- length(which(df_subset$Sex == 0))
age_avg <- round(mean(df_subset$Age),2)
age_sd <- round(sd(df_subset$Age),2)



