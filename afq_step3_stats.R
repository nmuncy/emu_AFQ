

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")
library("ez")
library("dplyr")
library("lme4")


### TODO:
#
# 1) Combine plot_diff functions
#
# 2) find group diffs from plot_diff when
#     group n=3



### Set Up
# Orienting paths - set globally
dataDir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
privateDir <- "/Users/nmuncy/Projects/emu_private/"
plotDir_gam <- paste0(dataDir, "plots_gam/")
statsDir_gam <- paste0(dataDir, "stats_gam/")
plotDir_lm <- paste0(dataDir, "plots_lm/")
statsDir_lm <- paste0(dataDir, "stats_lm/")
tableDir <- paste0(dataDir, "tables/")


# set lists
groupType <- 1:2
tractList <- c("UNC_L", "UNC_R", "FA")


# Functions
func_df_master <- function(g_type){
  
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
    
    # add pars
    h_anx <- df_full[ind_full,]$pars_6
    
    # determine group in one of two ways:
    #   a) 0 = con, 1 = anx
    #       skip subj when not anx/phobia/control
    #   b) 0 = con, 1 = gad, 2 = social/separation
    #       1 = GAD in dx.1, or dx GAD but SAD not dx.1
    #
    # ugly - maybe convert to case switch?
    
    if(g_type == 1){
      
      h_search <- c("Anxiety", "Phobia")
      
      if(sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis,])) != 0){
        h_group <- 1
      }else if(length(grep("None", df_adis[ind_adis,])) != 0){
        h_group <- 0
      }else{
        next
      }
      
    }else if(g_type == 2){
      
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
  outFile <- paste0(dataDir, "Master_dataframe_G", g_type,".csv")
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

func_plot_gam <- function(model, tract, g_type, df_tract){
  
  ### --- Notes:
  #
  # Will plot the GAM model of
  # dti data
  #
  # wrapped by func_stat_gam
  
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
  
  if(g_type == 1){
    h_cols <- c("0" = "blue", "1" = "darkred")
    h_breaks <- c("0", "1")
    h_labels <- c("Con", "Anx")
  }else if(g_type == 2){
    h_cols <- c("0" = "blue", "1" = "darkred", "2" = "black")
    h_breaks <- c("0", "1", "2")
    h_labels <- c("Con", "GAD", "SAD")
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
  
  ggsave(paste0(plotDir_gam, "Plot_GAM_", tract, "_", "G", g_type, ".png"))
}

func_plot_diff <- function(model, tract, g_type){
  
  ### --- Notes:
  #
  # This function will determine, plot differences
  # between two splines.
  #
  # A table of significant nodes is written.
  #
  # wrapped by func_stat_diff
  
  if(g_type == 2){
    
    png(filename = paste0(
        plotDir_gam, "Plot_Diff_", tract, "_", "G", g_type, ".png"
      ), width = 1800, height = 600
    )
    par(mfrow=c(1,3))
    par(mar=c(5,5,4,2))
    
    capture.output(plot_diff(model,
                       view="nodeID",
                       comp=list(Group=c("0", "1")),
                       rm.ranef = T,
                       main = "Difference Scores, Con-GAD",
                       ylab = "Est. FA difference",
                       xlab = "Tract Node",
                       cex.lab = 2,
                       cex.axis = 2,
                       cex.main = 2.5,
                       cex.sub = 1.5),
                   file = paste0(tableDir, 
                                 "Table_Diff_", 
                                 tract, "_", "G", 
                                 g_type, "_01.txt"
                                 )
                   )
    
    par(mar=c(5,3,4,2))
    
    capture.output(plot_diff(model,
                       view="nodeID",
                       comp=list(Group=c("0", "2")),
                       rm.ranef = T,
                       main = "Difference Scores, Con-SAD",
                       ylab = "",
                       xlab = "Tract Node",
                       cex.lab = 2,
                       cex.axis = 2,
                       cex.main = 2.5,
                       cex.sub = 2),
                   file = paste0(tableDir, 
                                 "Table_Diff_", 
                                 tract, "_", "G", 
                                 g_type, "_02.txt"
                                 )
                   )
    
    capture.output(plot_diff(model,
                             view="nodeID",
                             comp=list(Group=c("1", "2")),
                             rm.ranef = T,
                             main = "Difference Scores, GAD-SAD",
                             ylab = "",
                             xlab = "Tract Node",
                             cex.lab = 2,
                             cex.axis = 2,
                             cex.main = 2.5,
                             cex.sub = 2),
                   file = paste0(tableDir, 
                                 "Table_Diff_", 
                                 tract, "_", "G", 
                                 g_type, "_12.txt"
                                 )
                   )
    
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2))
    dev.off()
    
  }else if(g_type == 1){
    
    png(filename = paste0(
      plotDir_gam, "Plot_Diff_", tract, "_", "G", g_type, ".png"), 
      width = 600, height = 600
    )
    par(mar=c(5,5,4,2))
    capture.output(plot_diff(model,
                             view="nodeID",
                             comp=list(Group=c("0", "1")),
                             rm.ranef = T,
                             main = "Difference Scores, Con-Anx",
                             ylab = "Est. FA difference",
                             xlab = "Tract Node",
                             cex.lab = 2,
                             cex.axis = 2,
                             cex.main = 2.5,
                             cex.sub = 1.5),
                   file = paste0(tableDir, 
                                 "Table_Diff_", 
                                 tract, "_", "G",
                                 g_type, "_01.txt"
                                 )
                   )
    par(mar=c(5,4,4,2))
    dev.off()
  }
  
  # return(list(p01,p02,p12))
}

func_max_diff <- function(model, g_type){
  
  ### --- Notes:
  #
  # Return node of max diff per contrast
  
  if(g_type == 2){
    
    p01 <- plot_diff(model,
                     view="nodeID",
                     comp=list(Group=c("0", "1")),
                     rm.ranef = T,
                     plot = F)
    
    m01 <- p01[which(abs(p01$est) == max(abs(p01$est))),]$nodeID
    
    p02 <- plot_diff(model,
                     view="nodeID",
                     comp=list(Group=c("0", "2")),
                     rm.ranef = T,
                     plot = F)
    
    m02 <- p02[which(abs(p02$est) == max(abs(p02$est))),]$nodeID
    
    p12 <- plot_diff(model,
                     view="nodeID",
                     comp=list(Group=c("1", "2")),
                     rm.ranef = T, 
                     plot = F)
    
    m12 <- p12[which(abs(p12$est) == max(abs(p12$est))),]$nodeID
    
    df_out <- as.data.frame(matrix(NA, nrow=3, ncol=2))
    colnames(df_out) <- c("Comparison", "Node")
    df_out[,1] <- c("01", "02", "12")
    df_out[,2] <- c(m01, m02, m12)
    
  }else if(g_type == 1){
    
    p01 <- plot_diff(model,
                     view="nodeID",
                     comp=list(Group=c("0", "1")),
                     rm.ranef = T,
                     plot = F)
    
    m01 <- p01[which(abs(p01$est) == max(abs(p01$est))),]$nodeID
    
    df_out <- as.data.frame(matrix(NA, nrow=1, ncol=2))
    colnames(df_out) <- c("Comparison", "Node")
    df_out[,1] <- "01"
    df_out[,2] <- m01
  }
  
  return(df_out)
}

func_stat_gam <- function(tract, df_tract, g_type){
  
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
      statsDir_gam, "Stats_GAM-gamma_", tract, "_", "G", g_type, ".txt"
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
      statsDir_gam, "Stats_GAM-comp_", tract, "_", "G", g_type, ".txt"
    )
  )
  capture.output(
    summary(fit_cov_pds), 
    file = paste0(
      statsDir_gam, "Stats_GAM-cov_", tract, "_", "G", g_type, ".txt"
    )
  )
  
  # # plot
  # df_pred <- predict.bam(
  #   fit_cov_pds,
  #   exclude_terms = c("PDS", "Sex","subjectID"),
  #   values=list(PDS = NULL, Sex = NULL),
  #   se.fit=T,
  #   type="response")
  # 
  # df_pred <- data.frame(Group=df_tract$Group,
  #                       Sex=df_tract$Sex,
  #                       subjectID=df_tract$subjectID,
  #                       PDS=df_tract$PDS,
  #                       nodeID=df_tract$nodeID,
  #                       fit=df_pred$fit,
  #                       se.fit=df_pred$se.fit)
  # 
  # h_tract <- switch(
  #   tract,
  #   "UNC_L" = "L. Uncinate",
  #   "FA" = "A. Forceps"
  # )
  # 
  # plot_title = paste0("GAM Fit of ", h_tract," FA Values")
  # func_plot_gam(df_pred, plot_title, tract, g_type)
  return(fit_cov_pds)
}
  
func_stat_diff <- function(model, tract, g_type){
  
  ### --- Notes:
  #
  # Test for differences in GAM splines

  func_plot_diff(model, tract, g_type)

  # make table of sig regions
  df_out <- as.data.frame(matrix(NA, nrow=1, ncol=4))
  colnames(df_out) <- c("Comparison", "Section", "Start", "End")
  if(g_type == 1){
    compList <- "01"
  }else if(g_type == 2){
    compList <- c("01", "02", "12")
  }
  for(comp in compList){
    
    # get table made by func_plot_diff
    h_cmd = paste0(
      "tail -n +10 ", 
      tableDir, "Table_Diff_", tract, "_", "G", g_type, "_", comp, 
      ".txt | sed 's/-/,/g'"
    )
    
    h_lines <- system(h_cmd, intern = T)
    h_df <- read.table(
      text=paste(h_lines, collapse = "\n"), 
      header = F, 
      stringsAsFactors = F, 
      sep = ","
    )
    
    for(i in 1:dim(h_df)[1]){
      df_out <- rbind(df_out, c(comp, i, h_df[i,1], h_df[i,2]))
    }
  }
  df_out <- na.omit(df_out)
  return(df_out)
}

func_df_avg <- function(comp, df_tract, df_diff){
  
  ### --- Notes:
  # 
  # Make a dataframe of averaged FA
  #   values from all regions that differ
  #   between two splines
  
  grpA <- as.numeric(substr(comp, start=1, stop=1))
  grpB <- as.numeric(substr(comp, start=2, stop=2))
  df_comp <- df_diff[which(df_diff$Comparison == comp),]
  
  # make dataframe
  subjList <- unique(
    df_tract[which(
      df_tract$Group == grpA | df_tract$Group == grpB
    ),]$subjectID
  )
  
  df_lm <- as.data.frame(matrix(NA, nrow=length(subjList), ncol=9))
  colnames(df_lm) <- c("Comp", "Subj", "FAvalue", 
                       "Pars6", "Group", 
                       "NeuLGI", "NeuLDI", 
                       "NegLGI", "NegLDI")
  df_lm$Subj <- subjList
  
  for(subj in subjList){
    h_mean <- vector()
    for(i in 1:dim(df_comp)[1]){
      
      h_start <- which(
        df_tract$subjectID == subj & 
        df_tract$nodeID == df_comp[i,]$Start
      )
      
      h_end <- which(
        df_tract$subjectID == subj & 
        df_tract$nodeID == df_comp[i,]$End
      )
      
      h_mean <- c(h_mean, mean(df_tract[h_start:h_end,]$dti_fa))
    }
    
    ind_out <- which(df_lm$Subj == subj)
    df_lm[ind_out,]$FAvalue <- round(mean(h_mean), 4)
    
    ind_subj <- which(df_tract$subjectID == subj)[1]
    df_lm[ind_out,]$Pars6 <- df_tract[ind_subj,]$Pars6
    df_lm[ind_out,]$Group <- as.numeric(df_tract[ind_subj,]$Group) - 1
    
    df_lm[ind_out,]$NeuLGI <- df_tract[ind_subj,]$NeuLGI
    df_lm[ind_out,]$NeuLDI <- df_tract[ind_subj,]$NeuLDI
    df_lm[ind_out,]$NegLGI <- df_tract[ind_subj,]$NegLGI
    df_lm[ind_out,]$NegLDI <- df_tract[ind_subj,]$NegLDI
    
    df_lm[ind_out,]$Comp <- comp
  }
  df_lm$Group <- factor(df_lm$Group)
  return(df_lm)
}

func_df_max <- function(comp, df_tract, df_max){
  
  ### --- Notes:
  #
  # make a dataframe of FA values
  #   from single region showing showing
  #   maximal A-B difference
  
  grpA <- as.numeric(substr(comp, start=1, stop=1))
  grpB <- as.numeric(substr(comp, start=2, stop=2))
  node <- df_max[which(df_max$Comparison == comp),]$Node
  
  # make dataframe
  subjList <- unique(
    df_tract[which(
      df_tract$Group == grpA | df_tract$Group == grpB
    ),]$subjectID
  )
  df_lm <- as.data.frame(matrix(NA, nrow=length(subjList), ncol=9))
  colnames(df_lm) <- c("Comp","Subj", "FAvalue", 
                       "Pars6", "Group", 
                       "NeuLGI", "NeuLDI", 
                       "NegLGI", "NegLDI")
  df_lm$Subj <- subjList
  
  for(subj in subjList){
    
    ind_data <- which(df_tract$subjectID == subj & df_tract$nodeID == node)
    ind_out <- which(df_lm$Subj == subj)
    
    df_lm[ind_out,]$Comp <- comp
    df_lm[ind_out,]$FAvalue <- df_tract[ind_data,]$dti_fa
    df_lm[ind_out,]$Pars6 <- df_tract[ind_data,]$Pars6
    df_lm[ind_out,]$Group <- as.numeric(df_tract[ind_data,]$Group)-1
    df_lm[ind_out,]$NeuLGI <- df_tract[ind_data,]$NeuLGI
    df_lm[ind_out,]$NeuLDI <- df_tract[ind_data,]$NeuLDI
    df_lm[ind_out,]$NegLGI <- df_tract[ind_data,]$NegLGI
    df_lm[ind_out,]$NegLDI <- df_tract[ind_data,]$NegLDI
  }
  df_lm$Group <- factor(df_lm$Group)
  return(df_lm)
}

func_stat_lm <- function(df_lm, tract, gType, h_str, comp){
  
  ### --- Notes:
  #
  # Conduct linear model, test lines
  # then make plots.
  
  fit.int <- lm(NegLGI ~ FAvalue*Group, data = df_lm)
  
  capture.output(
    summary(fit.int), 
    file = paste0(statsDir_lm, 
                  "Stats_LM-", h_str, "_", tract, "_G", gType, ".txt"
                  )
  )
  capture.output(
    anova(fit.int), 
    file = paste0(statsDir_lm, 
                  "Stats_AN-", h_str, "_", tract, "_G", gType, ".txt"
                  )
  )
  
  # plot if group diff
  anova_stat <- anova(fit.int)
  if(anova_stat$`Pr(>F)`[2] < 0.05){
  
    comp1 <- substr(comp, start=1, stop=1)
    comp2 <- substr(comp, start=2, stop=2)
    
    if(gType == 1){
      h_cols = c(func_switch_g1(comp1)[[1]][1], func_switch_g1(comp2)[[1]][1])
      names(h_cols) <- c(comp1, comp2)
      h_breaks <- c(comp1, comp2)
      h_labels <- c(func_switch_g1(comp1)[[2]][1], func_switch_g1(comp2)[[2]][1])
      
    }else if(gType == 2){
      h_cols = c(func_switch_g2(comp1)[[1]][1], func_switch_g2(comp2)[[1]][1])
      names(h_cols) <- c(comp1, comp2)
      h_breaks <- c(comp1, comp2)
      h_labels <- c(func_switch_g2(comp1)[[2]][1], func_switch_g2(comp2)[[2]][1])
    }
    
    h_tract <- func_switch_name(tract)
    # h_insert <- paste(h_str, h_tract)
    # h_title <- paste0("Memory Index Predicted by ", h_insert, " Spline Difference")
    
    # p <- ggplot(df_lm) +
    #   aes(x=FAvalue, y=NegLGI, color=Group) +
    #   geom_point(aes(color=Group)) +
    #   geom_smooth(method = "lm") +
    #   ggtitle(h_title) +
    #   ylab("Neg LGI") +
    #   xlab("FA value")
    #   
    # p + scale_color_manual(
    #   values = h_cols,
    #   breaks = h_breaks, 
    #   labels = h_labels
    # )
    # 
    # ggsave(paste0(plotDir_lm, "Plot_LM-", h_str, "_", tract, "_G", gType, ".png"))
    
    # --- update, plot separately
    plot1.x <- df_lm[which(df_lm$Group == comp1),]$FAvalue
    plot1.y <- df_lm[which(df_lm$Group == comp1),]$NegLGI
    
    plot2.x <- df_lm[which(df_lm$Group == comp2),]$FAvalue
    plot2.y <- df_lm[which(df_lm$Group == comp2),]$NegLGI
    
    h_title <- paste(h_tract, "Spline Differences Predicting Memory Performance")
    x_title <- ifelse(h_str == "Avg", "Mean FA", "Max FA")
    
    png(filename = paste0(
        plotDir_lm, "Plot_LM-", h_str, "_", tract, "_G", gType, "_", comp, ".png"
      ),
      width = 8,
      height = 4,
      units = "in",
      res = 600
    )
    
    par(mfrow=c(1,2), oma=c(0,0,2,0), family="Times New Roman")
    
    plot(plot1.x, plot1.y, 
         xlab = x_title, 
         ylab = "Neg LGI", 
         ylim = c(min(df_lm$NegLGI), max(df_lm$NegLGI)),
         main = h_labels[1])
    abline(lm(plot1.y ~ plot1.x))
    
    plot(plot2.x, plot2.y, 
         xlab = x_title, 
         ylab = "",
         main = h_labels[2], 
         ylim = c(min(df_lm$NegLGI), max(df_lm$NegLGI)))
    abline(lm(plot2.y ~ plot2.x))
    
    mtext(h_title, outer = T, cex = 1.5)
    dev.off()
    par(mfrow=c(1,1))
  }
}

func_stat_lm_new <- function(df_lm, tract, gType, h_str, comp){
  
  ### --- Notes:
  #
  # Conduct linear model, test lines
  # then make plots.
  
  # "NeuLGI", "NeuLDI", 
  memList <- c("NegLGI", "NegLDI")
  
  for(mem in memList){
    
    ind_mem <- grep(mem, colnames(df_lm))
    
    df_mem <- as.data.frame(matrix(NA, nrow = dim(df_lm)[1], ncol=4))
    colnames(df_mem) <- c("Comp", "FAvalue", "Group", "MemScore")
    
    df_mem$Comp <- df_lm$Comp
    df_mem$FAvalue <- df_lm$FAvalue
    df_mem$Group <- df_lm$Group
    df_mem$MemScore <- df_lm[,ind_mem]
    
    fit.int <- lm(MemScore ~ FAvalue*Group, data = df_mem)
    
    capture.output(
      summary(fit.int),
      file = paste0(statsDir_lm,
        "Stats_LM-", h_str, "_", tract, "_G", gType, "_", comp, "_", mem, ".txt"
      )
    )
    capture.output(
      anova(fit.int),
      file = paste0(statsDir_lm,
        "Stats_AN-", h_str, "_", tract, "_G", gType, "_", comp, "_", mem, ".txt"
      )
    )
    
    # plot if group diff
    # anova_stat <- anova(fit.int)
    # if(anova_stat$`Pr(>F)`[2] < 0.05){
      
      comp1 <- substr(comp, start=1, stop=1)
      comp2 <- substr(comp, start=2, stop=2)
      
      if(gType == 1){
        h_cols = c(func_switch_g1(comp1)[[1]][1], func_switch_g1(comp2)[[1]][1])
        names(h_cols) <- c(comp1, comp2)
        h_breaks <- c(comp1, comp2)
        h_labels <- c(func_switch_g1(comp1)[[2]][1], func_switch_g1(comp2)[[2]][1])
        
      }else if(gType == 2){
        h_cols = c(func_switch_g2(comp1)[[1]][1], func_switch_g2(comp2)[[1]][1])
        names(h_cols) <- c(comp1, comp2)
        h_breaks <- c(comp1, comp2)
        h_labels <- c(func_switch_g2(comp1)[[2]][1], func_switch_g2(comp2)[[2]][1])
      }
      
      h_tract <- func_switch_name(tract)
      # h_insert <- paste(h_str, h_tract)
      # h_title <- paste0("Memory Index Predicted by ", h_insert, " Spline Difference")
      
      # p <- ggplot(df_lm) +
      #   aes(x=FAvalue, y=NegLGI, color=Group) +
      #   geom_point(aes(color=Group)) +
      #   geom_smooth(method = "lm") +
      #   ggtitle(h_title) +
      #   ylab("Neg LGI") +
      #   xlab("FA value")
      #   
      # p + scale_color_manual(
      #   values = h_cols,
      #   breaks = h_breaks, 
      #   labels = h_labels
      # )
      # 
      # ggsave(paste0(plotDir_lm, "Plot_LM-", h_str, "_", tract, "_G", gType, ".png"))
      
      # --- update, plot separately
      plot1.x <- df_mem[which(df_mem$Group == comp1),]$FAvalue
      plot1.y <- df_mem[which(df_mem$Group == comp1),]$MemScore
      
      plot2.x <- df_mem[which(df_mem$Group == comp2),]$FAvalue
      plot2.y <- df_mem[which(df_mem$Group == comp2),]$MemScore
      
      h_title <- paste(h_tract, "Spline Differences Predicting Memory Performance")
      x_title <- ifelse(h_str == "Avg", "Mean FA", "Max FA")
      
      png(filename = paste0(
          plotDir_lm, "Plot_LM-", h_str, "_", tract, "_G", gType, "_", comp, "_", mem, ".png"
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
           ylim = c(min(df_mem$MemScore), max(df_mem$MemScore)),
           main = h_labels[1])
      abline(lm(plot1.y ~ plot1.x))
      
      plot(plot2.x, plot2.y, 
           xlab = x_title, 
           ylab = "",
           main = h_labels[2], 
           ylim = c(min(df_mem$MemScore), max(df_mem$MemScore)))
      abline(lm(plot2.y ~ plot2.x))
      
      mtext(h_title, outer = T, cex = 1.5)
      dev.off()
      par(mfrow=c(1,1))
    # }
  }
}



func_plot_diff1 <- function(model, tract){
  
  ### --- Notes:
  #
  # This will draw plots and write tables of sig
  #   node differences for GAM when group style=1
  
  png(filename = paste0(
    plotDir_gam, "Plot_Diff_", tract, "_", "G", g_type, ".png"), 
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
                               g_type, "_01.txt"
                 )
  )
  par(mar=c(5,4,4,2))
  dev.off()
}

func_plot_diff2 <- function(model, tract){
  
  ### --- Notes:
  #
  # This will draw plots and write tables of sig
  #   node differences for GAM when group style=2
  
  png(filename = paste0(
      plotDir_gam, "Plot_Diff_", tract, "_", "G", g_type, ".png"
    ), width = 1800, height = 600
  )
  
  gA <- func_switch_g2("0")[[2]][1]
  gB <- func_switch_g2("1")[[2]][1]
  gC <- func_switch_g2("2")[[2]][1]
  
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
                               g_type, "_01.txt"
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
                               g_type, "_02.txt"
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
                               g_type, "_12.txt"
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


func_stat_diff_new <- function(model, tract, g_type){
  
  # 1) func_plot_diff
  #     make plots, write tables
  #     determine diff nodes
  #
  # 2) func_max_diff
  #     make dataframes
  #
  # 3) determine nodes where groups differ
  #     return df_diff.
  #     Account for 2/3 groups
  
  ## Step 1 - make plots and tables
  if(g_type == 1){
    func_plot_diff1(model, tract)
  }else if(g_type == 2){
    func_plot_diff2(model, tract)
  }
  
  # use tables to get sig nodes, save in df_nodeDiff
  df_nodeDiff <- as.data.frame(matrix(NA, nrow=1, ncol=4))
  colnames(df_nodeDiff) <- c("Comparison", "Section", "Start", "End")
  
  if(g_type == 1){
    compList <- "01"
  }else if(g_type == 2){
    compList <- c("01", "02", "12")
  }
  
  # get table made by func_plot_diff for e/comparison
  for(comp in compList){
    
    # write, use bash
    h_cmd = paste0(
      "tail -n +10 ", 
      tableDir, "Table_Diff_", tract, "_", "G", g_type, "_", comp, 
      ".txt | sed 's/-/,/g'"
    )
    h_lines <- system(h_cmd, intern = T)
    
    # make df
    h_df <- read.table(
      text=paste(h_lines, collapse = "\n"), 
      header = F, 
      stringsAsFactors = F, 
      sep = ","
    )
    
    # write to df_nodeDiff
    for(i in 1:dim(h_df)[1]){
      df_nodeDiff <- rbind(df_nodeDiff, c(comp, i, h_df[i,1], h_df[i,2]))
    }
  }
  df_nodeDiff <- na.omit(df_nodeDiff)
  
  
  ## Step 2 - get plot_diff dataframes
  if(g_type == 1){
    df_estDiff <- func_mkdf_diff1(model)
  }else if(g_type == 2){
    df_estDiff <- func_mkdf_diff2(model)
  }
  
  
  ## Step 3 - for grouping 2, determine
  #   intx nodes
  nodeList <- unique(df_allDiff$nodeID)
  
  
  
}



### Work
#   Two analyses (grouping types):
#     1) Con vs Anx
#     2) Con vs GAD vs SAD
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
    df_diff <- func_stat_diff(gam_model, tract, gType)
    compList <- unique(df_diff$Comparison)
    
    # pairwise lm, since plot_diff does pairwise spline tests
    for(comp in compList){
      
      # predict mem score from dti average diff
      df_lm <- func_df_avg(comp, df_tract, df_diff)
      # func_stat_lm(df_lm, tract, gType, "Avg", comp)
      func_stat_lm_new(df_lm, tract, gType, "Avg", comp)
      
      # predict mem score from dti max diff
      #   just rerun plot_diff to get table of max diff
      df_max <- func_max_diff(gam_model, gType)
      df_lm <- func_df_max(comp, df_tract, df_max)
      # func_stat_lm(df_lm, tract, gType, "Max", comp)
      func_stat_lm_new(df_lm, tract, gType, "Max", comp)
    }
  }
}



# # For working out df construction syntax
#
# df_adis <- read.delim(paste0(privateDir, "emuR01_adis.csv"), sep = ",", header=T)
# df_pds <- read.delim(paste0(privateDir, "emuR01_pds_latest.csv"), sep = ",", header=T)
# 
# subjList <- unique(df_tract$subjectID)
# df_practice <- as.data.frame(matrix(NA, nrow=length(subjList), ncol=9))
# colnames(df_practice) <- c("Subject", "Age", "Sex", "PDS", "Pars6", "ADIS.1", "ADIS.2", "ADIS.3", "ADIS.4")
# df_practice$Subject <- subjList
# 
# for(subj in subjList){
# 
#   ind_subj <- grep(subj, df_tract$subjectID)[1]
#   ind_adis <- grep(subj, df_adis$Participant.ID)
#   ind_pds <- grep(subj, df_pds$emu_study_id)
#   ind_out <- grep(subj, df_practice$Subject)
# 
#   df_practice[ind_out,]$Age <- df_tract[ind_subj,]$Age
#   df_practice[ind_out,]$Sex <- df_pds[ind_pds,]$pinf_gender
#   df_practice[ind_out,]$PDS <- as.integer(df_tract[ind_subj,]$PDS)
#   df_practice[ind_out,]$Pars6 <- df_tract[ind_subj,]$Pars6
#   df_practice[ind_out,6:9] <- df_adis[ind_adis,3:6]
# }
# # write.csv(df_practice, file="~/Desktop/dana_table.csv", quote = F, row.names=F, col.names = T)
# 
# g_type <- 2
# df_practice$Group <- NA
# for(i in 1:dim(df_practice)[1]){
#   
#   ind_adis <- which(df_adis$Participant.ID == df_practice[i,]$Subject)
#   
#   if(g_type == 1){
#     h_search <- c("Anxiety", "Phobia")
#     if(sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis,])) != 0){
#       df_practice[i,]$Group <- "Anx"
#     }else if(length(grep("None", df_adis[ind_adis,])) != 0){
#       df_practice[i,]$Group <- "Con"
#     }else{df_practice[i,]$Group <- "Excl"}
#     
#   }else if(g_type == 2){
#     
#     h_search <- c("Separation", "Social")
#     
#     # GAD in dx.1, or in dx.2 but SAD not dx.1
#     if(
#       grepl("Gen", df_adis[ind_adis,]$Diagnosis.1) == T |
#       (sum(grep("Gen", df_adis[ind_adis,])) != 0 &
#        sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis,])) == 0)
#       ){
#       df_practice[i,]$Group <- "GAD"
#     }else if(sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis,])) != 0){
#       df_practice[i,]$Group <- "SAD"
#     }else if(length(grep("None", df_adis[ind_adis,])) != 0){
#       df_practice[i,]$Group <- "Con"
#     }else{df_practice[i,]$Group <- "Excl"}
#       
#   }
# }






