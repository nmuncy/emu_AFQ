

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")
library("ez")
library("dplyr")
library("lme4")
library("ggpubr")
library("lsr")


### Set Up
dataDir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
privateDir <- "/Users/nmuncy/Projects/emu_data/emu_private/"

plotDir_gam <- paste0(dataDir, "plots_gam/")
statsDir_gam <- paste0(dataDir, "stats_gam/")
plotDir_lm <- paste0(dataDir, "plots_lm/")
statsDir_lm <- paste0(dataDir, "stats_lm/")
tableDir <- paste0(dataDir, "tables/")
memoryDir <- paste0(dataDir, "memory/")


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
    #       (not used)
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

func_switch_plot <- function(value, gType){
  
  ### --- Notes:
  #
  # Switch for determining group, coloring

  x_col <- switch(
    value,
    "0" = "blue",
    "1" = "darkred",
    "2" = "black"
  )
  
  if(gType == 2){
    x_label <- switch(
      value,
      "0" = "Con",
      "1" = "GAD",
      "2" = "SAD"
    )
  }else if(gType == 3){
    x_label <- switch(
      value,
      "0" = "Low",
      "1" = "Med",
      "2" = "High"
    )
  }
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

func_stat_mem <- function(df_afq, gType){
  
  ### --- Notes:
  #
  # This function will conduct a MANOVA
  #   and post-hoc on memory measures.
  
  # subset df for ease of use
  ind_mem <- which(df_afq$tractID == tractList[1] & df_afq$nodeID == 0)
  df_mem <- as.data.frame(df_afq[ind_mem,])
  subjList <- df_mem$subjectID
  
  # DF for Group * Neg/Neu Detect/Generalize indices
  #   update, only do LGI
  df_long <-  as.data.frame(matrix(NA,nrow=4*length(subjList), ncol=5))
  colnames(df_long) <- c("Subj", "Group", "Valence", "Mem", "Value")
  
  df_long$Subj <- rep(subjList, 4)
  df_long$Group <- rep(df_mem$Group, 4)
  df_long$Valence <- c(
    rep("Neg", 2*length(subjList)), 
    rep("Neu", 2*length(subjList))
  )
  df_long$Mem <- c(
    rep("LGI", length(subjList)),
    rep("LDI", length(subjList)),
    rep("LGI", length(subjList)),
    rep("LDI", length(subjList))
  )
  df_long$Value <- c(
    df_mem$NegLGI, 
    df_mem$NegLDI,
    df_mem$NeuLGI,
    df_mem$NeuLDI
  )

  # stats
  stats_mem <- ezANOVA(df_long,
   dv = Value,
   wid = Subj,
   between = Group,
   within = c(Valence, Mem),
   type = "III"
  )
  
  # write
  capture.output(
    stats_mem, 
    file = paste0(memoryDir, "Stats_G", gType, "_MANOVA.txt")
  )
  
  # # post-hoc
  # intx_p <- stats_mem$`Sphericity Corrections`$`p[GG]`[2]
  # if(intx_p < 0.05){
  #   h.man <- manova(
  #     cbind(NeuLGI, NegLGI, NeuLDI, NegLDI) ~ Group, 
  #     data = df_mem
  #   )
  #   capture.output(
  #     summary.aov(h.man), 
  #     file = paste0(memoryDir, "Stats_G", gType, "_post.txt")
  #   )
  # }
  # 
  # # Tukey for g3 NegLGI
  # if(gType == 3){
  #   
  #   group_fac <- as.numeric(df_mem$Group)
  #   group_value <- vector()
  #   for(val in group_fac){
  #     group_value <- c(group_value, func_switch_plot(val, gType)[[2]])
  #   }
  #   df_mem$Group <- group_value
  # 
  #   h.lm <- lm(NegLGI ~ Group, data = df_mem)
  #   h.av <- aov(h.lm)
  #   capture.output(
  #     TukeyHSD(h.av), 
  #     file = paste0(memoryDir, "Stats_G", gType, "_tuk.txt")
  #   )
  #   
  #   # make a plot
  #   h_colors <- c(func_switch_plot("2", gType)[[1]], 
  #     func_switch_plot("0", gType)[[1]], 
  #     func_switch_plot("1", gType)[[1]]
  #   )
  #   post_comp <- list(c("High", "Low"))
  #   ggboxplot(df_mem, 
  #     y="NegLGI", 
  #     x="Group", 
  #     color="Group", 
  #     palette = h_colors, 
  #     add = "jitter",
  #     title = "top"
  #   ) +
  #     stat_compare_means(comparisons = post_comp) +
  #     ggtitle("Memory Performance on Negative Stimuli") +
  #     rremove("legend") +
  #     theme(text=element_text(family="Times New Roman", face="bold", size=16))
  #   
  #   ggsave(
  #     paste0(memoryDir, "Plot_Box_G", gType, "_NegLGI.tiff"),
  #     units = "in",
  #     width = 6,
  #     height = 6,
  #     device = "tiff"
  #   )
  # }
}


# Functions - GAM
func_plot_gam <- function(model, tract, gType, df_tract){
  
  ### --- Notes:
  #
  # Will plot the GAM model of dti data
  
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
  
  h_cols = c(
    func_switch_plot("0", gType)[[1]][1], 
    func_switch_plot("1", gType)[[1]][1], 
    func_switch_plot("2", gType)[[1]][1]
  )
  names(h_cols) <- c("0", "1", "2")
  h_breaks <- c("0", "1", "2")
  h_labels <- c(
    func_switch_plot("0", gType)[[2]][1], 
    func_switch_plot("1", gType)[[2]][1], 
    func_switch_plot("2", gType)[[2]][1]
  )

  p <- ggplot(data = df_pred) +
    geom_smooth(mapping = aes(x=nodeID, y=fit, color=Group)) +
    ggtitle(h_title) +
    ylab("Fit FA") +
    xlab("Tract Node") +
    theme(text=element_text(family="Times New Roman", face="bold", size=14))
  
  p + scale_color_manual(
    values = h_cols,
    breaks = h_breaks, 
    labels = h_labels
  )
  
  ggsave(
    paste0(plotDir_gam, "Plot_GAM_", tract, "_", "G", gType, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

func_stat_gam <- function(tract, df_tract, gType){
  
  ### --- Notes:
  #
  # This function has 2 main steps
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
  # Note: this work can only be a function since
  #   it happened that all tracts had the same
  #   distribution, family, K, and best model.

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
  capture.output(
    summary(fit_gamma), 
    file = paste0(
      statsDir_gam, "Stats_GAM-gamma_", tract, "_", "G", gType, ".txt"
    )
  )

  fit_beta <- bam(dti_fa ~ Group +
                    Sex +
                    s(nodeID, by=Group, k=40) +
                    s(subjectID, bs="re"),
                  data = df_tract,
                  family = betar(link = "logit"),
                  method = "REML")

  # gam.check(fit_beta, rep = 500)
  capture.output(
    summary(fit_beta), 
    file = paste0(
      statsDir_gam, "Stats_GAM-beta_", tract, "_", "G", gType, ".txt"
    )
  )
  
  # infoMessages('on')
  # compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  capture.output(
    compareML(fit_gamma, fit_beta), 
    file = paste0(
      statsDir_gam, "Stats_GAM-comp_gam-beta_", tract, "_", "G", gType, ".txt"
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
  capture.output(
    summary(fit_cov_pds), 
    file = paste0(
      statsDir_gam, "Stats_GAM-cov_", tract, "_", "G", gType, ".txt"
    )
  )
  
  # Test cov model against gamma
  capture.output(
    compareML(fit_gamma, fit_cov_pds), 
    file = paste0(
      statsDir_gam, "Stats_GAM-comp_gam-cov_", tract, "_", "G", gType, ".txt"
    )
  )
  return(fit_cov_pds)
}

func_plot_diff_pair <- function(model, tract, gType, factorA, factorB){
  
  ### --- Notes:
  #
  # This will make plots and write tables of sig
  #   node differences for GAM splines bx 2 factors (groups)
  
  png(filename = paste0(
    plotDir_gam, "Plot_Diff_", tract, "_", "G", gType, "_pair.png"), 
    width = 600, height = 600
  )
  
  gA <- func_switch_plot(factorA, gType)[[2]][1]
  gB <- func_switch_plot(factorB, gType)[[2]][1]
  
  par(mar=c(5,5,4,2), family="Times New Roman")
  capture.output(plot_diff(model,
                           view="nodeID",
                           comp=list(Group=c(factorA, factorB)),
                           rm.ranef = T,
                           main = paste0("Difference Scores, ", gA, "-", gB),
                           ylab = "Est. FA difference",
                           xlab = "Tract Node",
                           cex.lab = 2,
                           cex.axis = 2,
                           cex.main = 2,
                           cex.sub = 1.5),
                 file = paste0(tableDir, 
                               "Table_Diff_", 
                               tract, "_", "G",
                               gType, "_",
                               factorA, factorB,
                               ".txt"
                 )
  )
  par(mar=c(5,4,4,2))
  dev.off()
}

func_plot_diff_anova <- function(model, tract, gType){
  
  ### --- Notes:
  #
  # This will draw plots and write tables of sig
  #   node differences for GAM when group style=2
  
  tiff(filename = paste0(
      plotDir_gam, "Plot_Diff_", 
      tract, "_", "G", gType, "_anova.tiff"
    ), 
    width = 2000, 
    height = 600
  )
  
  gA <- func_switch_plot("0", gType)[[2]][1]
  gB <- func_switch_plot("1", gType)[[2]][1]
  gC <- func_switch_plot("2", gType)[[2]][1]
  
  par(mfrow=c(1,3), mar=c(5,6,4,2), family="Times New Roman")
  capture.output(plot_diff(model,
                           view="nodeID",
                           comp=list(Group=c("0", "1")),
                           rm.ranef = T,
                           main = paste0("Difference Scores, ", gA, "-", gB),
                           ylab = "Est. FA difference",
                           xlab = "Tract Node",
                           cex.lab = 3,
                           cex.axis = 3,
                           cex.main = 3.5,
                           cex.sub = 2.5),
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
                           cex.lab = 3,
                           cex.axis = 3,
                           cex.main = 3.5,
                           cex.sub = 2.5),
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
                           cex.lab = 3,
                           cex.axis = 3,
                           cex.main = 3.5,
                           cex.sub = 2.5),
                 file = paste0(tableDir, 
                               "Table_Diff_", 
                               tract, "_", "G", 
                               gType, "_12.txt"
                 )
  )
  par(mfrow=c(1,1), mar=c(5,4,4,2))
  dev.off()
}

func_mkdf_diff_pair <- function(model, factorA, factorB){
  
  ### --- Notes:
  #
  # Returns a dataframe of difference
  #   scores for each node.
  
  df_pair <- plot_diff(model,
                   view="nodeID",
                   comp=list(Group=c(factorA, factorB)),
                   rm.ranef = T,
                   plot = F)
  
  colnames(df_pair) <- c(colnames(df_pair[,1:4]), "Comp")
  df_pair$Comp <- paste0(factorA, factorB)
  return(df_pair)
}

func_mkdf_diff_anova <- function(model){
  
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

func_stat_diff <- function(model, tract, gType, pairAn, comp_list){
  
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
  if(pairAn == "pair"){
    func_plot_diff_pair(model, tract, gType, comp_list[1], comp_list[2])
  }else if(pairAn == "anova"){
    func_plot_diff_anova(model, tract, gType)
  }
  
  # get plot_diff data frames
  if(pairAn == "pair"){
    df_estDiff <- func_mkdf_diff_pair(model, comp_list[1], comp_list[2])
  }else if(pairAn == "anova"){
    df_estDiff <- func_mkdf_diff_anova(model)
  }
  
  # determine where nodes differ, anova
  #   looks for nodes that differ bx ALL groups
  nodeList <- unique(df_estDiff$nodeID)
  diffList <- vector()
  for(node in nodeList){
    ind_node <- which(df_estDiff$nodeID == node)
    if(pairAn == "pair"){
      h_est <- abs(df_estDiff[ind_node[1],]$est)
      h_ci <- df_estDiff[ind_node[1],]$CI
      if((h_est - h_ci) > 0){
        diffList <- c(diffList, node)
      }
    }else if(pairAn == "anova"){
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
  
  # find node of max difference
  if(pairAn == "pair"){
    
    h_df <- subset(df_estDiff, nodeID %in% diffList)
    ind_max <- which(abs(h_df$est) == max(abs(h_df$est)))
    node_max <- h_df[ind_max,]$nodeID
    
  }else if(pairAn == "anova"){
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

func_plot_lm_pair <- function(df_plot, avg_max, mem, factorA, factorB){
  
  ### --- Notes:
  #
  # Plots LM of A and B
  
  # determine labels
  h_labels <- c(
    func_switch_plot(factorA, gType)[[2]][1], 
    func_switch_plot(factorB, gType)[[2]][1]
  )
  h_tract <- func_switch_name(tract)
  
  # set up for plot A
  plot1.x <- df_plot[which(df_plot$Group == as.numeric(factorA)),]$FAvalue
  plot1.y <- df_plot[which(df_plot$Group == as.numeric(factorA)),]$MemScore
  
  # set up for plot B
  plot2.x <- df_plot[which(df_plot$Group == as.numeric(factorB)),]$FAvalue
  plot2.y <- df_plot[which(df_plot$Group == as.numeric(factorB)),]$MemScore
  
  # omni title
  h_title <- paste(h_tract, "Spline Differences Predicting Memory Performance")
  x_title <- ifelse(avg_max == "Avg", "Mean FA", "Max FA")
  
  # set up, draw
  tiff(filename = paste0(
      plotDir_lm, 
      "Plot_LM-", 
      avg_max, "_", 
      tract, "_G",
      gType, "_", 
      mem, "_pair.tiff"
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

func_plot_lm_anova <- function(df_plot, avg_max, mem, gType){
  
  ### --- Notes:
  #
  # Plots A, B, and C
  
  h_labels <- c(
    func_switch_plot("0", gType)[[2]][1], 
    func_switch_plot("1", gType)[[2]][1], 
    func_switch_plot("2", gType)[[2]][1]
  )
  h_tract <- func_switch_name(tract)
  
  plot1.x <- df_plot[which(df_plot$Group == 0),]$FAvalue
  plot1.y <- df_plot[which(df_plot$Group == 0),]$MemScore
  
  plot2.x <- df_plot[which(df_plot$Group == 1),]$FAvalue
  plot2.y <- df_plot[which(df_plot$Group == 1),]$MemScore
  
  plot3.x <- df_plot[which(df_plot$Group == 2),]$FAvalue
  plot3.y <- df_plot[which(df_plot$Group == 2),]$MemScore
  
  h_title <- paste(h_tract, "Spline Differences Predicting Memory Performance")
  x_title <- ifelse(avg_max == "Avg", "Mean FA", "Max FA")
  
  tiff(filename = paste0(
      plotDir_lm, "Plot_LM-", 
      avg_max, "_", 
      tract, "_G", 
      gType, "_", 
      mem, "_anova.tiff"
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
       cex.lab = 1.5,
       ylab = mem, 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)),
       main = h_labels[1])
  abline(lm(plot1.y ~ plot1.x))
  
  plot(plot2.x, plot2.y, 
       xlab = x_title, 
       cex.lab = 1.5,
       ylab = "",
       main = h_labels[2], 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)))
  abline(lm(plot2.y ~ plot2.x))
  
  plot(plot3.x, plot3.y, 
       xlab = x_title, 
       cex.lab = 1.5,
       ylab = "",
       main = h_labels[3], 
       ylim = c(min(df_plot$MemScore), max(df_plot$MemScore)))
  abline(lm(plot3.y ~ plot3.x))
  
  mtext(h_title, outer = T, cex = 1.5)
  dev.off()
  par(mfrow=c(1,1))
}

func_stat_lm <- function(df_lm, tract, gType, avg_max, pairAn, comp_list){
  
  ### --- Notes:
  #
  # Conduct linear model for list of mem scores
  #   then make plots.
  
  # memList <- c("NeuLGI", "NeuLDI","NegLGI", "NegLDI")
  memList <- "NegLGI"
  
  for(mem in memList){
    
    # make dataframe for memory behavior
    df_mem <- as.data.frame(matrix(NA, nrow = dim(df_lm)[1], ncol=4))
    colnames(df_mem) <- c("Subj", "FAvalue", "Group", "MemScore")
    
    df_mem$Subj <- df_lm$Subj
    df_mem$FAvalue <- df_lm$FAvalue
    df_mem$Group <- df_lm$Group
    
    # get/add appropriate behavior score
    ind_mem <- grep(mem, colnames(df_lm))
    df_mem$MemScore <- df_lm[,ind_mem]
    
    # subset for pairwise comparison
    if(pairAn == "pair"){
      df_mem <- df_mem[which(
        df_mem$Group == comp_list[1] | df_mem$Group == comp_list[2]
      ),]
    }
    
    # linear regression
    fit.int <- lm(MemScore ~ FAvalue*Group, data = df_mem)
    capture.output(
      summary(fit.int),
      file = paste0(statsDir_lm,
                    "Stats_LM-", avg_max, "_", 
                    tract, "_G", gType, "_", 
                    mem, "_", pairAn, ".txt"
      )
    )
    etaSquared(fit.int)
    
    # make plots
    if(pairAn == "pair"){
      func_plot_lm_pair(df_mem, avg_max, mem, comp_list[1], comp_list[2])
    }else if(pairAn == "anova"){
      func_plot_lm_anova(df_mem, avg_max, mem, gType)
    }
  }
}



### --- Work:
#
# Two analyses (grouping types):
#   1) Con vs Anx (removed)
#   2) Con vs GAD vs SAD
#   3) PARS Low vs Med vs High
groupType <- 3
tractList <- c("UNC_L", "UNC_R", "FA")

# for spline comparisons, linear models,
#   determine whether to compare all (anova)
#   or a pair. If pair, set factors as string.
pairAn <- "pair"
g2_pairList <- c("0", "2")
g3_pairList <- c("0", "2")


for(gType in groupType){

  # make/get data, assign factor
  dataFile <- paste0(dataDir, "Master_dataframe_G", gType,".csv")
  
  if( ! file.exists(dataFile)){
    func_df_master(gType)
  }
  
  df_afq <- read.csv(dataFile)
  df_afq$Group <- factor(df_afq$Group)
  df_afq$Sex <- factor(df_afq$Sex)
  
  # Check Memory behavior
  func_stat_mem(df_afq, gType)
  
  # get comp list
  if(pairAn == "pair"){
    h_var <- paste0("g", gType, "_pairList")
    comp_list <- get(h_var)
  }
  
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
    
    # determine nodes of group differences
    #   Use "anova" to compare bx >2 groups, and note
    #     that differences between splines will be for
    #     nodes which differed between ALL groups
    nodeList <- func_stat_diff(gam_model, tract, gType, pairAn, comp_list)
    
    # deal w/no differences
    if(is.na(nodeList)){
      next
    }
    
    # avg lm
    avg_nList <- nodeList[[1]]
    df_avg <- func_mkdf_lm(df_tract, avg_nList, gType, "Avg")
    func_stat_lm(df_avg, tract, gType, "Avg", pairAn, comp_list)
    
    # # max lm
    # max_nList <- nodeList[[2]]
    # df_max <- func_mkdf_lm(df_tract, max_nList, gType, "Max")
    # func_stat_lm(df_max, tract, gType, "Max", pairAn, comp_list)
  }
}
