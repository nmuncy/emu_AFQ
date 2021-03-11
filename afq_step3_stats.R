

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")
library("ez")
library("dplyr")
library("lme4")


# Functions
func_makeDF <- function(){
  
  ### --- Notes
  #
  # This function will make a master dataframe that
  # contains node FA values, group membership, sex, PDS
  # and memory measures (LG/DI)
  #
  # Writes dataDir/Master_dataframe.csv
  
  # Get data
  df_afq <- read.delim(paste0(dataDir, "tract_profiles.csv"), sep = ",", header = T)
  colnames(df_afq) <- c("Counter", colnames(df_afq)[-1])
  df_full <- read.delim(paste0(privateDir, "emuR01_full_latest.csv"), sep = ",", header=T)
  df_pds <- read.delim(paste0(privateDir, "emuR01_pds_latest.csv"), sep = ",", header=T)
  
  # make lists
  subjList <- unique(df_afq$subjectID)
  tractList <- unique(df_afq$tractID)
  nodeList <- unique(df_afq$nodeID)
  
  # add group, pars6, pds, age
  #   add d-primes, sex
  df_afq$Group <- df_afq$Pars6 <- df_afq$PDS <- df_afq$Age <- NA
  # df_afq$NegDP1wk <- df_afq$NeuDP1wk <- NA
  df_afq$NegLDI <- df_afq$NeuLDI <- df_afq$NegLGI <- df_afq$NeuLGI <- NA
  df_afq$Sex <- NA
  
  for(subj in subjList){
    
    # determine subj indices
    ind_afq <- which(df_afq$subjectID == subj)
    ind_full <- which(df_full$src_subject_id == subj)
    ind_pds <- which(df_pds$emu_study_id == subj)
    
    # determine group
    #   low=0-3, med=4-12, high>12
    h_anx <- df_full[ind_full,]$pars_6
    
    if(h_anx <= 3){
      h_group <- 0
    }else if(h_anx > 3 & h_anx < 13){
      h_group <- 1
    }else if(h_anx > 12){
      h_group <- 2
    }
    
    # determine pds
    h_pds <- df_pds[ind_pds,]$pds_shirtcliff
    
    # determine age, sex
    h_age <- df_full[ind_full,]$pinf_age
    h_sex <- substr(df_full[ind_full,]$sex, 1, 1)
    if(h_sex == "f"){h_sexF <- 0}else if(h_sex == "m"){h_sexF <- 1}
    
    # Get Beh counts
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
        assign(check, 0.1)
      }
    }
    
    # Calculate Neg, Neu d'
    # dp_neg <- qnorm(neg_num_hit/(neg_num_hit + neg_num_miss)) - 
    #   qnorm(neg_num_Lfa/(neg_num_Lfa + neg_num_Lcr))
    # 
    # dp_neu <- qnorm(neu_num_hit/(neu_num_hit + neu_num_miss)) - 
    #   qnorm(neu_num_Lfa/(neu_num_Lfa + neu_num_Lcr))
    
    # Calculate Neg, Neu LDI
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
    
    # df_afq[ind_afq,]$NegDP1wk <- round(dp_neg, 2)
    # df_afq[ind_afq,]$NeuDP1wk <- round(dp_neu, 2)
    df_afq[ind_afq,]$NegLDI <- neg_LDI
    df_afq[ind_afq,]$NeuLDI <- neu_LDI
    df_afq[ind_afq,]$NegLGI <- neg_LGI
    df_afq[ind_afq,]$NeuLGI <- neu_LGI
  }
  
  # write csv
  outFile <- paste0(dataDir, "Master_dataframe.csv")
  write.csv(df_afq, file=outFile, quote=F, row.names = F)
  return(df_afq)
}

func_memStats <- function(){
  
  ### --- Notes
  #
  # This function will check for conduct a
  # WBRM ANOVA on memory measures.
  
  # get data, make lists
  df_afq <- read.csv(paste0(dataDir, "Master_dataframe.csv"))
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
  capture.output(stats_LDI, file = paste0(dataDir, "Stats_AN_LDI.txt"))
  
  
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
  
  capture.output(stats_LGI, file = paste0(dataDir, "Stats_AN_LGI.txt"))
}

func_ggplot_gam <- function(h_df, h_title, outDir, tract){
  
  ggplot(data = h_df) +
    geom_smooth(mapping = aes(x=nodeID, y=fit, color=Group)) +
    ggtitle(h_title) +
    ylab("Fit FA")
  
  ggsave(paste0(outDir, "Plot_GAM_", tract, ".png"))
}

func_plot_diff <- function(h_df, outDir, tract){
  
  ### --- Notes:
  #
  # This function will determine, plot differences
  # between two splines.
  # The difference values will be returned.
  
  png(filename = paste0(outDir, "Plot_Diff_", tract, ".png"), width = 1800, height = 600)
  par(mfrow=c(1,3))
  par(mar=c(5,5,4,2))
  
  p01 <- plot_diff(h_df,
                   view="nodeID",
                   comp=list(Group=c("0", "1")),
                   rm.ranef = T,
                   main = "Difference Scores, Low-Med",
                   ylab = "Est. FA difference",
                   xlab = "Tract Node",
                   cex.lab = 2,
                   cex.axis = 2,
                   cex.main = 2.5,
                   cex.sub = 1.5)
  
  par(mar=c(5,3,4,2))
  
  p02 <- plot_diff(h_df,
                   view="nodeID",
                   comp=list(Group=c("0", "2")),
                   rm.ranef = T,
                   main = "Difference Scores, Low-High",
                   ylab = "",
                   xlab = "Tract Node",
                   cex.lab = 2,
                   cex.axis = 2,
                   cex.main = 2.5,
                   cex.sub = 2)
  
  p12 <- plot_diff(h_df,
                   view="nodeID",
                   comp=list(Group=c("1", "2")),
                   rm.ranef = T,
                   main = "Difference Scores, Med-High",
                   ylab = "",
                   xlab = "Tract Node",
                   cex.lab = 2,
                   cex.axis = 2,
                   cex.main = 2.5,
                   cex.sub = 2)
  
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2))
  dev.off()
  
  return(list(p01,p02,p12))
}

func_gam <- function(tract, df, outDir){
  
  ### This function has 3 main steps
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
  
  # Subset data by tract
  df_tract <- df[which(df$tractID == tract), ]
  df_tract$dti_fa <- round(df_tract$dti_fa, 3)
  
  # plot mean data
  ggplot(data = df_tract) +
    geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

  ggplot(data = df_tract) +
    geom_point(mapping = aes(x=nodeID, y=dti_fa, color=Group),size=0.3) +
    geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

  # determine distribution
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
  capture.output(summary(fit_gamma), file = paste0(outDir, "Stats_GAM-gamma_", tract, ".txt"))
  
  
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
  capture.output(compareML(fit_gamma, fit_cov_pds), file = paste0(outDir, "Stats_GAM-comp_", tract, ".txt"))
  capture.output(summary(fit_cov_pds), file = paste0(outDir, "Stats_GAM-cov_", tract, ".txt"))
  
  # plot
  df_pred <- predict.bam(
    fit_cov_pds,
    exclude_terms = c("PDS", "Sex","subjectID"),
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
  
  if(tract == "UNC_L"){
    h_tract = "L. Uncinate"
  }else if(tract == "UNC_R"){
    h_tract = "R. Uncinate"
  }else if(tract == "CGC_L"){
    h_tract = "L. Cingulum"
  }else if(tract == "CGC_R"){
    h_tract = "R. Cingulum"
  }else if(tract == "ATR_L"){
    h_tract = "L. A. Thalamic Radiations"
  }else if(tract == "ARC_L"){
    h_tract = "L. Arcuate"
  }else if(tract == "ARC_R"){
    h_tract = "R. Arcuate"
  }
  plot_title = paste0("GAM Fit of ", h_tract," FA Values")
  func_ggplot_gam(df_pred, plot_title, outDir, tract)
  
  
  ### Test for differences
  plot_diff <- func_plot_diff(fit_cov_pds, outDir, tract)
  p01 <- plot_diff[[1]]
  p02 <- plot_diff[[2]]
  p12 <- plot_diff[[3]]
  
  # find biggest difference, node location
  df_est <- as.data.frame(matrix(NA, nrow=3*dim(p01)[1], ncol=dim(p01)[2]))
  colnames(df_est) <- colnames(p01)
  df_est[,1:5] <- rbind(p01, p02, p12)
  
  ind_max <- which(abs(df_est$est) == max(abs(df_est$est)))
  node_max <- df_est[ind_max,]$nodeID
  h_groups <- df_est[ind_max,]$comp
  groups <- stringr::str_extract_all(h_groups, "\\d+")
  gA <- as.numeric(groups[[1]][1])
  gB <- as.numeric(groups[[1]][2])
  
  # make df
  df_max <- as.data.frame(df_tract[which(
    df_tract$nodeID == node_max &
      (df_tract$Group == gA | df_tract$Group == gB)
  ),])
  
  return(df_max)
}


# Orienting vars
dataDir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
privateDir <- "/Users/nmuncy/Projects/emu_private/"

# Make dataset
func_makeDF()

# Check Memory behavior
func_memStats()


# Get data for GAMs
df_afq <- read.csv(paste0(dataDir, "Master_dataframe.csv"))
df_afq$Group <- factor(df_afq$Group)
df_afq$Sex <- factor(df_afq$Sex)

# L Unc
tract <- "UNC_L"
df_max <- func_gam(tract, df_afq, dataDir)

# linear models, plot sig
fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
summary(fit)
capture.output(summary(fit), file = paste0(dataDir, "Stats_LM_", tract, ".txt"))

ggplot(df_max, aes(x=dti_fa, y=NegLGI)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Group) +
  ggtitle("FA Values Predicting Memory Outcome")
ggsave(paste0(dataDir, "Plot_LM_", tract, ".png"))

fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
summary(fit)


# R Unc
df_max <- func_gam("UNC_R", df_afq, dataDir)

fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
summary(fit)
fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
summary(fit)


# L Arc
df_max <- func_gam("ARC_L", df_afq, dataDir)

fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
summary(fit)
fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
summary(fit)


# R Arc
df_max <- func_gam("ARC_R", df_afq, dataDir)

fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
summary(fit)
fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
summary(fit)


# L Cing
df_max <- func_gam("CGC_L", df_afq, dataDir)

fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
summary(fit)
fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
summary(fit)


# R Cing
df_max <- func_gam("CGC_R", df_afq, dataDir)

fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
summary(fit)
fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
summary(fit)


# L ATR
df_max <- func_gam("ATR_L", df_afq, dataDir)

fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
summary(fit)
fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
summary(fit)

