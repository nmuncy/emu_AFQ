

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")
library("ez")
library("dplyr")
library("lme4")


# Orienting vars
dataDir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
privateDir <- "/Users/nmuncy/Projects/emu_private/"

do_lunc = 0
do_runc = 0
do_lcing = 1
do_rcing = 1
do_latr = 1


# Plot functions
func_ggplot_gam <- function(h_df, h_title){
  
  ggplot(data = h_df) +
    geom_smooth(mapping = aes(x=nodeID, y=fit, color=Group)) +
    ggtitle(h_title) +
    ylab("Fit FA")
  
  ggsave(paste0(dataDir, "Plot_", tract, "_GAM.png"))
}

func_plot_diff <- function(h_df){
  
  png(filename = paste0(dataDir, "Plot_", tract, "_diff.png"), width = 1800, height = 600)
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


### --- Step 1: Make dataset
func_makeDF <- function(){

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
df_afq <- func_makeDF()


### --- Step 2: Check Memory behavior
func_memStats <- function(){
  
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
func_memStats()


# Get data
df_afq <- read.csv(paste0(dataDir, "Master_dataframe.csv"))
df_afq$Group <- factor(df_afq$Group)
df_afq$Sex <- factor(df_afq$Sex)


# L. Unc
if(do_lunc == 1){
  
  ### --- Step 3: Model tract, no covariates

  # Tract
  tract <- "UNC_L"
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$dti_fa <- round(df_tract$dti_fa, 3)
  
  # plot mean data
  ggplot(data = df_tract) +
    geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

  ggplot(data = df_tract) +
    geom_point(mapping = aes(x=nodeID, y=dti_fa, color=Group),size=0.3) +
    geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

  # determine distribution
  descdist(df_tract$dti_fa, discrete=F) # Could be beta or gamma

  fit.beta <- fitdist(df_tract$dti_fa, "beta")
  plot(fit.beta)
  fit.beta$aic

  fit.gamma <- fitdist(df_tract$dti_fa, "gamma")
  plot(fit.gamma)
  fit.gamma$aic
  
  #  determine k, compare families
  fit_gamma <- bam(dti_fa ~ Group +
                     Sex +
                     s(nodeID, by=Group, k=40) +
                     s(subjectID, bs="re"),
                   data = df_tract,
                   family = Gamma(link = "logit"),
                   method = "REML")
  
  gam.check(fit_gamma, rep = 500)
  
  fit_beta <- bam(dti_fa ~ Group +
                    Sex +
                    s(nodeID, by=Group, k=40) +
                    s(subjectID, bs="re"),
                  data = df_tract,
                  family = betar(link = "logit"),
                  method = "REML")
  
  gam.check(fit_beta, rep = 500)
  
  infoMessages('on')
  compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  
  # get stats
  summary(fit_gamma)  # R-sq = 0.819
  capture.output(summary(fit_gamma), file = paste0(dataDir, tract, "_GAM.txt"))
  
  
  ### --- Step 4: Model tract, covariates
  
  fit_cov_pds <- bam(dti_fa ~ Group +
                       Sex +
                       s(nodeID, by=Group, k=40) +
                       s(PDS, by=Sex) +
                       s(subjectID, bs="re"),
                     data = df_tract,
                     family = Gamma(link = "logit"),
                     method = "REML")
  
  gam.check(fit_cov_pds, rep = 500)
  
  # Test cov model against gamma
  # infoMessages('on')
  capture.output(compareML(fit_gamma, fit_cov_pds), file = paste0(dataDir, tract, "_GAM_comp.txt"))
  capture.output(summary(fit_cov_pds), file = paste0(dataDir, tract, "_GAM_cov.txt"))
  
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
  
  func_ggplot_gam(df_pred, "GAM Fit of L. Uncinate FA Values")
  
  
  ### --- Step 5: Test for differences
  #
  # Check for group differences in spline.
  plot_diff <- func_plot_diff(fit_cov_pds)
  p01 <- plot_diff[[1]]
  p02 <- plot_diff[[2]]
  p12 <- plot_diff[[3]]
  

  ### --- Step 6: Regress

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
  
  # linear models, plot sig
  fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
  # summary(fit)
  capture.output(summary(fit), file = paste0(dataDir, tract, "_lm.txt"))
  colnames(df_max$dti_fa) <- "FA Value"
  ggplot(df_max, aes(x=dti_fa, y=NegLGI)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~ Group) +
    ggtitle("FA Values Predicting Memory Outcome")
  
  ggsave(paste0(dataDir, "Plot_", tract, "_lm.png"))
  
  fit <- lmList(NegLDI ~ dti_fa | Group, data = df_max)
  summary(fit)
}


# R. UNC
if(do_runc == 1){
  
  # Tract
  tract <- "UNC_R"
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$dti_fa <- round(df_tract$dti_fa, 3)
  
  #  determine k, compare families
  fit_gamma <- bam(dti_fa ~ Group +
                     Sex +
                     s(nodeID, by=Group, k=40) +
                     s(subjectID, bs="re"),
                   data = df_tract,
                   family = Gamma(link = "logit"),
                   method = "REML")
  
  gam.check(fit_gamma, rep = 500)
  
  fit_beta <- bam(dti_fa ~ Group +
                    Sex +
                    s(nodeID, by=Group, k=40) +
                    s(subjectID, bs="re"),
                  data = df_tract,
                  family = betar(link = "logit"),
                  method = "REML")
  
  gam.check(fit_beta, rep = 500)
  
  infoMessages('on')
  compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  
  # get stats
  summary(fit_gamma)  # R-sq = 0.89
  capture.output(summary(fit_gamma), file = paste0(dataDir, tract, "_GAM.txt"))
  
  # covariates
  fit_cov_pds <- bam(dti_fa ~ Group +
                       Sex +
                       s(nodeID, by=Group, k=40) +
                       s(PDS, by=Sex) +
                       s(subjectID, bs="re"),
                     data = df_tract,
                     family = Gamma(link = "logit"),
                     method = "REML")
  
  gam.check(fit_cov_pds, rep = 500)
  
  # Test cov model against gamma
  # infoMessages('on')
  capture.output(compareML(fit_gamma, fit_cov_pds), file = paste0(dataDir, tract, "_GAM_comp.txt"))
  capture.output(summary(fit_cov_pds), file = paste0(dataDir, tract, "_GAM_cov.txt"))
  
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
  
  func_ggplot_gam(df_pred, "GAM Fit of R. Uncinate FA Values")
  
  # Check for group differences in spline.
  plot_diff <- func_plot_diff(fit_cov_pds)
  p01 <- plot_diff[[1]]
  p02 <- plot_diff[[2]]
  p12 <- plot_diff[[3]]
  
  # find biggest difference
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
  
  # beh x group
  fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
  fit <- lmList(NeuLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
}


# L. Cing
if(do_lcing == 1){
  
  # Tract
  tract <- "CGC_L"
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$dti_fa <- round(df_tract$dti_fa, 3)

  #  determine k, compare families
  fit_gamma <- bam(dti_fa ~ Group +
                     Sex +
                     s(nodeID, by=Group, k=40) +
                     s(subjectID, bs="re"),
                   data = df_tract,
                   family = Gamma(link = "logit"),
                   method = "REML")
  
  gam.check(fit_gamma, rep = 500)
  
  fit_beta <- bam(dti_fa ~ Group +
                    Sex +
                    s(nodeID, by=Group, k=40) +
                    s(subjectID, bs="re"),
                  data = df_tract,
                  family = betar(link = "logit"),
                  method = "REML")
  
  gam.check(fit_beta, rep = 500)
  
  infoMessages('on')
  compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  
  # get stats
  summary(fit_gamma)  # R-sq = 0.67
  capture.output(summary(fit_gamma), file = paste0(dataDir, tract, "_GAM.txt"))
  
  # covariates
  fit_cov_pds <- bam(dti_fa ~ Group +
                       Sex +
                       s(nodeID, by=Group, k=40) +
                       s(PDS, by=Sex) +
                       s(subjectID, bs="re"),
                     data = df_tract,
                     family = Gamma(link = "logit"),
                     method = "REML")
  
  gam.check(fit_cov_pds, rep = 500)
  
  # Test cov model against gamma
  # infoMessages('on')
  capture.output(compareML(fit_gamma, fit_cov_pds), file = paste0(dataDir, tract, "_GAM_comp.txt"))
  capture.output(summary(fit_cov_pds), file = paste0(dataDir, tract, "_GAM_cov.txt"))
  
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
  
  func_ggplot_gam(df_pred, "GAM Fit of L. Cingulate FA Values")
  
  # Check for group differences in spline.
  plot_diff <- func_plot_diff(fit_cov_pds)
  p01 <- plot_diff[[1]]
  p02 <- plot_diff[[2]]
  p12 <- plot_diff[[3]]
  
  # find biggest difference
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
  
  # negLGI x group
  fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
  # neuLGI x group
  fit <- lmList(NeuLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
}


# R. Cing
if(do_lcing == 1){

  # Tract
  tract <- "CGC_R"
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$dti_fa <- round(df_tract$dti_fa, 3)
  
  #  determine k, compare families
  fit_gamma <- bam(dti_fa ~ Group +
                     Sex +
                     s(nodeID, by=Group, k=40) +
                     s(subjectID, bs="re"),
                   data = df_tract,
                   family = Gamma(link = "logit"),
                   method = "REML")
  
  gam.check(fit_gamma, rep = 500)
  
  fit_beta <- bam(dti_fa ~ Group +
                    Sex +
                    s(nodeID, by=Group, k=40) +
                    s(subjectID, bs="re"),
                  data = df_tract,
                  family = betar(link = "logit"),
                  method = "REML")
  
  gam.check(fit_beta, rep = 500)
  
  infoMessages('on')
  compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  
  # get stats
  summary(fit_gamma)  # R-sq = 0.49
  capture.output(summary(fit_gamma), file = paste0(dataDir, tract, "_GAM.txt"))
  
  # covariates
  fit_cov_pds <- bam(dti_fa ~ Group +
                       Sex +
                       s(nodeID, by=Group, k=40) +
                       s(PDS, by=Sex) +
                       s(subjectID, bs="re"),
                     data = df_tract,
                     family = Gamma(link = "logit"),
                     method = "REML")
  
  gam.check(fit_cov_pds, rep = 500)
  
  # Test cov model against gamma
  # infoMessages('on')
  capture.output(compareML(fit_gamma, fit_cov_pds), file = paste0(dataDir, tract, "_GAM_comp.txt"))
  capture.output(summary(fit_cov_pds), file = paste0(dataDir, tract, "_GAM_cov.txt"))
  
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
  
  func_ggplot_gam(df_pred, "GAM Fit of R. Cingulate FA Values")
  
  
  # Check for group differences in spline.
  plot_diff <- func_plot_diff(fit_cov_pds)
  p01 <- plot_diff[[1]]
  p02 <- plot_diff[[2]]
  p12 <- plot_diff[[3]]
  
  
  
  # find biggest difference
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
  
  # negLGI x group
  fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
  # neuLGI x group
  fit <- lmList(NeuLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
}


# L. ATR
if(do_latr == 1){

  # Tract
  tract <- "ATR_L"
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$dti_fa <- round(df_tract$dti_fa, 3)
  
  #  determine k, compare families
  fit_gamma <- bam(dti_fa ~ Group +
                     Sex +
                     s(nodeID, by=Group, k=40) +
                     s(subjectID, bs="re"),
                   data = df_tract,
                   family = Gamma(link = "logit"),
                   method = "REML")
  
  gam.check(fit_gamma, rep = 500)
  
  fit_beta <- bam(dti_fa ~ Group +
                    Sex +
                    s(nodeID, by=Group, k=45) +
                    s(subjectID, bs="re"),
                  data = df_tract,
                  family = betar(link = "logit"),
                  method = "REML")
  
  gam.check(fit_beta, rep = 500)
  
  
  infoMessages('on')
  compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  
  # get stats
  summary(fit_gamma)  # R-sq = 0.864
  capture.output(summary(fit_gamma), file = paste0(dataDir, tract, "_GAM.txt"))
  
  # covariates
  fit_cov_pds <- bam(dti_fa ~ Group +
                       Sex +
                       s(nodeID, by=Group, k=40) +
                       s(PDS, by=Sex) +
                       s(subjectID, bs="re"),
                     data = df_tract,
                     family = Gamma(link = "logit"),
                     method = "REML")
  
  gam.check(fit_cov_pds, rep = 500)
  
  
  # Test cov model against gamma
  # infoMessages('on')
  capture.output(compareML(fit_gamma, fit_cov_pds), file = paste0(dataDir, tract, "_GAM_comp.txt"))
  capture.output(summary(fit_cov_pds), file = paste0(dataDir, tract, "_GAM_cov.txt"))
  
  
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
  
  func_ggplot_gam(df_pred, "GAM Fit of L. A. Thalamic Radiations FA Values")
  
  # Check for group differences in spline.
  plot_diff <- func_plot_diff(fit_cov_pds)
  p01 <- plot_diff[[1]]
  p02 <- plot_diff[[2]]
  p12 <- plot_diff[[3]]
  
  # find biggest difference
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
  
  # negLGI x group
  fit <- lmList(NegLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
  # neuLGI x group
  fit <- lmList(NeuLGI ~ dti_fa | Group, data = df_max)
  summary(fit)
  
}

