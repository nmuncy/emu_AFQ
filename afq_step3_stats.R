

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")
library("ez")
library("dplyr")
library("lme4")
# library("plotly")
# library("viridis")
# library("broom")


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
  df_adis <- read.delim(paste0(privateDir, "emuR01_adis.csv"), sep = ",", header=T)
  
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
    
    # if(h_anx <= 3){
    #   h_group <- 0
    # }else if(h_anx > 3 & h_anx < 13){
    #   h_group <- 1
    # }else if(h_anx > 12){
    #   h_group <- 2
    # }
    
    # determine group
    #   0 = con, 1 = anx, 2 = oth
    if(apply(
        df_adis[ind_adis,], 
        1, 
        function(r) any(r %in% "None")
      )){
      h_group <- 0
    }else if(apply(
        df_adis[ind_adis,], 
        1, 
        function(r) any(r %in% "Generalized Anxiety Disorder")
      )){
      h_group <- 1
    }else{
      h_group <- 2
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
        assign(check, 0.1)
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
  
  # write csv
  outFile <- paste0(dataDir, "Master_dataframe.csv")
  write.csv(df_afq, file=outFile, quote=F, row.names = F)
  return(df_afq)
}

func_memStats <- function(){
  
  ### --- Notes
  #
  # This function will conduct a
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
  
  capture.output(plot_diff(h_df,
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
                 file = paste0(outDir, "Table_Diff_", tract, "_01.txt"))
  
  par(mar=c(5,3,4,2))
  
  capture.output(plot_diff(h_df,
                     view="nodeID",
                     comp=list(Group=c("0", "2")),
                     rm.ranef = T,
                     main = "Difference Scores, Con-Oth",
                     ylab = "",
                     xlab = "Tract Node",
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2.5,
                     cex.sub = 2),
                 file = paste0(outDir, "Table_Diff_", tract, "_02.txt"))
  
  capture.output(plot_diff(h_df,
                     view="nodeID",
                     comp=list(Group=c("1", "2")),
                     rm.ranef = T,
                     main = "Difference Scores, Anx-Oth",
                     ylab = "",
                     xlab = "Tract Node",
                     cex.lab = 2,
                     cex.axis = 2,
                     cex.main = 2.5,
                     cex.sub = 2),
                 file = paste0(outDir, "Table_Diff_", tract, "_12.txt"))
  
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2))
  dev.off()
  
  # return(list(p01,p02,p12))
}

func_df_diff <- function(h_df){
  
  ### --- Notes:
  #
  # Return node of max diff per contrast
  
  p01 <- plot_diff(h_df,
                   view="nodeID",
                   comp=list(Group=c("0", "1")),
                   rm.ranef = T,
                   plot = F)
  m01 <- p01[which(abs(p01$est) == max(abs(p01$est))),]$nodeID
  
  p02 <- plot_diff(h_df,
                   view="nodeID",
                   comp=list(Group=c("0", "2")),
                   rm.ranef = T,
                   plot = F)
  m02 <- p02[which(abs(p02$est) == max(abs(p02$est))),]$nodeID
  
  p12 <- plot_diff(h_df,
                   view="nodeID",
                   comp=list(Group=c("1", "2")),
                   rm.ranef = T, 
                   plot = F)
  m12 <- p12[which(abs(p12$est) == max(abs(p12$est))),]$nodeID
  
  df_out <- as.data.frame(matrix(NA, nrow=3, ncol=2))
  colnames(df_out) <- c("Comparison", "Node")
  df_out[,1] <- c("01", "02", "12")
  df_out[,2] <- c(m01, m02, m12)
  return(df_out)
}

func_gam <- function(tract, df_tract, outDir){
  
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
  
  h_tract <- switch(
    tract,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "CGC_L" = "L. Cingulum",
    "CGC_R" = "R. Cingulum",
    "ATR_L" = "L. A. Thalamic Radiations",
    "ARC_L" = "L. Arcuate",
    "ARC_R" = "R. Arcuate"
  )
  
  plot_title = paste0("GAM Fit of ", h_tract," FA Values")
  func_ggplot_gam(df_pred, plot_title, outDir, tract)
  return(fit_cov_pds)
}
  
func_diff <- function(model, tract, outDir){
  
  ### Test for differences in GAM splines
  func_plot_diff(model, outDir, tract)

  # make table of sig regions
  df_out <- as.data.frame(matrix(NA, nrow=1, ncol=4))
  colnames(df_out) <- c("Comparison", "Section", "Start", "End")
  for(comp in c("01", "02", "12")){
    h_cmd = paste0("tail -n +10 ", dataDir, "Table_Diff_", tract, "_", comp, ".txt | sed 's/-/,/g'")
    h_lines <- system(h_cmd, intern = T)
    h_df <- read.table(text=paste(h_lines, collapse = "\n"), header = F, stringsAsFactors = F, sep = ",")
    for(i in 1:dim(h_df)[1]){
      df_out <- rbind(df_out, c(comp, i, h_df[i,1], h_df[i,2]))
    }
  }
  df_out <- na.omit(df_out)
  return(df_out)
}

func_dflm <- function(comp, df_tract, df_diff){
  
  grpA <- as.numeric(substr(comp, start=1, stop=1))
  grpB <- as.numeric(substr(comp, start=2, stop=2))
  df_comp <- df_diff[which(df_diff$Comparison == comp),]
  
  # make dataframe
  subjList <- unique(df_tract[which(df_tract$Group == grpA | df_tract$Group == grpB),]$subjectID)
  
  df_lm <- as.data.frame(matrix(NA, nrow=length(subjList), ncol=8))
  colnames(df_lm) <- c("Subj", "AvgFA", "Pars6", "Group", "NeuLGI", "NeuLDI", "NegLGI", "NegLDI")
  df_lm$Subj <- subjList
  
  for(subj in subjList){
    h_mean <- vector()
    for(i in 1:dim(df_comp)[1]){
      h_start <- which(df_tract$subjectID == subj & df_tract$nodeID == df_comp[i,]$Start)
      h_end <- which(df_tract$subjectID == subj & df_tract$nodeID == df_comp[i,]$End)
      h_mean <- c(h_mean, mean(df_tract[h_start:h_end,]$dti_fa))
    }
    ind_out <- which(df_lm$Subj == subj)
    df_lm[ind_out,]$AvgFA <- round(mean(h_mean), 4)
    
    ind_subj <- which(df_tract$subjectID == subj)[1]
    df_lm[ind_out,]$Pars6 <- df_tract[ind_subj,]$Pars6
    df_lm[ind_out,]$Group <- as.numeric(df_tract[ind_subj,]$Group)-1
    df_lm[ind_out,]$NeuLGI <- df_tract[ind_subj,]$NeuLGI
    df_lm[ind_out,]$NeuLDI <- df_tract[ind_subj,]$NeuLDI
    df_lm[ind_out,]$NegLGI <- df_tract[ind_subj,]$NegLGI
    df_lm[ind_out,]$NegLDI <- df_tract[ind_subj,]$NegLDI
  }
  df_lm$Group <- factor(df_lm$Group)
  return(df_lm)
}

func_dflm_max <- function(comp, df_tract, df_max){
  
  grpA <- as.numeric(substr(comp, start=1, stop=1))
  grpB <- as.numeric(substr(comp, start=2, stop=2))
  node <- df_max[which(df_max$Comparison == comp),]$Node
  
  # make dataframe
  subjList <- unique(df_tract[which(df_tract$Group == grpA | df_tract$Group == grpB),]$subjectID)
  df_lm <- as.data.frame(matrix(NA, nrow=length(subjList), ncol=8))
  colnames(df_lm) <- c("Subj", "MaxFA", "Pars6", "Group", "NeuLGI", "NeuLDI", "NegLGI", "NegLDI")
  df_lm$Subj <- subjList
  
  for(subj in subjList){
    
    ind_data <- which(df_tract$subjectID == subj & df_tract$nodeID == node)
    ind_out <- which(df_lm$Subj == subj)
    
    df_lm[ind_out,]$MaxFA <- df_tract[ind_data,]$dti_fa
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


### -- L Unc
# set up df
tract <- "UNC_L"
df_tract <- df_afq[which(df_afq$tractID == tract), ]
df_tract$dti_fa <- round(df_tract$dti_fa, 3)

# run gam, plot
gam_model <- func_gam(tract, df_tract, dataDir)
df_diff <- func_diff(gam_model, tract, dataDir)
df_max <- func_df_diff(gam_model)

# predict mem score from dti average diff
df_lm <- func_dflm("01", df_tract, df_diff)
fit.int <- lm(NegLGI ~ AvgFA*Group, data = df_lm)
summary(fit.int)
anova(fit.int)

ggplot(df_lm) +
  aes(x=AvgFA, y=NegLGI, shape=Group) +
  geom_point(aes(color=Group)) +
  geom_smooth(method = "lm")

fit.cov <- lm(NegLGI ~ AvgFA*Group + Pars6, data = df_lm)
summary(fit.cov)
anova(fit.cov)


# ggplot(df_lm, aes(x=AvgFA, y=NegLGI)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~ Group)

# df_long <- as.data.frame(matrix(NA, nrow=2*dim(df_lm)[1], ncol=4))
# colnames(df_long) <- c("Subj", "Group", "Meas", "Value")
# df_long$Subj <- rep(df_lm$Subj, 2)
# df_long$Group <- rep(df_lm$Group, 2)
# df_long$Meas <- c(rep("NegLGI", dim(df_lm)[1]), rep("AvgFA", dim(df_lm)[1]))
# df_long$Value <- c(df_lm$NegLGI, df_lm$AvgFA)

# predict mem score from dti max diff
df_lm <- func_dflm_max("01", df_tract, df_max)
fit.int <- lm(NegLGI ~ MaxFA*Group, data = df_lm)
summary(fit.int)
anova(fit.int)

ggplot(df_lm) +
  aes(x=MaxFA, y=NegLGI, shape=Group) +
  geom_point(aes(color=Group)) +
  geom_smooth(method = "lm")

fit.cov <- lm(NegLGI ~ MaxFA*Group + Pars6, data = df_lm)
summary(fit.cov)
anova(fit.cov)






