

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


### --- Step 1: Make dataset
#
# Add PDS, PARS, d-prime scores
#   add sex, age
#
# Convert PARS to factors
#   Low: <= 3
#   Med: >3 & <=12
#   Hi: >12
#
# Sex: 0 = female, 1 = male
#
# Writes analyses/Master_dataframe.csv

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
#
# Run ANOVA for LGI, LDI
#
# Writes to analyses/Stats_AN_L?I.txt
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



### --- Step 3: Model tract, no covariates
#
# Plot data, determine distribution,
#   compare model families.
#
# Plot best model

# Get data
df_afq <- read.csv(paste0(dataDir, "Master_dataframe.csv"))
df_afq$Group <- factor(df_afq$Group)
df_afq$Sex <- factor(df_afq$Sex)

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
compareML(fit_gamma, fit_cov_pds) # PDS wins
summary(fit_cov_pds) # R-sq = 0.85


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

ggplot(data = df_pred) +
  geom_smooth(mapping = aes(x=nodeID, y=fit, color=Group)) +
  ggtitle("GAM Fit of L. Uncinate FA Values") +
  ylab("Fit FA")



### --- Step 5: Test for differences
#
# Check for group differences in spline.

par(mfrow=c(1,3))

p01 <- plot_diff(fit_cov_pds,
    view="nodeID",
    comp=list(Group=c("0", "1")),
    rm.ranef = T,
    main = "Difference Scores, Low-Med",
    ylab = "Est. FA difference",
    xlab = "Tract Node")

p02 <- plot_diff(fit_cov_pds,
    view="nodeID",
    comp=list(Group=c("0", "2")),
    rm.ranef = T,
    main = "Difference Scores, Low-High",
    ylab = "Est. FA difference",
    xlab = "Tract Node")

p12 <- plot_diff(fit_cov_pds,
    view="nodeID",
    comp=list(Group=c("1", "2")),
    rm.ranef = T,
    main = "Difference Scores, Med-High",
    ylab = "Est. FA difference",
    xlab = "Tract Node")

par(mfrow=c(1,1))

# get location of maximum difference
max_node_01 <- as.numeric(which(abs(p01$est) == max(abs(p01$est)))) - 1
max_node_02 <- as.numeric(which(abs(p02$est) == max(abs(p02$est)))) - 1
max_node_12 <- as.numeric(which(abs(p12$est) == max(abs(p12$est)))) - 1



### --- Step 6: Regress
#
# regress node FA value on behavior

df_max <- as.data.frame(df_tract[which(
  df_tract$nodeID == max_node_01 |
  df_tract$nodeID == max_node_02 |
  df_tract$nodeID == max_node_12
),])
df_reduce <- as.data.frame(df_max[,-c(1:2, 5, 7:10, 14:16)])
df_reduce$nodeID <- factor(df_reduce$nodeID)

ind_pos02 <- which(df_reduce$nodeID == max_node_02 &
   (df_reduce$Group == 0 | df_reduce$Group == 2))
df_pos02 <- as.data.frame(df_reduce[ind_pos02,])

# negLGI x group
ggplot(df_pos02, aes(x=dti_fa, y=NegLGI)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Group)

# xyplot(NegLGI ~ dti_fa, groups=Group, data=df_pos02, type='l')
fit <- lmList(NegLGI ~ dti_fa | Group, data = df_pos02)
summary(fit)




