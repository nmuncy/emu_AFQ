

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")


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
# GroupSex: 0-2 = F-LMH, 3-5 = M-LMH
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
  df_afq$NegDP12h <- df_afq$NegDP1wk <- df_afq$NeuDP12h <- df_afq$NeuDP1wk <- NA
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

    # determine d-primes
    h_neg12 <- df_full[ind_full,]$dprime_neg_12H
    h_neu12 <- df_full[ind_full,]$dprime_neu_12H
    h_neg1wk <- df_full[ind_full,]$dprime_neg_1WK
    h_neu1wk <- df_full[ind_full,]$dprime_neu_1WK

    # fill
    df_afq[ind_afq,]$Group <- h_group
    df_afq[ind_afq,]$Pars6 <- h_anx
    df_afq[ind_afq,]$PDS <- h_pds
    df_afq[ind_afq,]$Age <- h_age
    df_afq[ind_afq,]$Sex <- h_sexF
    df_afq[ind_afq,]$NegDP12h <- h_neg12
    df_afq[ind_afq,]$NeuDP12h <- h_neu12
    df_afq[ind_afq,]$NegDP1wk <- h_neg1wk
    df_afq[ind_afq,]$NeuDP1wk <- h_neu1wk
  }

  # add GroupSex
  GroupSex <- df_afq$Sex

  ind0 <- which(df_afq$Group==0 & df_afq$Sex==0)
  ind1 <- which(df_afq$Group==1 & df_afq$Sex==0)
  ind2 <- which(df_afq$Group==2 & df_afq$Sex==0)
  ind3 <- which(df_afq$Group==0 & df_afq$Sex==1)
  ind4 <- which(df_afq$Group==1 & df_afq$Sex==1)
  ind5 <- which(df_afq$Group==2 & df_afq$Sex==1)

  GroupSex[ind0] <- 0
  GroupSex[ind1] <- 1
  GroupSex[ind2] <- 2
  GroupSex[ind3] <- 3
  GroupSex[ind4] <- 4
  GroupSex[ind5] <- 5

  df_afq$GroupSex <- as.factor(GroupSex)


  # write csv
  outFile <- paste0(dataDir, "Master_dataframe.csv")
  write.csv(df_afq, file=outFile, quote=F, row.names = F)
  return(df_afq)
}
df_afq <- func_makeDF()




### --- Step 2: Model tract, no covariates
#
# Plot data, determine distribution,
#   compare model families.
#
# Plot best model

# Get data
df_afq <- read.csv(paste0(dataDir, "Master_dataframe.csv"))
df_afq$Group <- factor(df_afq$Group)
df_afq$Sex <- factor(df_afq$Sex)
df_afq$GroupSex <- factor(df_afq$GroupSex)

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



### --- Step 3: Model tract, covariates

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
compareML(fit_gamma, fit_cov_pds) #PDS wins
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
  ggtitle("Tract: L. Uncinate") +
  ylab("Fit FA")



### --- Step 4: Test for differences
#
# Check for group differences in spline.

# par(mfrow=c(2,3))
# par(mar=c(5, 4, 4, 4))

# Male vs female on diff groups
p01 <- plot_diff(fit_cov_pds,
    view="nodeID",
    comp=list(Group=c("0", "1")),
    rm.ranef = T)

p02 <- plot_diff(fit_cov_pds,
    view="nodeID",
    comp=list(Group=c("0", "2")),
    rm.ranef = T)

p12 <- plot_diff(fit_cov_pds,
    view="nodeID",
    comp=list(Group=c("1", "2")),
    rm.ranef = T)

# plot_diff2(
#   fit_cov_pds,
#   view=c("nodeID", "PDS"),
#   comp=list(Group=c(0, 1)),
#   rm.ranef = T,
#   zlim=c(-0.12, 0.02),
#   dec=2,
#   se=0,
#   color="topo",
#   show.diff = T)


# par(mfrow=c(1,1))
# par(mar=c(5.1, 4.1, 4.1, 2.1))

# get location of maximum difference
max_node_01 <- as.numeric(which(abs(p01$est) == max(abs(p01$est)))) - 1
max_node_02 <- as.numeric(which(abs(p02$est) == max(abs(p02$est)))) - 1
max_node_12 <- as.numeric(which(abs(p12$est) == max(abs(p12$est)))) - 1

