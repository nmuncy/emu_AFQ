

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")


# Orienting vars
dataDir <- "/Users/nmuncy/Projects/emu_diff/analyses/"


### --- Step 1: Make dataset
#
# Add PDS, PARS scores
#
# Convert PARS to factors
#   Low: <= 3
#   Med: >3 & <=12
#   Hi: >12
#
# Writes analyses/Master_dataframe.csv

func_makeDF <- function(){
  
  # Get data
  df_afq <- read.delim(paste0(dataDir, "tract_profiles.csv"), sep = ",", header = T)
  colnames(df_afq) <- c("Counter", colnames(df_afq)[-1])
  df_full <- read.delim(paste0(dataDir, "emuR01_full_latest.csv"), sep = ",", header=T)
  df_pds <- read.delim(paste0(dataDir, "emuR01_pds_latest.csv"), sep = ",", header=T)
  
  # # Round FA/MDs
  # df_afq$dti_fa <- round(df_afq$dti_fa, 3)
  # df_afq$dti_md <- round(df_afq$dti_md, 3)
  
  # make lists
  subjList <- unique(df_afq$subjectID)
  tractList <- unique(df_afq$tractID)
  nodeList <- unique(df_afq$nodeID)
  
  # add group, pars6, pds
  df_afq$Group <- df_afq$Pars6 <- df_afq$PDS <- NA
  
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
    
    # fill
    df_afq[ind_afq,]$Group <- h_group
    df_afq[ind_afq,]$Pars6 <- h_anx
    df_afq[ind_afq,]$PDS <- h_pds
  }
  
  # convert Group to factor for mgcv
  df_afq$Group <- factor(df_afq$Group)
  # df_afq$PDS <- factor(df_afq$PDS)
  
  # write csv
  outFile <- paste0(dataDir, "Master_dataframe.csv")
  write.csv(df_afq, file=outFile, quote=F, row.names = F)
  return(df_afq)
}
df_afq <- func_makeDF()



### --- Step2: Tracts of Interest
#
# Determine distribution of data
#
# Find best GAM fit, determine K
#
# Use PARS as categorical, PDS as continuous
#
# Plot differences


# Get data
df_afq <- read.csv(paste0(dataDir, "Master_dataframe.csv"))
df_afq$Group <- factor(df_afq$Group)


### --- L. UNC
tract <- "UNC_L"
df_tract <- df_afq[which(df_afq$tractID == tract), ]
df_tract$dti_fa <- round(df_tract$dti_fa, 3)

# plot mean data
ggplot(data = df_tract) +
  geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

ggplot(data = df_tract) +
  geom_point(mapping = aes(x=nodeID, y=dti_fa, color=Group),size=0.3) +
  geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

# check shape
descdist(df_tract$dti_fa, discrete=F)

fit.beta <- fitdist(df_tract$dti_fa, "beta")
plot(fit.beta)
fit.beta$aic

fit.gamma <- fitdist(df_tract$dti_fa, "gamma")
plot(fit.gamma)
fit.gamma$aic


# compare different gam models
fit_gamma <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=35), data = df_tract, family = Gamma, method = "REML")
# gam.check(fit_gamma, rep = 500)
fit_beta <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# gam.check(fit_beta, rep = 500)
AIC(fit_gamma, fit_beta)


# plot better model
plot_smooths(model=fit_beta, series=nodeID, comparison=Group) + theme(legend.position = "top")
summary(g_beta)


# model with covariates
fit_cov <- gam(dti_fa ~ Group + PDS + s(nodeID, by = Group, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# gam.check(fit_cov, rep = 500)
summary(fit_cov)


# # plot all
# plot.gam(fit_cov, pages=1, ylab="Parameter Est.", all.terms=T)

# # plot single
# lunc_single <- getViz(fit_cov)
# p1 <- plot(sm(lunc_single, 1))
# p1 + l_fitLine(colour = "red") + 
#   l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
#   l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
#   l_points(shape = 19, size = 1, alpha = 0.1) + 
#   theme_classic()

# get pred controlling for PDS
df_pred <- predict_gam(fit_cov, exclude_terms = c("PDS"), values=list(PDS = NULL))
df_pred %>% 
  ggplot(aes(nodeID, fit)) +
  geom_smooth_ci(Group)

# # plot for e/cov
# plot_smooths(
#   model=fit_cov, 
#   series=nodeID, 
#   comparison=Group,
#   facet_terms = PDS
#   ) + theme(legend.position = "top")

# differences
# plot_difference(fit_beta, series=nodeID, difference = list(Group = c(0, 1, 2)))
# plot_difference(fit_beta, series=nodeID, difference = list(Group = c(0, 1)))
# plot_difference(fit_beta, series=nodeID, difference = list(Group = c(1, 2)))
# plot_difference(fit_beta, series=nodeID, difference = list(Group = c(0, 2)))

# test for differences between models
# anova(fit_beta, fit_cov, test="F")

infoMessages('on')
compareML(fit_beta, fit_cov)

par(mfrow=c(2,3))
plot_diff(fit_cov, view="nodeID", comp=list(Group=c(0, 1)), rm.ranef = T)
plot_diff(fit_cov, view="nodeID", comp=list(Group=c(0, 2)), rm.ranef = T)
plot_diff(fit_cov, view="nodeID", comp=list(Group=c(1, 2)), rm.ranef = T)

plot_diff2(
  fit_cov, 
  view=c("nodeID", "PDS"), 
  comp=list(Group=c(0, 1)), 
  rm.ranef = T, 
  zlim=c(-0.15, 0.06), 
  se=0, 
  show.diff = T)

plot_diff2(
  fit_cov, 
  view=c("nodeID", "PDS"), 
  comp=list(Group=c(0, 2)), 
  rm.ranef = T, 
  zlim=c(-0.15, 0.06), 
  se=0, 
  show.diff = T)

plot_diff2(
  fit_cov, 
  view=c("nodeID", "PDS"), 
  comp=list(Group=c(1, 2)), 
  rm.ranef = T, 
  zlim=c(-0.15, 0.06), 
  se=0, 
  show.diff = T)

par(mfrow=c(1,1))




### --- R. UNC
tract <- "UNC_R"
df_tract <- df_afq[which(df_afq$tractID == tract), ]
df_tract$dti_fa <- round(df_tract$dti_fa, 3)

# plot mean data
ggplot(data = df_tract) +
  geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

ggplot(data = df_tract) +
  geom_point(mapping = aes(x=nodeID, y=dti_fa, color=Group),size=0.3) +
  geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

# check shape
descdist(df_tract$dti_fa, discrete=F)

fit.beta <- fitdist(df_tract$dti_fa, "beta")
plot(fit.beta)
fit.beta$aic

fit.gamma <- fitdist(df_tract$dti_fa, "gamma")
plot(fit.gamma)
fit.gamma$aic

# compare different gam models
fit_gamma <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=35), data = df_tract, family = Gamma, method = "REML")
# gam.check(fit_gamma, rep = 500)
fit_beta <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# gam.check(fit_beta, rep = 500)
AIC(fit_gamma, fit_beta)

# plot better model
plot_smooths(model=fit_beta, series=nodeID, comparison=Group) + theme(legend.position = "top")
summary(g_beta)

# model with covariates
fit_cov <- gam(dti_fa ~ Group + PDS + s(nodeID, by = Group, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# gam.check(fit_cov, rep = 500)
summary(fit_cov)

# get pred controlling for PDS
df_pred <- predict_gam(fit_cov, exclude_terms = c("PDS"), values=list(PDS = NULL))
df_pred %>% 
  ggplot(aes(nodeID, fit)) +
  geom_smooth_ci(Group) + 
  ggtitle(paste0("Tract: ", tract))

# differences
infoMessages('on')
compareML(fit_beta, fit_cov)

par(mfrow=c(2,3))
plot_diff(fit_cov, view="nodeID", comp=list(Group=c(0, 1)), rm.ranef = T)
plot_diff(fit_cov, view="nodeID", comp=list(Group=c(0, 2)), rm.ranef = T)
plot_diff(fit_cov, view="nodeID", comp=list(Group=c(1, 2)), rm.ranef = T)

plot_diff2(
  fit_cov, 
  view=c("nodeID", "PDS"), 
  comp=list(Group=c(0, 1)), 
  rm.ranef = T, 
  zlim=c(-0.1, 0), 
  se=0, 
  show.diff = T)

plot_diff2(
  fit_cov, 
  view=c("nodeID", "PDS"), 
  comp=list(Group=c(0, 2)), 
  rm.ranef = T, 
  zlim=c(-0.15, 0.06), 
  se=0, 
  show.diff = T)

plot_diff2(
  fit_cov, 
  view=c("nodeID", "PDS"), 
  comp=list(Group=c(1, 2)), 
  rm.ranef = T, 
  zlim=c(-0.06, 0.06), 
  se=0, 
  show.diff = T)

par(mfrow=c(1,1))


