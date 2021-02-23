

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")


# Orienting vars
dataDir <- "/Users/nmuncy/Projects/emu_diff/analyses/"


### --- Make dataset

# Get data
df_afq <- read.delim(paste0(dataDir, "tract_profiles.csv"), sep = ",", header = T)
colnames(df_afq) <- c("Counter", colnames(df_afq)[-1])
df_full <- read.delim(paste0(dataDir, "emuR01_full_latest.csv"), sep = ",", header=T)
df_pds <- read.delim(paste0(dataDir, "emuR01_pds_latest.csv"), sep = ",", header=T)

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


# ### --- Run model for left uncinate
# 
# tract <- "UNC_L"
# df_tract <- df_afq[which(df_afq$tractID == tract), ]
# df_tract$dti_fa <- round(df_tract$dti_fa, 3)
# 
# # # plot mean data
# ggplot(data = df_tract) +
#   geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))
# 
# ggplot(data = df_tract) +
#   geom_point(mapping = aes(x=nodeID, y=dti_fa, color=Group),size=0.3) +
#   geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))
# 
# 
# # check shape
# descdist(df_tract$dti_fa, discrete=F)
# 
# fit.beta <- fitdist(df_tract$dti_fa, "beta")
# plot(fit.beta)
# fit.beta$aic
# 
# fit.gamma <- fitdist(df_tract$dti_fa, "gamma")
# plot(fit.gamma)
# fit.gamma$aic
# 
# 
# # separate model for each group, include diff intercepts
# #   each group gets own smoothness & wiggly parameters
# #   
# g_gam <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=35), data = df_tract, family = Gamma, method = "REML")
# gam.check(g_gam, rep = 500)
# 
# g_beta <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# gam.check(g_beta, rep = 500)
# 
# AIC(g_gam, g_beta)
# 
# plot_smooths(model=g_beta, series=nodeID, comparison=Group) + theme(legend.position = "top")
# summary(g_beta)
# 
# # # use single basis (fs) - likely the better model
# # #   assumes the groups will have a similar wiggliness,
# # #   closer to group as random effect
# # g <- gam(dti_fa ~ s(nodeID, Group, bs = 'fs'), data = df_tract, method = "REML")
# # plot(g,pages=1,scheme=1,unconditional=TRUE) 
# # summary(g)
# # gam.check(g, rep = 1000)
# 
# # differences
# plot_difference(g_beta, series=nodeID, difference = list(Group = c(0, 1, 2)))
# plot_difference(g_beta, series=nodeID, difference = list(Group = c(0, 1)))
# plot_difference(g_beta, series=nodeID, difference = list(Group = c(1, 2)))
# plot_difference(g_beta, series=nodeID, difference = list(Group = c(0, 2)))
# 
# 
# # covariate
# g_simp <- gam(dti_fa ~ s(nodeID, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# plot(g_simp,pages=1,scheme=1,unconditional=TRUE) 
# gam.check(g_simp, rep = 500)
# 
# # g_b <- gam(dti_fa ~ s(nodeID, k=35), data = df_tract, family = Gamma, method = "REML")
# # gam.check(g_b, rep = 500)
# # AIC(g_simp, g_b)
# 
# g_pars <- gam(dti_fa ~ Pars6 + s(nodeID, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# plot(g_pars,pages=1,scheme=1,unconditional=TRUE) 
# gam.check(g_pars, rep = 500)
# 
# g_pds <- gam(dti_fa ~ PDS + s(nodeID, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# plot(g_pds,pages=1,scheme=1,unconditional = T)
# gam.check(g_pds, rep=500)
# 
# g_pars_pds <- gam(dti_fa ~ PDS + Pars6 + s(nodeID, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# plot(g_pars_pds,pages=1,scheme=1,unconditional = T)
# gam.check(g_pars_pds, rep=500)
# 
# AIC(g_simp, g_beta)
# AIC(g_simp, g_pars)
# AIC(g_simp, g_pds)
# AIC(g_pars, g_pds)
# AIC(g_simp, g_pars_pds)


### --- Run model for right uncinate

tract <- "UNC_R"
df_tract <- df_afq[which(df_afq$tractID == tract), ]
# df_tract$dti_fa <- round(df_tract$dti_fa, 3)

# # plot mean data
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


# group model
g_gam <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=35), data = df_tract, family = Gamma, method = "REML")
gam.check(g_gam, rep = 500)

g_beta <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
gam.check(g_beta, rep = 500)

AIC(g_gam, g_beta)

plot_smooths(model=g_beta, series=nodeID, comparison=Group) + theme(legend.position = "top")
summary(g_beta)


# differences
plot_difference(g_beta, series=nodeID, difference = list(Group = c(0, 1, 2)))
plot_difference(g_beta, series=nodeID, difference = list(Group = c(0, 1)))
plot_difference(g_beta, series=nodeID, difference = list(Group = c(1, 2)))
plot_difference(g_beta, series=nodeID, difference = list(Group = c(0, 2)))


g_pds <- gam(dti_fa ~ Group + PDS + s(nodeID, by = Group, k=35), data = df_tract, family = betar(link = "logit"), method = "REML")
gam.check(g_pds, rep = 500)
# plot_smooths(
#   model=g_cov, 
#   series=nodeID, 
#   comparison=Group,
#   facet_terms=PDS) + theme(legend.position = "top")
# plot(g_cov,pages=1,scheme=0,unconditional=T)


# g_pds <- gam(dti_fa ~  PDS + s(nodeID, by=PDS, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
# gam.check(g_pds, rep = 500)
# plot_smooths(model=g_pds, series=nodeID, comparison=PDS) + theme(legend.position = "top")


