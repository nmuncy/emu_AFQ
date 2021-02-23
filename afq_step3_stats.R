

library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")


# Orienting vars
dataDir <- "/Users/nmuncy/Projects/emu_diff/analyses/"


### --- Make dataset

# Get group data
df_all <- read.delim(paste0(dataDir, "tract_profiles.csv"), sep = ",", header = T)
subjList <- unique(df_all$subjectID)
tractList <- unique(df_all$tractID)

# add anxiety
df_group <- read.delim(paste0(dataDir, "Group_data.txt"), sep = "\t", header=T)
df_all$Group <- NA
for(subj in subjList){

  # determine subj indices
  ind_master <- which(df_all$subjectID == subj)
  ind_anx <- which(df_group$Subj == subj)

  # determine group
  if(df_group[ind_anx,]$Pars7 > 13){
    h_group <- 1
  }else if(df_group[ind_anx,]$Pars7 <= 13){
    h_group <- 0
  }

  # fill
  df_all[ind_master,]$Group <- h_group
}
df_all$Group <- factor(df_all$Group)


### --- Run model for each tract

# subset
# for(tract in tractList){
#   
# }
tract <- "UNC_L"
df_tract <- df_all[which(df_all$tractID == tract), ]
df_tract$dti_fa <- round(df_tract$dti_fa, 3)

# # plot mean data
ggplot(data = df_tract) +
  geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

ggplot(data = df_tract) +
  geom_point(mapping = aes(x=nodeID, y=dti_fa, color=Group),size=0.3) +
  geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=Group))

# nodeList <- unique(df_tract$nodeID)
# df_plot <- as.data.frame(matrix(NA, nrow=2*length(nodeList), ncol=3))
# colnames(df_plot) <- c("Node", "Group", "FA")
# df_plot$Node <- rep(nodeList,2)
# df_plot$Group <- c(rep(0,length(nodeList)), rep(1,length(nodeList)))
# 
# for(i in 1:length(nodeList)){
#   node <- i-1
#   ind_con <- which(df_plot$Node == node & df_plot$Group == 0)
#   ind_anx <- which(df_plot$Node == node & df_plot$Group == 1)
#   df_plot[ind_con, ]$FA <- round(mean(df_tract[which(df_tract$nodeID == node & df_tract$Group == 0), ]$dti_fa), 3)
#   df_plot[ind_anx, ]$FA  <- round(mean(df_tract[which(df_tract$nodeID == node & df_tract$Group == 1), ]$dti_fa), 3)
# }
# 
# data_con <- df_plot[which(df_plot$Group == 0),]$FA
# data_anx <- df_plot[which(df_plot$Group == 1),]$FA
# 
# plot(data_con, type="l", col="blue")
# lines(data_anx, type="l", col="red")

# check shape
descdist(df_tract$dti_fa, discrete=F)

fit.beta <- fitdist(df_tract$dti_fa, "beta")
plot(fit.beta)
fit.beta$aic

fit.gamma <- fitdist(df_tract$dti_fa, "gamma")
plot(fit.gamma)
fit.gamma$aic


# separate model for each group, include diff intercepts
#   each group gets own smoothness & wiggly parameters
g <- gam(dti_fa ~ Group + s(nodeID, by = Group), data = df_tract, family = betar(link = "logit"), method = "REML")
plot(g,pages=1,scheme=1,unconditional=TRUE)
summary(g)
gam.check(g, rep = 500)

# # use single basis (fs) - likely the better model
# #   assumes the groups will have a similar wiggliness,
# #   closer to group as random effect
# g <- gam(dti_fa ~ s(nodeID, Group, bs = 'fs'), data = df_tract, method = "REML")
# plot(g,pages=1,scheme=1,unconditional=TRUE) 
# summary(g)
# gam.check(g, rep = 1000)



### --- Diagnose

# check K (k-index should be ~1)
g_k <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=30), data = df_tract, family = betar(link = "logit"), method = "REML")
gam.check(g_k, rep = 1000)
# plot(g_k,pages=1,scheme=1,unconditional=TRUE)
plot_smooths(model=g_k, series=nodeID, comparison=Group) + theme(legend.position = "top")
summary(g_k)

# check fam
g_f <- gam(dti_fa ~ Group + s(nodeID, by = Group, k=35), data = df_tract, family = Gamma, method = "REML")
gam.check(g_f, rep = 1000)
plot_smooths(model=g_f, series=nodeID, comparison=Group) + theme(legend.position = "top")
summary(g_f)

AIC(g_k, g_f)

# differences
plot_difference(g_k, series=nodeID, difference = list(Group = c(0, 1)))




