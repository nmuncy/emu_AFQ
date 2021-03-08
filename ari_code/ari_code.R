
# Intx
GroupSex <- df_tract$Sex

ind0 <- which(df_tract$Group==0 & df_tract$Sex==0)
ind1 <- which(df_tract$Group==0 & df_tract$Sex==1)
ind2 <- which(df_tract$Group==1 & df_tract$Sex==0)
ind3 <- which(df_tract$Group==1 & df_tract$Sex==1)
ind4 <- which(df_tract$Group==2 & df_tract$Sex==0)
ind5 <- which(df_tract$Group==2 & df_tract$Sex==1)

GroupSex[ind0] <- 0
GroupSex[ind1] <- 1
GroupSex[ind2] <- 2
GroupSex[ind3] <- 3
GroupSex[ind4] <- 4
GroupSex[ind5] <- 5

df_tract$GroupSex <- as.factor(GroupSex)

fit2 <- bam(dti_fa ~ Group + 
    Sex +
    s(nodeID,by=GroupSex,k=40) +
    s(PDS, by=Sex) +
    s(subjectID,bs="re")  ,
    data = df_tract,
    family = Gamma(link = "logit"),
    method = "REML")

plot_diff(fit2, view="nodeID", comp=list(GroupSex=c("0", "1")), rm.ranef = T)
plot_diff(fit2, view="nodeID", comp=list(GroupSex=c("2", "3")), rm.ranef = T, col='blue',add=T)
plot_diff(fit2, view="nodeID", comp=list(GroupSex=c("4", "5")), rm.ranef = T, col='red',add=T)


# Plots
df_pred <- predict.bam(
  fit_cov_pds,
  exclude_terms = c("PDS", "Sex","subjectID"),
  values=list(PDS = NULL, Sex = NULL),
  se.fit=T,
  type="response")

names(df_pred)
head(df_pred$fit)

df_pred <- data.frame(Group=df_tract$Group,
    Sex=df_tract$Sex,
    subjectID=df_tract$subjectID,
    PDS=df_tract$PDS,
    nodeID=df_tract$nodeID,
    fit=df_pred$fit,
    se.fit=df_pred$se.fit)

ggplot(data = df_pred) +
  geom_smooth(mapping = aes(x=nodeID, y=fit, color=Group))
