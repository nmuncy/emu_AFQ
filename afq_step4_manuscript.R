

library("lsr")

### --- Notes:
#
# Quick stats for the manuscript

dataDir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"

# For group 3 (L v M v H)
gType <- 3
tract <- "UNC_L"

dataFile <- paste0(dataDir, "Master_dataframe_G", gType,".csv")
df_afq <- read.csv(dataFile)
df_subset <- df_afq[which(df_afq$nodeID == 0 & df_afq$tractID == tract),]

num_subj <- dim(df_subset)[1]
num_female <- length(which(df_subset$Sex == 0))
age_avg <- round(mean(df_subset$Age),2)
age_sd <- round(sd(df_subset$Age),2)

stat_age <- aov(Age ~ as.factor(Group), data = df_subset)
summary(stat_age)
etaSquared(stat_age)

stat_pds <- aov(PDS ~ as.factor(Group), data = df_subset)
summary(stat_pds)
etaSquared(stat_pds)
