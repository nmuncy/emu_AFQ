

### --- Notes:
#
# Quick stats for the manuscript

dataDir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"

# For group 2 (Con v GAD v SAD)
gType <- 2
tract <- "UNC_L"

dataFile <- paste0(dataDir, "Master_dataframe_G", gType,".csv")
df_afq <- read.csv(dataFile)
df_subset <- df_afq[which(df_afq$nodeID == 0 & df_afq$tractID == tract),]

num_subj <- dim(df_subset)[1]
num_female <- length(which(df_subset$Sex == 0))
age_avg <- round(mean(df_subset$Age),2)
age_sd <- round(sd(df_subset$Age),2)

stat_age <- aov(Age ~ Group, data = df_subset)
summary(stat_age)
stat_pds <- aov(PDS ~ Group, data = df_subset)
summary(stat_pds)
