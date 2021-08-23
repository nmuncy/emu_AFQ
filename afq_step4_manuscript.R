library("lsr")


# Notes:
#
# Quick stats for the manuscript


# Set Up -------
data_dir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
tract <- "UNC_L"

# Use group 3 (L v M v H) data
g_type <- 3

# read in data, subset
data_file <- paste0(data_dir, "Master_dataframe_G", g_type, ".csv")
df_afq <- read.csv(data_file)
df_subset <- df_afq[which(df_afq$nodeID == 0 & df_afq$tractID == tract), ]


# Get Demographics ------
num_subj <- dim(df_subset)[1]
num_female <- length(which(df_subset$Sex == 0))
age_avg <- round(mean(df_subset$Age), 2)
age_sd <- round(sd(df_subset$Age), 2)

stat_age <- aov(Age ~ as.factor(Group), data = df_subset)
summary(stat_age)
etaSquared(stat_age)

stat_pds <- aov(PDS ~ as.factor(Group), data = df_subset)
summary(stat_pds)
etaSquared(stat_pds)
