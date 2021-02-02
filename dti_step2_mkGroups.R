data_dir <- getwd()

# import, reduce data
df_raw <- read.csv(paste0(data_dir, "/emuR01_pars_latest.csv"))
df_reduced <- as.data.frame(matrix(NA, nrow=dim(df_raw)[1], ncol=5))
colnames(df_reduced) <- c("Subj", "Sleep", "Anxious", "Pars6", "Pars7")
df_reduced[,1:5] <- c(df_raw$src_subject_id, df_raw$pinf_random, df_raw$pinf_group, df_raw$pars_6, df_raw$pars_7)

# write
out_file <- paste0(data_dir, "/Group_data.csv")
write.table(df_reduced, out_file, quote = F, sep = ",", row.names = F)