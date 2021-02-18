
# Orienting vars
afqDir <- "/Users/nmuncy/Projects/emu_diff/emu_data/derivatives/afq/"
subjList <- list.files(path=afqDir, pattern="sub-*")
sess <- "ses-S2"


### --- Step One
#
# Make master data frame in long form
# 
# Add group membership

df_master <- as.data.frame(matrix(NA, ncol=5, nrow=1))
colnames(df_master) <- c("Subj", "Tract", "Node", "FA", "MD")

for(subj in subjList){
  # subj <- as.character(subjList[1])
  
  # get data
  subjFile <- paste0(afqDir, subj, "/", sess, "/", subj, "_", sess, "_dwi_profiles.csv")
  df_subj <- read.csv(subjFile, sep = ",", header=T)
  
  # prepare, add to master
  colnames(df_subj) <- colnames(df_master)
  df_subj$Subj <- subj
  df_master <- rbind(df_master, df_subj)
  
  # make tract plots
  tractList <- unique(df_subj$tractID)
  for(tract in tractList){
    plot(df_subj[which(df_subj$tractID == tract), ]$dti_fa, type = "l", lwd=2)
    title(main=tract)
  }
}

# clean
df_master <- df_master[-1,]
row.names(df_master) <- 1:nrow(df_master)

# add anxiety
refFile <- read.delim("Group_Data.txt", sep="\t", header=T)
df_master$Group <- NA
for(subj in subjList){
  
  # determine subj indices
  ind_master <- which(df_master$Subj == subj)
  ind_anx <- which(refFile$Subj == gsub(".*-", "", subj))
  
  # determine group
  if(refFile[ind_anx,]$Pars7 > 3){
    h_group <- "Anxious"
  }else if(refFile[ind_anx,]$Pars7 <= 3){
    h_group <- "Control"
  }
  
  # fill
  df_master[ind_master,]$Group <- h_group
}






