

# Notes ----
#
# Something about traditional approach


# Setup ----
data_dir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
data_file <- paste0(data_dir, "Master_dataframe_G3.csv")
tract <- "UNC_L"

df_afq <- as.data.frame(read.csv(data_file))
df_tract <- df_afq[which(df_afq$tractID == tract), ]
df_tract$sex <- factor(df_tract$sex)
df_tract$group <- factor(df_tract$group)

p_capture <- vector()
node_list <- unique(df_tract$nodeID)
for (node in node_list) {

  # subset dataframe for data pertaining to node
  df_node <- df_tract[which(
    df_tract$nodeID == node
  ), ]

  # model with GLM, control for sex, sex by PDS interaction
  fit_glm <- glm(dti_fa ~ group + sex + sex * pds,
    family = gaussian(link = "logit"),
    data = df_node
  )

  # capture p-value of High (group2) vs Low (intercept)
  p_capture <- c(p_capture, round(coef(summary(fit_glm))[3, 4], 6))
}

# FDR adjust
fdr_list <- p.adjust(p_capture, method = "fdr")

# find nodes which differ at adjusted alpha
sig_nodes <- which(fdr_list < 0.05)
