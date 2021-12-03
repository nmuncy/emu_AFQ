

# Notes ----
#
# Something about traditional approach


# Setup ----
data_dir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
out_dir <- paste0(data_dir, "traditional/")
data_file <- paste0(data_dir, "Master_dataframe_G3.csv")
tract_list <- c("UNC_L", "UNC_R", "FA")
df_afq <- as.data.frame(read.csv(data_file))


# Test Nodes ----
#
# Conduct a GLM at each node of each tract,
# and FDR correct for each tract.

stat_switch <- function(df_node, tract){
  h_stats <- switch(
    tract,
    "FA" = glm(dti_fa ~ group + sex + sex * pds,
               family = gaussian(link = "logit"),
               data = df_node
    ),
    glm(dti_fa ~ group + sex + sex * pds,
                  family = Gamma(link = "logit"),
                  data = df_node
    )
  )
  return(h_stats)
}


for (tract in tract_list) {
  
  # subset df_afq for tract
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$sex <- factor(df_tract$sex)
  df_tract$group <- factor(df_tract$group)
  
  # start empty vector, get list of nodes
  node_list <- unique(df_tract$nodeID)
  col_names <- c("Node", "Est.", "SE", "t.stat", "p.value", "LB", "UB")
  stat_table <- as.data.frame(
    matrix(NA, ncol=length(col_names), nrow=length(node_list))
  )
  colnames(stat_table) <- col_names
  stat_table$Node <- node_list
  
  # Do "point-wise" GLMs
  for (node in node_list) {
    
    # subset dataframe for data pertaining to node
    df_node <- df_tract[which(
      df_tract$nodeID == node
    ), ]
    
    # model with GLM, control for sex, sex by PDS interaction
    fit_glm <- stat_switch(df_node, tract)
    
    # get H vs L coefs, conf intervals, fill table
    h_coef <- round(coef(summary(fit_glm))[3, 4], 4)
    h_conf <- suppressMessages(round(confint(fit_glm)[3,], 2))
    ind_node <- which(stat_table$Node == node)
    stat_table[ind_node, 2:7] <- c(h_coef, h_conf)
  }
  
  # FDR adjust
  stat_table$FDR <- round(p.adjust(stat_table$p.value, method = "fdr"), 4)
  out_file <- paste0(out_dir, "stats_", tract, ".csv")
  write.table(stat_table, file = out_file, col.names = T, row.names = F, sep = ",")
}

