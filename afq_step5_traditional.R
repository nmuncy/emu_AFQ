library("ggplot2")

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

switch_tract_name <- function(tract) {
  # Switch for decoding AFQ tract names
  #
  # Arguments:
  #   tract = AFQ tract string
  #
  # Returns:
  #   x_tract = str, reformatted tract name
  
  x_tract <- switch(
    tract,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "FA" = "A. Forceps",
  )
  return(x_tract)
}


for (tract in tract_list) {
  
  # subset df_afq for tract
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$sex <- factor(df_tract$sex)
  df_tract$group <- factor(df_tract$group)
  
  # get list of nodes, start empty stat, plot tables
  node_list <- unique(df_tract$nodeID)
  
  # start empty stat table
  col_names <- c("Node", "Est.", "SE", "t.stat", "p.value", "LB", "UB")
  stat_table <- as.data.frame(
    matrix(NA, ncol=length(col_names), nrow=length(node_list))
  )
  colnames(stat_table) <- col_names
  stat_table$Node <- node_list
  
  # start empty plot df
  plot_names <- c("Node", "Group", "Avg", "SE")
  group_list <- unique(df_tract$group)
  df_plot <- as.data.frame(
    matrix(
      NA, 
      nrow = length(group_list) * length(node_list), 
      ncol = length(plot_names)
    )
  )
  colnames(df_plot) <- plot_names
  df_plot$Node <- rep(node_list, 3)
  df_plot$Group <- rep(group_list, each = length(node_list))
  
  # Do "point-wise" GLMs
  for (node in node_list) {
    
    # subset dataframe for data pertaining to node
    df_node <- df_tract[which(
      df_tract$nodeID == node
    ), ]
    
    # model with GLM, use correct family for each tract
    fit_glm <- stat_switch(df_node, tract)
    
    # get H vs L coefs, conf intervals, fill table
    h_coef <- round(coef(summary(fit_glm))[3, 4], 4)
    h_conf <- suppressMessages(round(confint(fit_glm)[3,], 2))
    ind_node <- which(stat_table$Node == node)
    stat_table[ind_node, 2:7] <- c(h_coef, h_conf)
    
    # get avg, se for plotting
    for (group in group_list) {
      num_group <- length(which(df_node$group == group))
      ind_plot <- which(df_plot$Node == node & df_plot$Group == group)
      df_plot[ind_plot, ]$Avg <- round(
        mean(df_node[which(df_node$group == group), ]$dti_fa), 4
      )
      h_sd <- sd(df_node[which(df_node$group == group), ]$dti_fa)
      df_plot[ind_plot, ]$SE <- round(h_sd/(sqrt(num_group)), 4)
    }
    
  }
  
  # FDR adjust
  stat_table$FDR <- round(p.adjust(stat_table$p.value, method = "fdr"), 4)
  out_file <- paste0(out_dir, "stats_", tract, ".csv")
  write.table(
    stat_table, file = out_file, col.names = T, row.names = F, sep = ","
  )
  
  # set plot standard colors
  h_cols <- c("blue", "darkred", "black")
  names(h_cols) <- c("0", "1", "2")
  
  # plot lines, confidence
  h_title <- paste(switch_tract_name(tract), "Tract FA Values")
  p <- ggplot(data = df_plot, aes(x = Node, y = Avg, group = Group)) +
    geom_line(aes(color = Group)) +
    geom_ribbon(aes(ymin = Avg - SE, ymax = Avg + SE), alpha = 0.2) +
    ggtitle(h_title) +
    ylab("Avg FA") +
    xlab("Tract Node") +
    theme(text = element_text(
      family = "Times New Roman", face = "bold", size = 14
    ))
  
  # draw sig box, apply colors & labels
  sig_nodes <- which(stat_table$FDR < 0.05)
  if (length(sig_nodes) > 0) {
    p + annotate(
      "rect", 
      xmin = sig_nodes[1], 
      xmax = sig_nodes[length(sig_nodes)],
      ymin = min(df_plot$Avg),
      ymax = max(df_plot$Avg) + 0.025,
      alpha = 0.1,
      fill = "red"
    ) + scale_color_manual(
      values = h_cols,
      breaks = c("0", "1", "2"),
      labels = c("Low", "Med", "High")
    )
  } else {
    p + scale_color_manual(
      values = h_cols,
      breaks = c("0", "1", "2"),
      labels = c("Low", "Med", "High")
    )
  }
  
  
  
  
  

}





