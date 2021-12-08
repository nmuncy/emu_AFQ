library("ggplot2")
library("fitdistrplus")
library("mgcv")
library("itsadug")


# Notes ----
#
# This file contains code to address reviewer questions
#  and concerns, in addition to updates to step3. Namely
#  a "point-wise" of l_unc FA values has been added, 
#  and a GAM of l_unc MD values.


# Set Switches ----
stat_switch <- function(df_node, tract) {
  # Switch for using appropriate family with tracts.
  #
  # Defaults to Gamma family, but allows for other
  # tracts (FA) to use a GLM with another family.
  # Used for conducting traditional "point-wise"
  # comparisons of tract FA values at each node.
  #
  # Arguments:
  #   df_node = dataframe of data for single node
  #   tract = AFQ tract string
  #
  # Returns:
  #   h_stats = GLM model fit

  h_stats <- switch(tract,
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
  # Switch for decoding AFQ tract names.
  #
  # Arguments:
  #   tract = AFQ tract string
  #
  # Returns:
  #   x_tract = str, reformatted tract name

  x_tract <- switch(tract,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "FA" = "A. Forceps",
  )
  return(x_tract)
}

switch_plot_values <- function(value) {
  # Switch for determining group, coloring
  #
  # Used for keeping colors, labels consistent
  # across plots.
  #
  # Arguments:
  #   value = group factor (0-2), string
  #
  # Returns:
  #   x_col, x_label = indexed lists
  #     [[1]] = color, [[2]] = group label

  x_col <- switch(value,
    "0" = "blue",
    "1" = "darkred",
    "2" = "black"
  )

  x_label <- switch(value,
    "0" = "Low",
    "1" = "Med",
    "2" = "High"
  )

  return(list(x_col, x_label))
}


# Test Nodes ----
#
# Conduct a GLM at each node of each tract,
#   and FDR correct for number of node comparisons.
#   Then produce a plot of averaged FA values by
#   group, and indicate where tracts differ (Low vs High)
#   after an FDR correction.

data_dir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
out_dir <- paste0(data_dir, "traditional/")
data_file <- paste0(data_dir, "Master_dataframe_G3.csv")
tract_list <- c("UNC_L", "UNC_R", "FA")
df_afq <- as.data.frame(read.csv(data_file))

for (tract in tract_list) {

  # subset df_afq for tract
  df_tract <- df_afq[which(df_afq$tractID == tract), ]
  df_tract$sex <- factor(df_tract$sex)
  df_tract$group <- factor(df_tract$group)

  # get list of nodes
  node_list <- unique(df_tract$nodeID)

  # start empty stat table
  col_names <- c("Node", "Est.", "SE", "t.stat", "p.value", "LB", "UB")
  stat_table <- as.data.frame(
    matrix(NA, ncol = length(col_names), nrow = length(node_list))
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
    h_conf <- suppressMessages(round(confint(fit_glm)[3, ], 2))
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
      df_plot[ind_plot, ]$SE <- round(h_sd / (sqrt(num_group)), 4)
    }
  }

  # FDR adjust
  stat_table$FDR <- round(p.adjust(stat_table$p.value, method = "fdr"), 4)
  out_file <- paste0(out_dir, "stats_", tract, ".csv")
  write.table(
    stat_table,
    file = out_file, col.names = T, row.names = F, sep = ","
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

  # draw sig box, apply colors & labels - note that
  # this approach assumes sig_nodes are serial
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

  ggsave(
    paste0(out_dir, "Plot_", tract, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}



# Model MD ----
#
# Use the proposed workflow to investigate whether differences 
# exist between High and Lows PARS-6 groups in spline fits of 
# mean diffusivity.

md_dir <- paste0(data_dir, "model_md/")
tract <- "UNC_L"
df_afq$group <- factor(df_afq$group)
df_afq$sex <- factor(df_afq$sex)
df_tract <- df_afq[which(df_afq$tractID == tract), ]

# view data
ggplot(data = df_tract) +
  geom_smooth(mapping = aes(x = nodeID, y = dti_md, color = group))
descdist(df_tract$dti_md, discrete = F) # lognormal or gamma, gamma recommended

# model w/gamma
fit_gamma <- bam(dti_md ~ group +
  sex +
  s(nodeID, by = group, k = 10) +
  s(subjectID, bs = "re"),
data = df_tract,
family = Gamma(link = "logit"),
method = "REML"
)

# determine k
gam.check(fit_gamma, rep = 500)

capture.output(
  summary(fit_gamma),
  file = paste0(
    md_dir, "Stats_GAM-gamma_", tract, "_G3.txt"
  )
)

# model tract with covariates
fit_cov_pds <- bam(dti_md ~ group +
  sex +
  s(nodeID, by = group, k = 10) +
  s(pds, by = sex) +
  s(subjectID, bs = "re"),
data = df_tract,
family = Gamma(link = "logit"),
method = "REML"
)

# determine k
gam.check(fit_cov_pds, rep = 500)
capture.output(
  summary(fit_cov_pds),
  file = paste0(
    md_dir,
    "Stats_GAM-cov_",
    tract, "_G3.txt"
  )
)

# compare head-to-head
capture.output(
  compareML(fit_gamma, fit_cov_pds),
  file = paste0(
    md_dir,
    "Stats_GAM-comp_gam-cov_",
    tract, "_G3.txt"
  )
)

# generate predictions
df_pred <- predict.bam(
  fit_cov_pds,
  exclude_terms = c("pds", "sex", "subjectID"),
  values = list(pds = NULL, sex = NULL),
  se.fit = T,
  type = "response"
)

# convert predictions to dataframe
df_pred <- data.frame(
  Group = df_tract$group,
  sex = df_tract$sex,
  subjectID = df_tract$subjectID,
  pds = df_tract$pds,
  nodeID = df_tract$nodeID,
  fit = df_pred$fit,
  se.fit = df_pred$se.fit
)

# set up for plot
h_tract <- switch_tract_name(tract)
h_title <- paste0("GAM Fit of ", h_tract, " MD Values")

h_cols <- c(
  switch_plot_values("0")[[1]][1],
  switch_plot_values("1")[[1]][1],
  switch_plot_values("2")[[1]][1]
)
names(h_cols) <- c("0", "1", "2")
h_breaks <- c("0", "1", "2")
h_labels <- c(
  switch_plot_values("0")[[2]][1],
  switch_plot_values("1")[[2]][1],
  switch_plot_values("2")[[2]][1]
)

# draw plot
p <- ggplot(data = df_pred) +
  geom_smooth(mapping = aes(x = nodeID, y = fit, color = Group)) +
  ggtitle(h_title) +
  ylab("Fit MD") +
  xlab("Tract Node") +
  theme(text = element_text(
    family = "Times New Roman", face = "bold", size = 14
  ))

p + scale_color_manual(
  values = h_cols,
  breaks = h_breaks,
  labels = h_labels
)

ggsave(
  paste0(md_dir, "Plot_GAM_", tract, "_G3.png"),
  units = "in",
  width = 6,
  height = 6,
  device = "png"
)

# setup for plotting
factor_a <- "0"
factor_b <- "2"
group_a <- switch_plot_values(factor_a)[[2]][1]
group_b <- switch_plot_values(factor_b)[[2]][1]

# determine sig nodes
gam_model <- fit_cov_pds

p_summary <- capture.output(plot_diff(gam_model,
  view = "nodeID",
  comp = list(group = c(factor_a, factor_b)),
  rm.ranef = T
))
sig_regions <- p_summary[10:length(p_summary)]
sig_regions <- gsub("\\t", "", sig_regions)

# make list of start and end nodes, for shading
sig_list <- as.list(strsplit(sig_regions, " - "))
start_list <- as.numeric(sapply(sig_list, "[[", 1))
end_list <- as.numeric(sapply(sig_list, "[[", 2))

# determine bottom of plot
p_est <- plot_diff(gam_model,
  view = "nodeID",
  comp = list(group = c(factor_a, factor_b)),
  rm.ranef = T,
  plot = F
)
h_min <- min(p_est$est)
h_ci <- p_est[which(p_est$est == h_min), ]$CI
min_val <- h_min - h_ci

# set output
png(
  filename = paste0(
    md_dir, "Plot_Diff_", tract, "_G3_pair.png"
  ),
  width = 600, height = 600
)

# draw plot
par(mar = c(5, 5, 4, 2), family = "Times New Roman")
plot_diff(gam_model,
  view = "nodeID",
  comp = list(group = c(factor_a, factor_b)),
  rm.ranef = T,
  main = paste0(
    "Difference Scores, ", group_a, "-", group_b
  ),
  ylab = "Est. MD difference",
  xlab = "Tract Node",
  cex.lab = 2,
  cex.axis = 2,
  cex.main = 2,
  cex.sub = 1.5,
  col.diff = "red"
)

# shade significant regions
for (h_ind in 1:length(start_list)) {
  polygon(
    x = c(rep(start_list[h_ind], 2), rep(end_list[h_ind], 2)),
    y = c(0.035, min_val, min_val, 0.35),
    col = rgb(1, 0, 0, 0.2),
    border = NA
  )
}

par(mar = c(5, 4, 4, 2))
dev.off()
