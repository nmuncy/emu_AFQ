library("mgcv")
library("tidymv")
library("fitdistrplus")
library("ggplot2")
library("itsadug")
library("mgcViz")
library("ez")
library("dplyr")
library("lme4")
library("ggpubr")
library("lsr")


# Notes:
#
# This script reads in pyAFQ data and conducts a series of 
#   GAMs and pre-determined tracts.
#
# Group differences in spline fits are then investigated,
#   and the average of differential nodes are then used
#   to investigate LGI performance.
#
# Potential groupings:
#   1) Control vs Anxious
#   2) Low vs Med vs High Pars6 Score (currently used)
#   3) Control vs GAD vs SAD
#
# For testing group spline differences, one pair can be
#   tested with a single test ("pair"), or all pair 
#   permutations can be compaired ("anova"). "Pair" is
#   currently used.


# Set Up ------
#
# Various input, output paths/directories

data_dir <- "/Users/nmuncy/Projects/emu_AFQ/analyses/"
private_dir <- "/Users/nmuncy/Projects/emu_data/emu_private/"

gam_plot_dir <- paste0(data_dir, "plots_gam/")
gam_stats_dir <- paste0(data_dir, "stats_gam/")
lm_plot_dir <- paste0(data_dir, "plots_lm/")
lm_stats_dir <- paste0(data_dir, "stats_lm/")
table_dir <- paste0(data_dir, "tables/")
memory_dir <- paste0(data_dir, "memory/")


# General Functions ------
make_master_df <- function(g_type) {

  ### --- Notes:
  #
  # This function will make a master dataframe that
  # contains node FA values, group membership, sex, pds
  # and memory measures (LG/DI)
  #
  # Writes data_dir/Master_dataframe.csv

  # Get data
  df_afq <- read.delim(
    paste0(data_dir, "tract_profiles.csv"),
    sep = ",", header = T
  )
  colnames(df_afq) <- c("Counter", colnames(df_afq)[-1])

  df_full <- read.delim(
    paste0(private_dir, "emuR01_full_latest.csv"),
    sep = ",", header = T
  )

  df_pds <- read.delim(
    paste0(private_dir, "emuR01_pds_latest.csv"),
    sep = ",", header = T
  )

  df_adis <- read.delim(
    paste0(private_dir, "emuR01_adis.csv"),
    sep = ",", header = T
  )

  # make lists
  subj_list <- unique(df_afq$subjectID)
  tract_list <- unique(df_afq$tractID)
  node_list <- unique(df_afq$nodeID)

  # add group (Adis), pars6, pds
  #   add LG/DI, sex, age
  df_afq$sex <- df_afq$age <- NA
  df_afq$group <- df_afq$pars6 <- df_afq$pds <- NA
  df_afq$neg_ldi <- df_afq$neu_ldi <- df_afq$neg_lgi <- df_afq$neu_lgi <- NA

  for (subj in subj_list) {

    # determine subj indices
    ind_afq <- which(df_afq$subjectID == subj)
    ind_full <- which(df_full$src_subject_id == subj)
    ind_pds <- which(df_pds$emu_study_id == subj)
    ind_adis <- which(df_adis$Participant.ID == subj)

    # get pars
    h_anx <- df_full[ind_full, ]$pars_6

    # Determine group in one of three ways:
    #   a) 0 = con, 1 = anx
    #       skip subj when not anx/phobia/control
    #       (not used)
    #   b) 0 = con, 1 = gad, 2 = social/separation
    #       1 = GAD in dx.1, or dx GAD but SAD not dx.1
    #   c) 0 = low, 1 = med, 2 = high pars
    #       low=0-3, med=4-12, high>12
    if (g_type == 1) {
      h_search <- c("Anxiety", "Phobia")

      if (
        sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis, ])) != 0
      ) {
        h_group <- 1
      } else if (length(grep("None", df_adis[ind_adis, ])) != 0) {
        h_group <- 0
      } else {
        next
      }
    } else if (g_type == 2) {
      h_search <- c("Separation", "Social")

      if (
        grepl("Gen", df_adis[ind_adis, ]$Diagnosis.1) == T |
          (sum(grep("Gen", df_adis[ind_adis, ])) != 0 &
            sum(
              grep(paste(h_search, collapse = "|"), df_adis[ind_adis, ])
            ) == 0
          )
      ) {
        h_group <- 1
      } else if (
        sum(grep(paste(h_search, collapse = "|"), df_adis[ind_adis, ])) != 0
      ) {
        h_group <- 2
      } else if (length(grep("None", df_adis[ind_adis, ])) != 0) {
        h_group <- 0
      } else {
        next
      }
    } else if (g_type == 3) {
      if (h_anx <= 3) {
        h_group <- 0
      } else if (h_anx > 3 & h_anx < 13) {
        h_group <- 1
      } else if (h_anx > 12) {
        h_group <- 2
      }
    }

    # get pds
    h_pds <- df_pds[ind_pds, ]$pds_shirtcliff

    # get age, sex
    h_age <- df_full[ind_full, ]$pinf_age
    h_sex <- substr(df_full[ind_full, ]$sex, 1, 1)
    if (h_sex == "f") {
      h_sex_f <- 0
    } else if (h_sex == "m") {
      h_sex_f <- 1
    }

    # get Beh counts
    neg_num_hit <- df_full[ind_full, ]$negtarght_cnt_1WK
    neg_num_miss <- df_full[ind_full, ]$negtargms_cnt_1WK
    neg_num_lcr <- df_full[ind_full, ]$neglurecr_cnt_1WK
    neg_num_lfa <- df_full[ind_full, ]$neglurefa_cnt_1WK
    neg_num_fcr <- df_full[ind_full, ]$negfoilcr_cnt_1WK
    neg_num_ffa <- df_full[ind_full, ]$negfoilfa_cnt_1WK

    neu_num_hit <- df_full[ind_full, ]$neutarght_cnt_1WK
    neu_num_miss <- df_full[ind_full, ]$neutargms_cnt_1WK
    neu_num_lcr <- df_full[ind_full, ]$neulurecr_cnt_1WK
    neu_num_lfa <- df_full[ind_full, ]$neulurefa_cnt_1WK
    neu_num_fcr <- df_full[ind_full, ]$neufoilcr_cnt_1WK
    neu_num_ffa <- df_full[ind_full, ]$neufoilfa_cnt_1WK

    # adjust 0 counts
    for (check in c(
      "neg_num_hit",
      "neg_num_miss",
      "neg_num_lcr",
      "neg_num_lfa",
      "neg_num_fcr",
      "neg_num_ffa",
      "neu_num_hit",
      "neu_num_miss",
      "neu_num_lcr",
      "neu_num_lfa",
      "neu_num_fcr",
      "neu_num_ffa"
    )) {
      check_val <- get(check)
      if (check_val == 0) {
        assign(check, 0.001)
      }
    }

    # Calculate Neg/Neu LD/GI
    #   ldi = p(N|L) - p(N|T)
    #   lgi = p(O|L) - p(O|F)
    neg_ldi <- round(
      (neg_num_lcr / (neg_num_lcr + neg_num_lfa)) -
        (neg_num_miss / (neg_num_miss + neg_num_hit)),
      2
    )

    neu_ldi <- round(
      (neu_num_lcr / (neu_num_lcr + neu_num_lfa)) -
        (neu_num_miss / (neu_num_miss + neu_num_hit)),
      2
    )

    neg_lgi <- round(
      (neg_num_lfa / (neg_num_lcr + neg_num_lfa)) -
        (neg_num_ffa / (neg_num_ffa + neg_num_fcr)),
      2
    )

    neu_lgi <- round(
      (neu_num_lfa / (neu_num_lcr + neu_num_lfa)) -
        (neu_num_ffa / (neu_num_ffa + neu_num_fcr)),
      2
    )

    # fill
    df_afq[ind_afq, ]$group <- h_group
    df_afq[ind_afq, ]$pars6 <- h_anx
    df_afq[ind_afq, ]$pds <- h_pds
    df_afq[ind_afq, ]$age <- h_age
    df_afq[ind_afq, ]$sex <- h_sex_f

    df_afq[ind_afq, ]$neg_ldi <- neg_ldi
    df_afq[ind_afq, ]$neu_ldi <- neu_ldi
    df_afq[ind_afq, ]$neg_lgi <- neg_lgi
    df_afq[ind_afq, ]$neu_lgi <- neu_lgi
  }

  # clean NA (from group skip), write csv
  df_out <- df_afq[complete.cases(df_afq$group), ]
  out_file <- paste0(data_dir, "Master_dataframe_G", g_type, ".csv")
  write.csv(df_out, file = out_file, quote = F, row.names = F)
  # return(df_out)
}

switch_plot_values <- function(value, g_type) {

  ### --- Notes:
  #
  # Switch for determining group, coloring

  x_col <- switch(
    value,
    "0" = "blue",
    "1" = "darkred",
    "2" = "black"
  )

  if (g_type == 2) {
    x_label <- switch(
      value,
      "0" = "Con",
      "1" = "GAD",
      "2" = "SAD"
    )
  } else if (g_type == 3) {
    x_label <- switch(
      value,
      "0" = "Low",
      "1" = "Med",
      "2" = "High"
    )
  }
  return(list(x_col, x_label))
}

switch_tract_name <- function(tract) {

  ### --- Notes:
  #
  # Switch for determining name
  #   of tracts

  x_tract <- switch(
    tract,
    "UNC_L" = "L. Uncinate",
    "UNC_R" = "R. Uncinate",
    "FA" = "A. Forceps",
  )
  return(x_tract)
}

calc_memory_stats <- function(df_afq, g_type) {

  ### --- Notes:
  #
  # This function will conduct a MANOVA
  #   and post-hoc on memory measures.

  # subset df for ease of use
  ind_mem <- which(df_afq$tractID == tract_list[1] & df_afq$nodeID == 0)
  df_mem <- as.data.frame(df_afq[ind_mem, ])
  subj_list <- df_mem$subjectID

  # DF for group * Neg/Neu Detect/Generalize indices
  #   update, only do lgi
  df_long <- as.data.frame(matrix(NA, nrow = 4 * length(subj_list), ncol = 5))
  colnames(df_long) <- c("subj", "group", "valence", "memory", "value")

  df_long$subj <- rep(subj_list, 4)
  df_long$group <- rep(df_mem$group, 4)
  df_long$valence <- c(
    rep("neg", 2 * length(subj_list)),
    rep("neu", 2 * length(subj_list))
  )
  df_long$memory <- c(
    rep("lgi", length(subj_list)),
    rep("ldi", length(subj_list)),
    rep("lgi", length(subj_list)),
    rep("ldi", length(subj_list))
  )
  df_long$value <- c(
    df_mem$neg_lgi,
    df_mem$neg_ldi,
    df_mem$neu_lgi,
    df_mem$neu_ldi
  )

  # stats
  stats_mem <- ezANOVA(df_long,
    dv = value,
    wid = subj,
    between = group,
    within = c(valence, memory),
    type = "III"
  )

  # write
  capture.output(
    stats_mem,
    file = paste0(memory_dir, "Stats_G", g_type, "_MANOVA.txt")
  )

  # # post-hoc
  # intx_p <- stats_mem$`Sphericity Corrections`$`p[GG]`[2]
  # if(intx_p < 0.05){
  #   h.man <- manova(
  #     cbind(neu_lgi, neg_lgi, neu_ldi, neg_ldi) ~ group,
  #     data = df_mem
  #   )
  #   capture.output(
  #     summary.aov(h.man),
  #     file = paste0(memory_dir, "Stats_G", g_type, "_post.txt")
  #   )
  # }
  #
  # # Tukey for g3 neg_lgi
  # if(g_type == 3){
  #
  #   group_fac <- as.numeric(df_mem$group)
  #   group_value <- vector()
  #   for(val in group_fac){
  #     group_value <- c(group_value, switch_plot_values(val, g_type)[[2]])
  #   }
  #   df_mem$group <- group_value
  #
  #   h.lm <- lm(neg_lgi ~ group, data = df_mem)
  #   h.av <- aov(h.lm)
  #   capture.output(
  #     TukeyHSD(h.av),
  #     file = paste0(memory_dir, "Stats_G", g_type, "_tuk.txt")
  #   )
  #
  #   # make a plot
  #   h_colors <- c(switch_plot_values("2", g_type)[[1]],
  #     switch_plot_values("0", g_type)[[1]],
  #     switch_plot_values("1", g_type)[[1]]
  #   )
  #   post_comp <- list(c("High", "Low"))
  #   ggboxplot(df_mem,
  #     y="neg_lgi",
  #     x="group",
  #     color="group",
  #     palette = h_colors,
  #     add = "jitter",
  #     title = "top"
  #   ) +
  #     stat_compare_means(comparisons = post_comp) +
  #     ggtitle("Memory Performance on Negative Stimuli") +
  #     rremove("legend") +
  #     theme(text=element_text(family="Times New Roman", face="bold", size=16))
  #
  #   ggsave(
  #     paste0(memory_dir, "Plot_Box_G", g_type, "_neg_lgi.tiff"),
  #     units = "in",
  #     width = 6,
  #     height = 6,
  #     device = "tiff"
  #   )
  # }
}


# GAM Functions ------
help_plot_gam <- function(model, tract, g_type, df_tract) {

  ### --- Notes:
  #
  # Will plot the GAM model of dti data

  # plot
  df_pred <- predict.bam(
    model,
    exclude_terms = c("pds", "sex", "subjectID"),
    values = list(pds = NULL, sex = NULL),
    se.fit = T,
    type = "response"
  )

  df_pred <- data.frame(
    group = df_tract$group,
    sex = df_tract$sex,
    subjectID = df_tract$subjectID,
    pds = df_tract$pds,
    nodeID = df_tract$nodeID,
    fit = df_pred$fit,
    se.fit = df_pred$se.fit
  )

  h_tract <- switch_tract_name(tract)
  h_title <- paste0("GAM Fit of ", h_tract, " FA Values")

  h_cols <- c(
    switch_plot_values("0", g_type)[[1]][1],
    switch_plot_values("1", g_type)[[1]][1],
    switch_plot_values("2", g_type)[[1]][1]
  )
  names(h_cols) <- c("0", "1", "2")
  h_breaks <- c("0", "1", "2")
  h_labels <- c(
    switch_plot_values("0", g_type)[[2]][1],
    switch_plot_values("1", g_type)[[2]][1],
    switch_plot_values("2", g_type)[[2]][1]
  )

  p <- ggplot(data = df_pred) +
    geom_smooth(mapping = aes(x = nodeID, y = fit, color = group)) +
    ggtitle(h_title) +
    ylab("Fit FA") +
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
    paste0(gam_plot_dir, "Plot_GAM_", tract, "_", "G", g_type, ".png"),
    units = "in",
    width = 6,
    height = 6,
    device = "png"
  )
}

calc_gam_stats <- function(tract, df_tract, g_type) {

  ### --- Notes:
  #
  # This function has 2 main steps
  #
  # 1) Determine best GAM family.
  #   This is based off the distribution of the
  #   data. It was determined that a gamma
  #   fit better than beta for all tracts.
  #   It was also determined that k=40 worked
  #   for all tracts.
  #
  # 2) Model with covariates, compare models.
  #   All tracts had better model fit when pds
  #   was added. The covariate model is then
  #   plotted.
  #
  # Note: this work can only be a function since
  #   it happened that all tracts had the same
  #   distribution, family, K, and best model.

  ### Model tract, no covariates
  # # plot mean data
  # ggplot(data = df_tract) +
  #   geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=group))
  #
  # ggplot(data = df_tract) +
  #   geom_point(mapping = aes(x=nodeID, y=dti_fa, color=group),size=0.3) +
  #   geom_smooth(mapping = aes(x=nodeID, y=dti_fa, color=group))

  # # determine distribution
  # descdist(df_tract$dti_fa, discrete=F) # Could be beta or gamma
  #
  # fit.beta <- fitdist(df_tract$dti_fa, "beta")
  # plot(fit.beta)
  # fit.beta$aic
  #
  # fit.gamma <- fitdist(df_tract$dti_fa, "gamma")
  # plot(fit.gamma)
  # fit.gamma$aic

  #  determine k, compare families
  fit_gamma <- bam(dti_fa ~ group +
    sex +
    s(nodeID, by = group, k = 40) +
    s(subjectID, bs = "re"),
  data = df_tract,
  family = Gamma(link = "logit"),
  method = "REML"
  )

  # gam.check(fit_gamma, rep = 500)
  capture.output(
    summary(fit_gamma),
    file = paste0(
      gam_stats_dir, "Stats_GAM-gamma_", tract, "_", "G", g_type, ".txt"
    )
  )

  fit_beta <- bam(dti_fa ~ group +
    sex +
    s(nodeID, by = group, k = 40) +
    s(subjectID, bs = "re"),
  data = df_tract,
  family = betar(link = "logit"),
  method = "REML"
  )

  # gam.check(fit_beta, rep = 500)
  capture.output(
    summary(fit_beta),
    file = paste0(
      gam_stats_dir, "Stats_GAM-beta_", tract, "_", "G", g_type, ".txt"
    )
  )

  if (tract == "FA") {
    fit_gaussian <- bam(dti_fa ~ group +
      sex +
      s(nodeID, by = group, k = 40) +
      s(subjectID, bs = "re"),
    data = df_tract,
    family = gaussian(link = "logit"),
    method = "REML"
    )

    # gam.check(fit_beta, rep = 500)
    capture.output(
      summary(fit_gaussian),
      file = paste0(
        gam_stats_dir,
        "Stats_GAM-gaussian_",
        tract, "_", "G", g_type, ".txt"
      )
    )
  }


  # infoMessages('on')
  # compareML(fit_gamma, fit_beta)  # fit_gamma recommended
  capture.output(
    compareML(fit_gamma, fit_beta),
    file = paste0(
      gam_stats_dir,
      "Stats_GAM-comp_gam-beta_",
      tract, "_", "G", g_type, ".txt"
    )
  )

  if (tract == "FA") {
    capture.output(
      compareML(fit_gamma, fit_gaussian),
      file = paste0(
        gam_stats_dir,
        "Stats_GAM-comp_gam-gaus_",
        tract, "_", "G", g_type, ".txt"
      )
    )
  }


  ### Model tract with covariates
  if (tract == "FA") {
    fit_cov_pds <- bam(dti_fa ~ group +
      sex +
      s(nodeID, by = group, k = 40) +
      s(pds, by = sex) +
      s(subjectID, bs = "re"),
    data = df_tract,
    family = gaussian(link = "logit"),
    method = "REML"
    )
  } else {
    fit_cov_pds <- bam(dti_fa ~ group +
      sex +
      s(nodeID, by = group, k = 40) +
      s(pds, by = sex) +
      s(subjectID, bs = "re"),
    data = df_tract,
    family = Gamma(link = "logit"),
    method = "REML"
    )
  }


  # gam.check(fit_cov_pds, rep = 500)
  capture.output(
    summary(fit_cov_pds),
    file = paste0(
      gam_stats_dir,
      "Stats_GAM-cov_",
      tract, "_", "G", g_type, ".txt"
    )
  )

  # Test cov model against gamma
  if (tract == "FA") {
    capture.output(
      compareML(fit_gaussian, fit_cov_pds),
      file = paste0(
        gam_stats_dir,
        "Stats_GAM-comp_gaus-cov_",
        tract, "_", "G", g_type, ".txt"
      )
    )
  } else {
    capture.output(
      compareML(fit_gamma, fit_cov_pds),
      file = paste0(
        gam_stats_dir,
        "Stats_GAM-comp_gam-cov_",
        tract, "_", "G", g_type, ".txt"
      )
    )
  }
  return(fit_cov_pds)
}

plot_spline_diff_pair <- function(model, tract, g_type, factor_a, factor_b) {

  ### --- Notes:
  #
  # This will make plots and write tables of sig
  #   node differences for GAM splines bx 2 factors (groups)

  png(
    filename = paste0(
      gam_plot_dir, "Plot_Diff_", tract, "_", "G", g_type, "_pair.png"
    ),
    width = 600, height = 600
  )

  group_a <- switch_plot_values(factor_a, g_type)[[2]][1]
  group_b <- switch_plot_values(factor_b, g_type)[[2]][1]

  par(mar = c(5, 5, 4, 2), family = "Times New Roman")
  capture.output(plot_diff(model,
    view = "nodeID",
    comp = list(group = c(factor_a, factor_b)),
    rm.ranef = T,
    main = paste0(
      "Difference Scores, ", group_a, "-", group_b
    ),
    ylab = "Est. FA difference",
    xlab = "Tract Node",
    cex.lab = 2,
    cex.axis = 2,
    cex.main = 2,
    cex.sub = 1.5
  ),
  file = paste0(
    table_dir,
    "Table_Diff_",
    tract, "_", "G",
    g_type, "_",
    factor_a, factor_b,
    ".txt"
  )
  )
  par(mar = c(5, 4, 4, 2))
  dev.off()
}

plot_split_diff_mult <- function(model, tract, g_type) {

  ### --- Notes:
  #
  # This will draw plots and write tables of sig
  #   node differences for GAM when group style=2

  tiff(
    filename = paste0(
      gam_plot_dir, "Plot_Diff_",
      tract, "_", "G", g_type, "_anova.tiff"
    ),
    width = 2000,
    height = 600
  )

  group_a <- switch_plot_values("0", g_type)[[2]][1]
  group_b <- switch_plot_values("1", g_type)[[2]][1]
  group_c <- switch_plot_values("2", g_type)[[2]][1]

  par(mfrow = c(1, 3), mar = c(5, 6, 4, 2), family = "Times New Roman")
  capture.output(plot_diff(model,
    view = "nodeID",
    comp = list(group = c("0", "1")),
    rm.ranef = T,
    main = paste0(
      "Difference Scores, ", group_a, "-", group_b
    ),
    ylab = "Est. FA difference",
    xlab = "Tract Node",
    cex.lab = 3,
    cex.axis = 3,
    cex.main = 3.5,
    cex.sub = 2.5
  ),
  file = paste0(
    table_dir,
    "Table_Diff_",
    tract, "_", "G",
    g_type, "_01.txt"
  )
  )

  par(mar = c(5, 3, 4, 2))
  capture.output(plot_diff(model,
    view = "nodeID",
    comp = list(group = c("0", "2")),
    rm.ranef = T,
    main = paste0(
      "Difference Scores, ", group_a, "-", group_c
    ),
    ylab = "",
    xlab = "Tract Node",
    cex.lab = 3,
    cex.axis = 3,
    cex.main = 3.5,
    cex.sub = 2.5
  ),
  file = paste0(
    table_dir,
    "Table_Diff_",
    tract, "_", "G",
    g_type, "_02.txt"
  )
  )

  capture.output(plot_diff(model,
    view = "nodeID",
    comp = list(group = c("1", "2")),
    rm.ranef = T,
    main = paste0(
      "Difference Scores, ", group_b, "-", group_c
    ),
    ylab = "",
    xlab = "Tract Node",
    cex.lab = 3,
    cex.axis = 3,
    cex.main = 3.5,
    cex.sub = 2.5
  ),
  file = paste0(
    table_dir,
    "Table_Diff_",
    tract, "_", "G",
    g_type, "_12.txt"
  )
  )
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
  dev.off()
}

make_spline_diff_pair_df <- function(model, factor_a, factor_b) {

  ### --- Notes:
  #
  # Returns a dataframe of difference
  #   scores for each node.

  df_pair <- plot_diff(model,
    view = "nodeID",
    comp = list(group = c(factor_a, factor_b)),
    rm.ranef = T,
    plot = F
  )

  colnames(df_pair) <- c(colnames(df_pair[, 1:4]), "Comp")
  df_pair$Comp <- paste0(factor_a, factor_b)
  return(df_pair)
}

make_spline_diff_mult_df <- function(model) {

  ### --- Notes:
  #
  # Returns a dataframe of difference
  #   scores for each node.

  p01 <- plot_diff(model,
    view = "nodeID",
    comp = list(group = c("0", "1")),
    rm.ranef = T,
    plot = F
  )
  colnames(p01) <- c(colnames(p01[, 1:4]), "Comp")
  p01$Comp <- "01"

  p02 <- plot_diff(model,
    view = "nodeID",
    comp = list(group = c("0", "2")),
    rm.ranef = T,
    plot = F
  )
  colnames(p02) <- c(colnames(p02[, 1:4]), "Comp")
  p02$Comp <- "02"

  p12 <- plot_diff(model,
    view = "nodeID",
    comp = list(group = c("1", "2")),
    rm.ranef = T,
    plot = F
  )
  colnames(p12) <- c(colnames(p12[, 1:4]), "Comp")
  p12$Comp <- "12"

  df_out <- as.data.frame(
    matrix(NA, nrow = 3 * dim(p01)[1], ncol = dim(p01)[2])
  )
  colnames(df_out) <- colnames(p01)
  df_out <- rbind(p01, p02, p12)
  return(df_out)
}

calc_spline_diff <- function(model, tract, g_type, pair_anov, comp_list) {

  ### --- Notes:
  #
  # 1) Makes plots and tables
  #     Tables are not currently used,
  #     func_plot_diff and func_mkdf_diff could
  #     be combined.
  #
  # 2) Make dataframes of difference estimates
  #
  # 3) Determine nodes that differ in difference
  #     estimation. Repeats what is done in tables.
  #     Also gets max node.
  #
  # Returns list of nodes, max node

  # make plots and tables
  if (pair_anov == "pair") {
    plot_spline_diff_pair(model, tract, g_type, comp_list[1], comp_list[2])
  } else if (pair_anov == "anova") {
    plot_split_diff_mult(model, tract, g_type)
  }

  # get plot_diff data frames
  if (pair_anov == "pair") {
    df_est_diff <- make_spline_diff_pair_df(model, comp_list[1], comp_list[2])
  } else if (pair_anov == "anova") {
    df_est_diff <- make_spline_diff_mult_df(model)
  }

  # determine where nodes differ, "anova" looks
  #   for nodes that differ bx ALL groups
  node_list <- unique(df_est_diff$nodeID)
  diff_list <- vector()

  for (node in node_list) {
    ind_node <- which(df_est_diff$nodeID == node)

    if (pair_anov == "pair") {
      h_est <- abs(df_est_diff[ind_node[1], ]$est)
      h_ci <- df_est_diff[ind_node[1], ]$CI

      if ((h_est - h_ci) > 0) {
        diff_list <- c(diff_list, node)
      }
    } else if (pair_anov == "anova") {
      if (
        (abs(df_est_diff[ind_node[1], ]$est) - df_est_diff[ind_node[1], ]$CI) > 0 &
          (abs(df_est_diff[ind_node[2], ]$est) - df_est_diff[ind_node[2], ]$CI) > 0 &
          (abs(df_est_diff[ind_node[3], ]$est) - df_est_diff[ind_node[3], ]$CI) > 0
      ) {
        diff_list <- c(diff_list, node)
      }
    }
  }

  if (length(diff_list) == 0) {
    return(NA)
    stop
  }

  # find node of max difference
  if (pair_anov == "pair") {
    h_df <- subset(df_est_diff, nodeID %in% diff_list)
    ind_max <- which(abs(h_df$est) == max(abs(h_df$est)))
    node_max <- h_df[ind_max, ]$nodeID
  } else if (pair_anov == "anova") {
    f_sum <- function(x) {
      sum(abs(df_est_diff[which(df_est_diff$nodeID == x), ]$est))
    }

    h_sum <- sapply(diff_list, f_sum)
    names(h_sum) <- diff_list
    h_max <- which(h_sum == max(h_sum))
    node_max <- as.numeric(names(h_sum[h_max]))
  }

  return(list(diff_list, node_max))
}


# LM Functions ------
make_diff_df <- function(df_tract, node_list, g_type, avg_max) {

  ### --- Notes:
  #
  # Returns a data frame for linear models
  #   with either average FA or max FA
  #   difference

  subj_list <- unique(df_tract$subjectID)
  df_out <- as.data.frame(matrix(NA, nrow = length(subj_list), ncol = 8))
  colnames(df_out) <- c(
    "subj", "fa_value", "pars6", "group",
    "neu_lgi", "neu_ldi", "neg_lgi", "neg_ldi"
  )
  df_out$subj <- subj_list

  for (subj in subj_list) {
    ind_out <- which(df_out$subj == subj)
    df_subj <- df_tract[which(df_tract$subjectID == subj), ]

    df_out[ind_out, ]$pars6 <- df_subj[1, ]$pars6
    df_out[ind_out, ]$group <- as.numeric(df_subj[1, ]$group) - 1
    df_out[ind_out, ]$neu_lgi <- df_subj[1, ]$neu_lgi
    df_out[ind_out, ]$neu_ldi <- df_subj[1, ]$neu_ldi
    df_out[ind_out, ]$neg_lgi <- df_subj[1, ]$neg_lgi
    df_out[ind_out, ]$neg_ldi <- df_subj[1, ]$neg_ldi

    if (avg_max == "avg") {
      df_node <- subset(df_subj, nodeID %in% node_list)
      df_out[ind_out, ]$fa_value <- round(mean(df_node$dti_fa), 4)
    } else if (avg_max == "max") {
      ind_max <- which(df_subj$nodeID == node_list)
      df_out[ind_out, ]$fa_value <- round(df_subj[ind_max, ]$dti_fa, 4)
    }
  }
  return(df_out)
}

plot_lm_pair <- function(df_plot, avg_max, mem, factor_a, factor_b) {

  ### --- Notes:
  #
  # Plots LM of A and B

  # determine labels
  h_labels <- c(
    switch_plot_values(factor_a, g_type)[[2]][1],
    switch_plot_values(factor_b, g_type)[[2]][1]
  )
  h_tract <- switch_tract_name(tract)

  # set up for plot A
  plot1.x <- df_plot[which(df_plot$group == as.numeric(factor_a)), ]$fa_value
  plot1.y <- df_plot[which(df_plot$group == as.numeric(factor_a)), ]$mem_score

  # set up for plot B
  plot2.x <- df_plot[which(df_plot$group == as.numeric(factor_b)), ]$fa_value
  plot2.y <- df_plot[which(df_plot$group == as.numeric(factor_b)), ]$mem_score

  # omni title
  h_title <- paste(
    h_tract, "Spline Differences Predicting Memory Performance"
  )
  x_title <- ifelse(avg_max == "avg", "Mean FA", "Max FA")

  # set up, draw
  tiff(
    filename = paste0(
      lm_plot_dir,
      "Plot_LM-",
      avg_max, "_",
      tract, "_G",
      g_type, "_",
      mem, "_pair.tiff"
    ),
    width = 8,
    height = 4,
    units = "in",
    res = 600
  )

  par(
    mfrow = c(1, 2),
    oma = c(0, 0, 2, 0),
    family = "Times New Roman"
  )

  plot(plot1.x, plot1.y,
    xlab = x_title,
    ylab = mem,
    ylim = c(min(df_plot$mem_score), max(df_plot$mem_score)),
    main = h_labels[1]
  )
  abline(lm(plot1.y ~ plot1.x))

  plot(plot2.x, plot2.y,
    xlab = x_title,
    ylab = "",
    main = h_labels[2],
    ylim = c(min(df_plot$mem_score), max(df_plot$mem_score))
  )
  abline(lm(plot2.y ~ plot2.x))

  mtext(h_title, outer = T, cex = 1.5)
  dev.off()
  par(mfrow = c(1, 1))
}

plot_lm_mult <- function(df_plot, avg_max, mem, g_type) {

  ### --- Notes:
  #
  # Plots A, B, and C

  h_labels <- c(
    switch_plot_values("0", g_type)[[2]][1],
    switch_plot_values("1", g_type)[[2]][1],
    switch_plot_values("2", g_type)[[2]][1]
  )
  h_tract <- switch_tract_name(tract)

  plot1.x <- df_plot[which(df_plot$group == 0), ]$fa_value
  plot1.y <- df_plot[which(df_plot$group == 0), ]$mem_score

  plot2.x <- df_plot[which(df_plot$group == 1), ]$fa_value
  plot2.y <- df_plot[which(df_plot$group == 1), ]$mem_score

  plot3.x <- df_plot[which(df_plot$group == 2), ]$fa_value
  plot3.y <- df_plot[which(df_plot$group == 2), ]$mem_score

  h_title <- paste(
    h_tract, "Spline Differences Predicting Memory Performance"
  )
  x_title <- ifelse(avg_max == "avg", "Mean FA", "Max FA")

  tiff(
    filename = paste0(
      lm_plot_dir, "Plot_LM-",
      avg_max, "_",
      tract, "_G",
      g_type, "_",
      mem, "_anova.tiff"
    ),
    width = 8,
    height = 4,
    units = "in",
    res = 600
  )

  par(
    mfrow = c(1, 3),
    oma = c(0, 0, 2, 0),
    family = "Times New Roman"
  )

  plot(plot1.x, plot1.y,
    xlab = x_title,
    cex.lab = 1.5,
    ylab = mem,
    ylim = c(min(df_plot$mem_score), max(df_plot$mem_score)),
    main = h_labels[1]
  )
  abline(lm(plot1.y ~ plot1.x))

  plot(plot2.x, plot2.y,
    xlab = x_title,
    cex.lab = 1.5,
    ylab = "",
    main = h_labels[2],
    ylim = c(min(df_plot$mem_score), max(df_plot$mem_score))
  )
  abline(lm(plot2.y ~ plot2.x))

  plot(plot3.x, plot3.y,
    xlab = x_title,
    cex.lab = 1.5,
    ylab = "",
    main = h_labels[3],
    ylim = c(min(df_plot$mem_score), max(df_plot$mem_score))
  )
  abline(lm(plot3.y ~ plot3.x))

  mtext(h_title, outer = T, cex = 1.5)
  dev.off()
  par(mfrow = c(1, 1))
}

calc_lm_stats <- function(
                          df_lm, tract, g_type, avg_max, pair_anov, comp_list) {

  ### --- Notes:
  #
  # Conduct linear model for list of mem scores
  #   then make plots.

  # mem_list <- c("neu_lgi", "neu_ldi","neg_lgi", "neg_ldi")
  mem_list <- "neg_lgi"

  for (mem in mem_list) {

    # make dataframe for memory behavior
    df_mem <- as.data.frame(matrix(NA, nrow = dim(df_lm)[1], ncol = 4))
    colnames(df_mem) <- c("subj", "fa_value", "group", "mem_score")

    df_mem$subj <- df_lm$subj
    df_mem$fa_value <- df_lm$fa_value
    df_mem$group <- df_lm$group

    # get/add appropriate behavior score
    ind_mem <- grep(mem, colnames(df_lm))
    df_mem$mem_score <- df_lm[, ind_mem]

    # subset for pairwise comparison
    if (pair_anov == "pair") {
      df_mem <- df_mem[which(
        df_mem$group == comp_list[1] | df_mem$group == comp_list[2]
      ), ]
    }

    # linear regression
    fit.int <- lm(mem_score ~ fa_value * group, data = df_mem)
    capture.output(
      summary(fit.int),
      file = paste0(
        lm_stats_dir,
        "Stats_LM-", avg_max, "_",
        tract, "_G", g_type, "_",
        mem, "_", pair_anov, ".txt"
      )
    )
    etaSquared(fit.int)

    # make plots
    if (pair_anov == "pair") {
      plot_lm_pair(df_mem, avg_max, mem, comp_list[1], comp_list[2])
    } else if (pair_anov == "anova") {
      plot_lm_mult(df_mem, avg_max, mem, g_type)
    }
  }
}



# Work ------
#
# Two analyses (grouping types):
#   1) Con vs Anx (removed)
#   2) Con vs GAD vs SAD
#   3) PARS Low vs Med vs High

# group_type <- c(1, 2, 3)
group_type <- 3
tract_list <- c("UNC_L", "UNC_R", "FA")

# for spline comparisons, linear models,
#   determine whether to compare all (anova)
#   or a pair. If pair, set factors as string.

pair_anov <- "pair"
# g2_pair_list <- c("0", "2")
g3_pair_list <- c("0", "2")


for (g_type in group_type) {

  # make/get data, assign factor
  data_file <- paste0(data_dir, "Master_dataframe_G", g_type, ".csv")

  if (!file.exists(data_file)) {
    make_master_df(g_type)
  }

  df_afq <- read.csv(data_file)
  df_afq$group <- factor(df_afq$group)
  df_afq$sex <- factor(df_afq$sex)

  # Check Memory behavior
  calc_memory_stats(df_afq, g_type)

  # get comp list
  if (pair_anov == "pair") {
    h_var <- paste0("g", g_type, "_pair_list")
    comp_list <- get(h_var)
  }

  for (tract in tract_list) {

    # subset df_afq for tract
    df_tract <- df_afq[which(df_afq$tractID == tract), ]
    df_tract$dti_fa <- round(df_tract$dti_fa, 3)

    # run gam, plot
    gamFile <- paste0(data_dir, "G", g_type, "_gam_", tract, ".Rda")

    if (!file.exists(gamFile)) {
      h_gam <- calc_gam_stats(tract, df_tract, g_type)
      saveRDS(h_gam, file = gamFile)
      rm(h_gam)
    }

    gam_model <- readRDS(gamFile)
    help_plot_gam(gam_model, tract, g_type, df_tract)

    # determine nodes of group differences
    #   Use "anova" to compare bx >2 groups, and note
    #     that differences between splines will be for
    #     nodes which differed between ALL groups
    node_list <- calc_spline_diff(
      gam_model, tract, g_type, pair_anov, comp_list
    )

    # deal w/no differences
    if (is.na(node_list)) {
      next
    }

    # avg lm
    avg_nList <- node_list[[1]]
    df_avg <- make_diff_df(df_tract, avg_nList, g_type, "avg")
    calc_lm_stats(df_avg, tract, g_type, "avg", pair_anov, comp_list)

    # # max lm
    # max_nList <- node_list[[2]]
    # df_max <- make_diff_df(df_tract, max_nList, g_type, "max")
    # calc_lm_stats(df_max, tract, g_type, "max", pair_anov, comp_list)
  }
}
