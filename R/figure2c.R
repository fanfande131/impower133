#' Plot Figure 2C: Subgroup forest plot
#'
#' Calculate hazard ratios for each subgroup based on simulated
#' baseline data merged with reconstructed OS data, and generate
#' a forest plot in the original IMpower133 publication style.
#'
#' @param data Data frame from simulate_impower133() with os_path specified
#' @param output_path Path to output image file, default "figure2c.png"
#' @return No return value, saves image to output_path
#' @export
plot_figure2c <- function(data, output_path = "figure2c.png") {
  requireNamespace("survival", quietly = TRUE)
  requireNamespace("survminer", quietly = TRUE)
  requireNamespace("forestplot", quietly = TRUE)

  # Placebo as reference
  data$treatment <- factor(data$treatment, levels = c("Placebo", "Atezolizumab"))

  # ========== Helper functions ==========
  fmt_med <- function(x) ifelse(is.na(x), "NR", sprintf("%.1f", x))

  calc_hr <- function(dat, var, level) {
    sub <- dat[dat[[var]] == level & !is.na(dat[[var]]), ]
    if (nrow(sub) < 10) return(c(NA, NA, NA, NA, NA, NA))

    fit <- survival::survfit(
      survival::Surv(os_time, os_status) ~ treatment, data = sub)
    med_tab <- survminer::surv_median(fit)

    med_placebo <- med_tab$median[med_tab$strata == "treatment=Placebo"]
    med_atezo <- med_tab$median[med_tab$strata == "treatment=Atezolizumab"]

    fit_cox <- survival::coxph(
      survival::Surv(os_time, os_status) ~ treatment, data = sub)
    hr <- exp(stats::coef(fit_cox))[1]
    ci <- exp(stats::confint(fit_cox))[1, ]

    c(nrow(sub), med_placebo, med_atezo, hr, ci[1], ci[2])
  }

  # ========== Calculate subgroups ==========
  subgroup_list <- list(
    c("sex", "Male"), c("sex", "Female"),
    c("age_group", "<65 yr"), c("age_group", ">=65 yr"),
    c("ecog", "0"), c("ecog", "1"),
    c("brain_mets", "Yes"), c("brain_mets", "No"),
    c("tmb_10", "<10 mut/Mb"), c("tmb_10", ">=10 mut/Mb"),
    c("tmb_16", "<16 mut/Mb"), c("tmb_16", ">=16 mut/Mb")
  )

  results <- do.call(rbind,lapply(subgroup_list,function(s) calc_hr(data,s[1],s[2])))

  # Overall population
  fit_all <- survival::survfit(
    survival::Surv(os_time, os_status) ~ treatment, data = data)
  med_all <- survminer::surv_median(fit_all)
  med_placebo_all <- med_all$median[med_all$strata == "treatment=Placebo"]
  med_atezo_all <- med_all$median[med_all$strata == "treatment=Atezolizumab"]

  cox_all <- survival::coxph(
    survival::Surv(os_time, os_status) ~ treatment, data = data)
  total <- c(
    nrow(data),
    med_placebo_all,
    med_atezo_all,
    exp(stats::coef(cox_all))[1],
    exp(stats::confint(cox_all))[1, 1],
    exp(stats::confint(cox_all))[1, 2]
  )

  # ========== Build data frame ==========
  forest_df <- data.frame(
    subgroup = c(
      "Sex", "  Male", "  Female",
      "Age", "  <65 yr", "  >=65 yr",
      "ECOG score", "  0", "  1",
      "Brain metastases", "  Yes", "  No",
      "Tumor mutational burden",
      "  <10 mutations/Mb", "  >=10 mutations/Mb",
      "  <16 mutations/Mb", "  >=16 mutations/Mb",
      "Intention-to-treat population", "  All patients"
    ),
    n_pct = c(
      "",
      paste0(results[1, 1], " (", round(100 * results[1, 1] / nrow(data)), ")"),
      paste0(results[2, 1], " (", round(100 * results[2, 1] / nrow(data)), ")"),
      "",
      paste0(results[3, 1], " (", round(100 * results[3, 1] / nrow(data)), ")"),
      paste0(results[4, 1], " (", round(100 * results[4, 1] / nrow(data)), ")"),
      "",
      paste0(results[5, 1], " (", round(100 * results[5, 1] / nrow(data)), ")"),
      paste0(results[6, 1], " (", round(100 * results[6, 1] / nrow(data)), ")"),
      "",
      paste0(results[7, 1], " (", round(100 * results[7, 1] / nrow(data)), ")"),
      paste0(results[8, 1], " (", round(100 * results[8, 1] / nrow(data)), ")"),
      "",
      paste0(results[9, 1], " (", round(100 * results[9, 1] / nrow(data)), ")"),
      paste0(results[10, 1], " (", round(100 * results[10, 1] / nrow(data)), ")"),
      paste0(results[11, 1], " (", round(100 * results[11, 1] / nrow(data)), ")"),
      paste0(results[12, 1], " (", round(100 * results[12, 1] / nrow(data)), ")"),
      "",
      paste0(total[1], " (100)")
    ),
    med_atezo = c(
      "", fmt_med(results[1, 3]), fmt_med(results[2, 3]),
      "", fmt_med(results[3, 3]), fmt_med(results[4, 3]),
      "", fmt_med(results[5, 3]), fmt_med(results[6, 3]),
      "", fmt_med(results[7, 3]), fmt_med(results[8, 3]),
      "", fmt_med(results[9, 3]), fmt_med(results[10, 3]),
      fmt_med(results[11, 3]), fmt_med(results[12, 3]),
      "", fmt_med(total[3])
    ),
    med_placebo = c(
      "", fmt_med(results[1, 2]), fmt_med(results[2, 2]),
      "", fmt_med(results[3, 2]), fmt_med(results[4, 2]),
      "", fmt_med(results[5, 2]), fmt_med(results[6, 2]),
      "", fmt_med(results[7, 2]), fmt_med(results[8, 2]),
      "", fmt_med(results[9, 2]), fmt_med(results[10, 2]),
      fmt_med(results[11, 2]), fmt_med(results[12, 2]),
      "", fmt_med(total[2])
    ),
    hr = c(
      NA, results[1, 4], results[2, 4],
      NA, results[3, 4], results[4, 4],
      NA, results[5, 4], results[6, 4],
      NA, results[7, 4], results[8, 4],
      NA, results[9, 4], results[10, 4], results[11, 4], results[12, 4],
      NA, total[4]
    ),
    lower = c(
      NA, results[1, 5], results[2, 5],
      NA, results[3, 5], results[4, 5],
      NA, results[5, 5], results[6, 5],
      NA, results[7, 5], results[8, 5],
      NA, results[9, 5], results[10, 5], results[11, 5], results[12, 5],
      NA, total[5]
    ),
    upper = c(
      NA, results[1, 6], results[2, 6],
      NA, results[3, 6], results[4, 6],
      NA, results[5, 6], results[6, 6],
      NA, results[7, 6], results[8, 6],
      NA, results[9, 6], results[10, 6], results[11, 6], results[12, 6],
      NA, total[6]
    ),
    stringsAsFactors = FALSE
  )

  forest_df$hr_ci <- ifelse(
    is.na(forest_df$hr), "",
    sprintf("%.2f (%.2f-%.2f)", forest_df$hr, forest_df$lower, forest_df$upper)
  )

  # ========== Build table text ==========
  tabletext <- rbind(
    c("Subgroup", "No. of Patients (%)", "Median Overall Survival (mo)", "",
      "Hazard Ratio for Death\n(95% CI)"),
    c("", "", "Atezolizumab", "Placebo", ""),
    cbind(forest_df$subgroup, forest_df$n_pct,
          forest_df$med_atezo, forest_df$med_placebo, forest_df$hr_ci)
  )

  is_header <- is.na(forest_df$hr)

  grDevices::png(output_path, width = 2400, height = 1200, res = 200)
  # ========== Draw forest plot ==========
  p1<-forestplot::forestplot(
    labeltext = tabletext,
    mean      = c(NA, NA, forest_df$hr),
    lower     = c(NA, NA, forest_df$lower),
    upper     = c(NA, NA, forest_df$upper),
    is.summary= c(TRUE,TRUE,is_header),
    zero      = 1,
    clip      = c(0.1, 2.5),
    xlog      = TRUE,
    graph.pos = 5,
    boxsize   = 0.12,
    col       = forestplot::fpColors(box = "black", line = "gray30", summary = "black"),
    hrzl_lines = list(
      "3"  = grid::gpar(lwd = 1, col = "gray80"),
      "5"  = grid::gpar(lwd = 1, col = "gray85"),
      "8"  = grid::gpar(lwd = 1, col = "gray85"),
      "11" = grid::gpar(lwd = 1, col = "gray85"),
      "14" = grid::gpar(lwd = 1, col = "gray85"),
      "19" = grid::gpar(lty = 2, col = "gray60")
    ),
    txt_gp = forestplot::fpTxtGp(
      label  = grid::gpar(fontsize = 11),
      ticks  = grid::gpar(fontsize = 10),
      xlab   = grid::gpar(fontsize = 11),
      title  = grid::gpar(fontsize = 14, fontface = "bold")
    ),
    title = "C  Overall Survival According to Baseline Characteristics",
    xlab  = "Hazard Ratio"
  )
  print(p1)
  grDevices::dev.off()
  message("Figure 2C saved to ", output_path)
  invisible(NULL)
}
