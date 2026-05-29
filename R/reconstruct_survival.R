#' Reconstruct individual patient data from K-M curves
#'
#' Use IPDfromKM to extract coordinates from published Kaplan-Meier
#' curves and reconstruct individual patient data (IPD).
#'
#' @param img_path Path to K-M curve image
#' @param trisk Vector of risk time points
#' @param nrisk_control Number at risk for control group
#' @param nrisk_treatment Number at risk for treatment group
#' @param x1,x2 X-axis range (actual values)
#' @param y1,y2 Y-axis range (actual values)
#' @param output_csv Path to output CSV file
#' @param type "OS" or "PFS"
#' @return A data frame with columns time, status, arm (invisibly)
#' @export
reconstruct_survival <- function(img_path,
                                 trisk,
                                 nrisk_control,
                                 nrisk_treatment,
                                 x1 = 0, x2 = 21, y1 = 0, y2 = 1,
                                 output_csv = tempfile(fileext = ".csv"),
                                 type = c("OS", "PFS")) {

  type <- match.arg(type)

  requireNamespace("IPDfromKM", quietly = TRUE)

  message("Extracting ", type, " coordinates from image...")
  points_control   <- IPDfromKM::getpoints(img_path, x1 = x1, x2 = x2, y1 = y1, y2 = y2)
  points_treatment <- IPDfromKM::getpoints(img_path, x1 = x1, x2 = x2, y1 = y1, y2 = y2)

  message("Preprocessing coordinates...")
  pre_control   <- IPDfromKM::preprocess(dat = points_control,
                                         trisk = trisk,
                                         nrisk = nrisk_control,
                                         maxy = 1)
  pre_treatment <- IPDfromKM::preprocess(dat = points_treatment,
                                         trisk = trisk,
                                         nrisk = nrisk_treatment,
                                         maxy = 1)

  message("Reconstructing IPD...")
  ipd_control   <- IPDfromKM::getIPD(prep = pre_control, armID = 0)
  ipd_treatment <- IPDfromKM::getIPD(prep = pre_treatment, armID = 1)

  ipd_all <- rbind(ipd_control$IPD, ipd_treatment$IPD)
  names(ipd_all) <- c("Survival_time", "Status", "Treatment_group")
  ipd_all$Treatment_group <- factor(ipd_all$Treatment_group,
                                    levels = c(0, 1),
                                    labels = c("Placebo", "Atezolizumab"))

  if (is.null(output_csv)) {
    output_csv <- paste0("reconstructed_", tolower(type), ".csv")
  }
  utils::write.csv(ipd_all, output_csv, row.names = FALSE, fileEncoding = "UTF-8")
  message("Reconstructed ", type, " IPD saved to ", output_csv)

  return(invisible(ipd_all))
}


#' Analyze survival data and print key statistics
#'
#' @param ipd_data Reconstructed data frame or CSV path
#' @param type "OS" or "PFS"
#' @param times Landmark time point in months, default 12
#' @return A list of key statistics (invisibly)
#' @export
analyze_survival <- function(ipd_data,
                             type = c("OS", "PFS"),
                             times = 12) {

  type <- match.arg(type)
  requireNamespace("survival", quietly = TRUE)

  if (is.character(ipd_data)) {
    ipd_data <- utils::read.csv(ipd_data, header = TRUE)
  }

  ipd_data$Treatment_group <- factor(ipd_data$Treatment_group,
                                     levels = c("Placebo", "Atezolizumab"))

  km_fit <- survival::survfit(
    survival::Surv(Survival_time, Status) ~ Treatment_group,
    data = ipd_data
  )

  message("\n========== ", type, " Median Survival ==========")
  print(km_fit)

  summary_t <- summary(km_fit, times = times)
  message(sprintf("\n========== %d-Month %s Rate ==========", times, type))
  message(sprintf("Placebo: %.1f%% (95%% CI, %.1f%%-%.1f%%)",
                  summary_t$surv[1] * 100, summary_t$lower[1] * 100, summary_t$upper[1] * 100))
  message(sprintf("Atezolizumab: %.1f%% (95%% CI, %.1f%%-%.1f%%)",
                  summary_t$surv[2] * 100, summary_t$lower[2] * 100, summary_t$upper[2] * 100))

  cox_fit <- survival::coxph(
    survival::Surv(Survival_time, Status) ~ Treatment_group,
    data = ipd_data
  )
  hr <- exp(stats::coef(cox_fit))
  ci <- exp(stats::confint(cox_fit))
  pval <- summary(cox_fit)$coefficients[5]

  message(sprintf("\n========== %s Hazard Ratio ==========", type))
  message(sprintf("HR = %.2f (95%% CI, %.2f-%.2f), P = %.3f", hr, ci[1], ci[2], pval))

  invisible(list(
    km_fit = km_fit,
    median = summary(km_fit)$table[, "median"],
    surv_at_t = summary_t,
    hr = hr,
    hr_ci = ci,
    hr_pval = pval
  ))
}


#' Plot survival K-M curve
#'
#' @param ipd_data Reconstructed data frame or CSV path
#' @param type "OS" or "PFS"
#' @param output_path Path to output image file
#' @param width,height Image dimensions in inches
#' @return No return value, called for side effects (saves a plot to file).
#' @export
plot_survival <- function(ipd_data,
                          type = c("OS", "PFS"),
                          output_path = tempfile(fileext = ".png"),
                          width = 10, height = 6) {

  type <- match.arg(type)
  requireNamespace("survival", quietly = TRUE)
  requireNamespace("survminer", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  if (is.character(ipd_data)) {
    ipd_data <- utils::read.csv(ipd_data, header = TRUE)
  }

  ipd_data$Treatment_group <- factor(ipd_data$Treatment_group,
                                     levels = c("Placebo", "Atezolizumab"))

  km_fit <- survival::survfit(
    survival::Surv(Survival_time, Status) ~ Treatment_group,
    data = ipd_data
  )

  cox_fit <- survival::coxph(
    survival::Surv(Survival_time, Status) ~ Treatment_group,
    data = ipd_data
  )
  hr <- exp(stats::coef(cox_fit))
  ci <- exp(stats::confint(cox_fit))
  pval <- summary(cox_fit)$coefficients[5]

  title_map <- c(OS = "A Overall Survival", PFS = "B Progression-Free Survival")

  p <- survminer::ggsurvplot(
    fit = km_fit,
    data = ipd_data,
    risk.table = TRUE,
    risk.table.title = "",
    risk.table.y.text.col = FALSE,
    risk.table.y.text = TRUE,
    risk.table.height = 0.25,
    title = title_map[type],
    palette = c("#6c8fc5", "#efa23f"),
    linetype = "solid",
    linewidth = 1.0,
    legend.labs = c("Atezolizumab", "Placebo"),
    xlab = "Months",
    ylab = "Patients Who Survived(%)",
    xlim = c(0, 22),
    ylim = c(0, 1),
    break.x.by = 1,
    break.y.by = 0.1,
    surv.median.line = "hv",
    pval = TRUE,
    pval.size = 4,
    ggtheme = ggplot2::theme_classic() +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        aspect.ratio = 0.3,
        axis.line = ggplot2::element_line(color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.text = ggplot2::element_text(color = "black")
      )
  )

  p$table <- p$table +
    ggplot2::scale_x_continuous(limits = c(0, 21), breaks = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10, hjust = 1, color = "black"),
      axis.text.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 5),
      aspect.ratio = 0.1,
      plot.margin = ggplot2::margin(t = 5, b = 0, l = 10, r = 10)
    )

  if (is.null(output_path)) {
    output_path <- paste0("figure2", ifelse(type == "OS", "a", "b"), ".png")
  }

  grDevices::png(output_path, width = width, height = height, units = "in", res = 300)
  print(p)
  grDevices::dev.off()

  message("Figure 2", ifelse(type == "OS", "a", "b"), " saved to ", output_path)
  print(p)
  invisible(p)
}
