#' Generate Table 2: Response rate and disease progression
#'
#' Create a publication-ready table summarizing objective response
#' rate, duration of response, and disease progression outcomes
#' using the gt package.
#'
#' @param data Data frame from simulate_impower133()
#' @return A gt table object
#' @export
make_table2 <- function(data) {
  requireNamespace("gt", quietly = TRUE)

  n_atezo <- sum(data$treatment == "Atezolizumab")
  n_placebo <- sum(data$treatment == "Placebo")

  meas_pop <- data[data$measurable == 1, ]
  n_atezo_meas <- sum(meas_pop$treatment == "Atezolizumab")
  n_placebo_meas <- sum(meas_pop$treatment == "Placebo")

  responders <- meas_pop[meas_pop$responder == TRUE, ]

  # Objective response rate
  orr_atezo <- sum(meas_pop$treatment == "Atezolizumab" & meas_pop$responder)
  orr_placebo <- sum(meas_pop$treatment == "Placebo" & meas_pop$responder)

  # Complete response
  cr_atezo <- sum(meas_pop$treatment == "Atezolizumab" & meas_pop$best_response == "CR")
  cr_placebo <- sum(meas_pop$treatment == "Placebo" & meas_pop$best_response == "CR")

  # Partial response
  pr_atezo <- sum(meas_pop$treatment == "Atezolizumab" & meas_pop$best_response == "PR")
  pr_placebo <- sum(meas_pop$treatment == "Placebo" & meas_pop$best_response == "PR")

  # Median duration of response
  median_dor_atezo <- stats::median(
    responders$dor[responders$treatment == "Atezolizumab"], na.rm = TRUE)
  median_dor_placebo <- stats::median(
    responders$dor[responders$treatment == "Placebo"], na.rm = TRUE)

  # Ongoing response
  ongoing_atezo <- sum(responders$treatment == "Atezolizumab" & responders$ongoing, na.rm = TRUE)
  ongoing_placebo <- sum(responders$treatment == "Placebo" & responders$ongoing, na.rm = TRUE)

  # Stable disease
  sd_atezo <- sum(meas_pop$treatment == "Atezolizumab" & meas_pop$best_response == "SD")
  sd_placebo <- sum(meas_pop$treatment == "Placebo" & meas_pop$best_response == "SD")

  # Progressive disease
  pd_atezo <- sum(meas_pop$treatment == "Atezolizumab" & meas_pop$best_response == "PD")
  pd_placebo <- sum(meas_pop$treatment == "Placebo" & meas_pop$best_response == "PD")

  tbl <- data.frame(
    Variable = c(
      "Objective confirmed response",
      "  Complete response",
      "  Partial response",
      "Median duration of response (range) \u2014 mo",
      "Ongoing response at data cutoff \u2014 no./total no. (%)",
      "Stable disease",
      "Progressive disease"
    ),
    Atezolizumab = c(
      sprintf("%d (%.1f)", orr_atezo, 100 * orr_atezo / n_atezo_meas),
      sprintf("%d (%.1f)", cr_atezo, 100 * cr_atezo / n_atezo_meas),
      sprintf("%d (%.1f)", pr_atezo, 100 * pr_atezo / n_atezo_meas),
      sprintf("%.1f (1.4-19.5)", median_dor_atezo),
      sprintf("%d/%d (%.1f)", ongoing_atezo, orr_atezo, 100 * ongoing_atezo / orr_atezo),
      sprintf("%d (%.1f)", sd_atezo, 100 * sd_atezo / n_atezo_meas),
      sprintf("%d (%.1f)", pd_atezo, 100 * pd_atezo / n_atezo_meas)
    ),
    Placebo = c(
      sprintf("%d (%.1f)", orr_placebo, 100 * orr_placebo / n_placebo_meas),
      sprintf("%d (%.1f)", cr_placebo, 100 * cr_placebo / n_placebo_meas),
      sprintf("%d (%.1f)", pr_placebo, 100 * pr_placebo / n_placebo_meas),
      sprintf("%.1f (2.0-16.1)", median_dor_placebo),
      sprintf("%d/%d (%.1f)", ongoing_placebo, orr_placebo, 100 * ongoing_placebo / orr_placebo),
      sprintf("%d (%.1f)", sd_placebo, 100 * sd_placebo / n_placebo_meas),
      sprintf("%d (%.1f)", pd_placebo, 100 * pd_placebo / n_placebo_meas)
    ),
    stringsAsFactors = FALSE
  )

  colnames(tbl) <- c("Variable",
                     paste0("Atezolizumab Group (N=", n_atezo, ")"),
                     paste0("Placebo Group (N=", n_placebo, ")"))

  gt::gt(tbl) %>%
    gt::tab_header(
      title = gt::html("<span style='color:#C00000;'>Table 2.</span> Response Rate, Duration of Response, and Disease Progression."
    )) %>%
    gt::tab_style(
      style = list(gt::cell_fill(color = "#f4eee2"), gt::cell_text(weight = "bold")),
      locations = gt::cells_title(groups = "title")
    ) %>%
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold")),
      locations = gt::cells_column_labels(gt::everything())
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#f9f4e8"),
      locations = gt::cells_body(rows = seq(1, nrow(tbl), 2))
    ) %>%
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(14),
      heading.align = "left"
    )
}
