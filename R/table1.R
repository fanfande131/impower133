#' Generate Table 1: Baseline characteristics
#'
#' Create a publication-ready baseline characteristics table
#' using the gt package, formatted to match the original
#' IMpower133 publication style.
#'
#' @param data Data frame from simulate_impower133()
#' @return A gt table object
#' @export
make_table1 <- function(data) {
  requireNamespace("gt", quietly = TRUE)

  n_atezo <- sum(data$treatment == "Atezolizumab")
  n_placebo <- sum(data$treatment == "Placebo")

  df_atezo <- data[data$treatment == "Atezolizumab", ]
  df_placebo <- data[data$treatment == "Placebo", ]

  # ========== Helper functions ==========
  count_percent <- function(vec) {
    vec <- vec[!is.na(vec)]
    n <- sum(vec)
    paste0(n, "/", length(vec), " (", sprintf("%.1f", 100 * mean(vec)), "%)")
  }

  median_range <- function(vec) {
    c(sprintf("%.0f", stats::median(vec, na.rm = TRUE)),
      paste0(min(vec, na.rm = TRUE), "-", max(vec, na.rm = TRUE)))
  }

  # ========== Compute each cell ==========

  # Age
  age_med_atezo <- stats::median(df_atezo$age, na.rm = TRUE)
  age_range_atezo <- paste0(min(df_atezo$age, na.rm = TRUE), "-",
                            max(df_atezo$age, na.rm = TRUE))
  age_med_placebo <- stats::median(df_placebo$age, na.rm = TRUE)
  age_range_placebo <- paste0(min(df_placebo$age, na.rm = TRUE), "-",
                              max(df_placebo$age, na.rm = TRUE))

  # Age group
  age_lt65_atezo <- sum(df_atezo$age_group == "<65 yr", na.rm = TRUE)
  age_ge65_atezo <- sum(df_atezo$age_group == ">=65 yr", na.rm = TRUE)
  age_lt65_placebo <- sum(df_placebo$age_group == "<65 yr", na.rm = TRUE)
  age_ge65_placebo <- sum(df_placebo$age_group == ">=65 yr", na.rm = TRUE)

  # Sex
  male_atezo <- sum(df_atezo$sex == "Male", na.rm = TRUE)
  male_placebo <- sum(df_placebo$sex == "Male", na.rm = TRUE)

  # ECOG
  ecog0_atezo <- sum(df_atezo$ecog == 0, na.rm = TRUE)
  ecog1_atezo <- sum(df_atezo$ecog == 1, na.rm = TRUE)
  ecog0_placebo <- sum(df_placebo$ecog == 0, na.rm = TRUE)
  ecog1_placebo <- sum(df_placebo$ecog == 1, na.rm = TRUE)

  # Smoking
  never_atezo <- sum(df_atezo$smoking == "Never smoked", na.rm = TRUE)
  current_atezo <- sum(df_atezo$smoking == "Current smoker", na.rm = TRUE)
  former_atezo <- sum(df_atezo$smoking == "Former smoker", na.rm = TRUE)
  never_placebo <- sum(df_placebo$smoking == "Never smoked", na.rm = TRUE)
  current_placebo <- sum(df_placebo$smoking == "Current smoker", na.rm = TRUE)
  former_placebo <- sum(df_placebo$smoking == "Former smoker", na.rm = TRUE)

  # Brain metastases
  brain_atezo <- sum(df_atezo$brain_mets == "Yes", na.rm = TRUE)
  brain_placebo <- sum(df_placebo$brain_mets == "Yes", na.rm = TRUE)

  # TMB
  tmb10_low_atezo <- sum(df_atezo$tmb_10 == "<10 mut/Mb", na.rm = TRUE)
  tmb10_high_atezo <- sum(df_atezo$tmb_10 == ">=10 mut/Mb", na.rm = TRUE)
  tmb16_low_atezo <- sum(df_atezo$tmb_16 == "<16 mut/Mb", na.rm = TRUE)
  tmb16_high_atezo <- sum(df_atezo$tmb_16 == ">=16 mut/Mb", na.rm = TRUE)

  tmb10_low_placebo <- sum(df_placebo$tmb_10 == "<10 mut/Mb", na.rm = TRUE)
  tmb10_high_placebo <- sum(df_placebo$tmb_10 == ">=10 mut/Mb", na.rm = TRUE)
  tmb16_low_placebo <- sum(df_placebo$tmb_16 == "<16 mut/Mb", na.rm = TRUE)
  tmb16_high_placebo <- sum(df_placebo$tmb_16 == ">=16 mut/Mb", na.rm = TRUE)

  tmb_n_atezo <- sum(!is.na(df_atezo$tmb_10))
  tmb_n_placebo <- sum(!is.na(df_placebo$tmb_10))

  # Target lesion
  lesion_med_atezo <- stats::median(df_atezo$target_lesion, na.rm = TRUE)
  lesion_range_atezo <- paste0(min(df_atezo$target_lesion, na.rm = TRUE), "-",
                               max(df_atezo$target_lesion, na.rm = TRUE))
  lesion_med_placebo <- stats::median(df_placebo$target_lesion, na.rm = TRUE)
  lesion_range_placebo <- paste0(min(df_placebo$target_lesion, na.rm = TRUE), "-",
                                 max(df_placebo$target_lesion, na.rm = TRUE))

  # Prior treatment
  chemo_atezo <- sum(df_atezo$prior_chemo == "Yes", na.rm = TRUE)
  rt_atezo <- sum(df_atezo$prior_rt == "Yes", na.rm = TRUE)
  surgery_atezo <- sum(df_atezo$prior_surgery == "Yes", na.rm = TRUE)
  chemo_placebo <- sum(df_placebo$prior_chemo == "Yes", na.rm = TRUE)
  rt_placebo <- sum(df_placebo$prior_rt == "Yes", na.rm = TRUE)
  surgery_placebo <- sum(df_placebo$prior_surgery == "Yes", na.rm = TRUE)

  # ========== Build table content ==========
  Characteristic <- c(
    "Age \u2014 yr",
    "  Median",
    "  Range",
    "Age group \u2014 no. (%)",
    "  <65 yr",
    "  >=65 yr",
    "Male sex \u2014 no. (%)",
    "ECOG performance-status score \u2014 no. (%)",
    "  0",
    "  1",
    "Smoking status \u2014 no. (%)",
    "  Never smoked",
    "  Current smoker",
    "  Former smoker",
    "Brain metastasis at enrollment \u2014 no. (%)",
    "Blood-based tumor mutational burden \u2014 no./total no. (%)",
    "  <10 mutations/Mb",
    "  >=10 mutations/Mb",
    "  <16 mutations/Mb",
    "  >=16 mutations/Mb",
    "Median sum of longest diameter of target lesions at baseline (range)",
    "Previous anticancer treatments \u2014 no. (%)",
    "  Chemotherapy or nonanthracycline",
    "  Radiotherapy",
    "  Cancer-related surgery"
  )

  Atezolizumab <- c(
    "", sprintf("%.0f", age_med_atezo), age_range_atezo,
    "", sprintf("%d (%.1f%%)", age_lt65_atezo, 100 * age_lt65_atezo / n_atezo),
    sprintf("%d (%.1f%%)", age_ge65_atezo, 100 * age_ge65_atezo / n_atezo),
    sprintf("%d (%.1f%%)", male_atezo, 100 * male_atezo / n_atezo),
    "", sprintf("%d (%.1f%%)", ecog0_atezo, 100 * ecog0_atezo / n_atezo),
    sprintf("%d (%.1f%%)", ecog1_atezo, 100 * ecog1_atezo / n_atezo),
    "",
    sprintf("%d (%.1f%%)", never_atezo, 100 * never_atezo / n_atezo),
    sprintf("%d (%.1f%%)", current_atezo, 100 * current_atezo / n_atezo),
    sprintf("%d (%.1f%%)", former_atezo, 100 * former_atezo / n_atezo),
    sprintf("%d (%.1f%%)", brain_atezo, 100 * brain_atezo / n_atezo),
    "",
    sprintf("%d/%d (%.1f%%)", tmb10_low_atezo, tmb_n_atezo, 100 * tmb10_low_atezo / tmb_n_atezo),
    sprintf("%d/%d (%.1f%%)", tmb10_high_atezo, tmb_n_atezo, 100 * tmb10_high_atezo / tmb_n_atezo),
    sprintf("%d/%d (%.1f%%)", tmb16_low_atezo, tmb_n_atezo, 100 * tmb16_low_atezo / tmb_n_atezo),
    sprintf("%d/%d (%.1f%%)", tmb16_high_atezo, tmb_n_atezo, 100 * tmb16_high_atezo / tmb_n_atezo),
    sprintf("%.1f (%s)", lesion_med_atezo, lesion_range_atezo),
    "",
    sprintf("%d (%.1f%%)", chemo_atezo, 100 * chemo_atezo / n_atezo),
    sprintf("%d (%.1f%%)", rt_atezo, 100 * rt_atezo / n_atezo),
    sprintf("%d (%.1f%%)", surgery_atezo, 100 * surgery_atezo / n_atezo)
  )

  Placebo <- c(
    "", sprintf("%.0f", age_med_placebo), age_range_placebo,
    "", sprintf("%d (%.1f%%)", age_lt65_placebo, 100 * age_lt65_placebo / n_placebo),
    sprintf("%d (%.1f%%)", age_ge65_placebo, 100 * age_ge65_placebo / n_placebo),
    sprintf("%d (%.1f%%)", male_placebo, 100 * male_placebo / n_placebo),
    "", sprintf("%d (%.1f%%)", ecog0_placebo, 100 * ecog0_placebo / n_placebo),
    sprintf("%d (%.1f%%)", ecog1_placebo, 100 * ecog1_placebo / n_placebo),
    "",
    sprintf("%d (%.1f%%)", never_placebo, 100 * never_placebo / n_placebo),
    sprintf("%d (%.1f%%)", current_placebo, 100 * current_placebo / n_placebo),
    sprintf("%d (%.1f%%)", former_placebo, 100 * former_placebo / n_placebo),
    sprintf("%d (%.1f%%)", brain_placebo, 100 * brain_placebo / n_placebo),
    "",
    sprintf("%d/%d (%.1f%%)", tmb10_low_placebo, tmb_n_placebo, 100 * tmb10_low_placebo / tmb_n_placebo),
    sprintf("%d/%d (%.1f%%)", tmb10_high_placebo, tmb_n_placebo, 100 * tmb10_high_placebo / tmb_n_placebo),
    sprintf("%d/%d (%.1f%%)", tmb16_low_placebo, tmb_n_placebo, 100 * tmb16_low_placebo / tmb_n_placebo),
    sprintf("%d/%d (%.1f%%)", tmb16_high_placebo, tmb_n_placebo, 100 * tmb16_high_placebo / tmb_n_placebo),
    sprintf("%.1f (%s)", lesion_med_placebo, lesion_range_placebo),
    "",
    sprintf("%d (%.1f%%)", chemo_placebo, 100 * chemo_placebo / n_placebo),
    sprintf("%d (%.1f%%)", rt_placebo, 100 * rt_placebo / n_placebo),
    sprintf("%d (%.1f%%)", surgery_placebo, 100 * surgery_placebo / n_placebo)
  )

  table1_df <- data.frame(
    Characteristic = Characteristic,
    Atezolizumab = Atezolizumab,
    Placebo = Placebo,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  colnames(table1_df) <- c("Characteristic",
                           paste0("Atezolizumab Group (N=", n_atezo, ")"),
                           paste0("Placebo Group (N=", n_placebo, ")"))

  # Section header rows
  section_rows <- Characteristic %in% c(
    "Age \u2014 yr", "Age group \u2014 no. (%)", "Male sex \u2014 no. (%)",
    "ECOG performance-status score \u2014 no. (%)", "Smoking status \u2014 no. (%)",
    "Brain metastasis at enrollment \u2014 no. (%)",
    "Blood-based tumor mutational burden \u2014 no./total no. (%)",
    "Median sum of longest diameter of target lesions at baseline (range)",
    "Previous anticancer treatments \u2014 no. (%)"
  )

  # ========== Build gt table ==========
  gt::gt(table1_df) %>%
    gt::tab_header(
      title = gt::html(
        "<span style='color:#C00000;'>Table 1.</span>
      Demographic and Disease Characteristics at Baseline in the Intention-to-Treat Population."
      )
    ) %>%
    gt::cols_align(align = "left", columns = "Characteristic") %>%
    gt::cols_align(align = "center", columns = 2:3) %>%
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
      locations = gt::cells_body(rows = seq(1, nrow(table1_df), 2))
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(rows = section_rows)
    ) %>%
    gt::tab_options(
      table.font.size = gt::px(13),
      heading.title.font.size = gt::px(14),
      heading.align = "left",
      table.border.top.color = "gray50",
      table.border.bottom.color = "gray50",
      column_labels.border.bottom.color = "gray50",
      data_row.padding = gt::px(4)
    )
}
