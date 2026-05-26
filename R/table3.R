#' Generate Table 3: Adverse events
#'
#' Create a publication-ready adverse events table summarizing
#' treatment-related adverse events by grade in the as-treated
#' population, using the gt package.
#'
#' @param data Data frame from simulate_impower133()
#' @return A gt table object
#' @export
make_table3 <- function(data) {
  requireNamespace("gt", quietly = TRUE)

  safety_pop <- data[data$received_tx == TRUE, ]
  n_atezo_safety <- sum(safety_pop$treatment == "Atezolizumab")
  n_placebo_safety <- sum(safety_pop$treatment == "Placebo")

  summarize_ae <- function(dat, ae_col, group) {
    sub <- dat[dat$treatment == group, ]
    n <- nrow(sub)
    g12_count <- sum(sub[[ae_col]] == "G1-2")
    g34_count <- sum(sub[[ae_col]] == "G3-4")
    g5_count  <- sum(sub[[ae_col]] == "G5")
    c(g12_count, g34_count, g5_count, n)
  }

  format_ae <- function(counts) {
    c(sprintf("%d (%.1f)", counts[1], 100 * counts[1] / counts[4]),
      sprintf("%d (%.1f)", counts[2], 100 * counts[2] / counts[4]),
      sprintf("%d (%.1f)", counts[3], 100 * counts[3] / counts[4]))
  }

  # Any adverse event
  atezo_any <- summarize_ae(safety_pop, "grade", "Atezolizumab")
  placebo_any <- summarize_ae(safety_pop, "grade", "Placebo")
  row1 <- c(format_ae(atezo_any), format_ae(placebo_any))

  # Subtitle row
  row2 <- c("", "", "", "", "", "")

  # Specific adverse events
  ae_names <- c("Neutropenia","Anemia","Alopecia","Nausea","Fatigue",
                "Decreased neutrophil count","Decreased appetite","Thrombocytopenia",
                "Decreased platelet count","Vomiting","Constipation","Leukopenia",
                "Decreased white-cell count","Diarrhea","Febrile neutropenia",
                "Infusion-related reaction")

  rows_ae <- do.call(rbind, lapply(ae_names, function(name) {
    ae_col <- paste0("AE_", name)
    atezo <- summarize_ae(safety_pop, ae_col, "Atezolizumab")
    placebo <- summarize_ae(safety_pop, ae_col, "Placebo")
    c(format_ae(atezo), format_ae(placebo))
  }))

  table3_values <- rbind(row1, row2, rows_ae)

  table3_final <- data.frame(
    Event = c(
      "Any adverse event",
      "Adverse events with an incidence of >=10% in any grade category or events of grade 3 or 4 with an incidence of >=2% in either group",
      ae_names
    ),
    G12_Atezo = table3_values[, 1],
    G34_Atezo = table3_values[, 2],
    G5_Atezo  = table3_values[, 3],
    G12_Placebo = table3_values[, 4],
    G34_Placebo = table3_values[, 5],
    G5_Placebo  = table3_values[, 6],
    stringsAsFactors = FALSE
  )

  gt::gt(table3_final) %>%
    gt::tab_header(
      title = gt::html("<span style='color:#C00000;'>Table 3.</span> Adverse Events in the As-Treated Population."
    )) %>%
    gt::tab_spanner(
      label = paste0("Atezolizumab Group (N=", n_atezo_safety, ")"),
      columns = c("G12_Atezo", "G34_Atezo", "G5_Atezo")
    ) %>%
    gt::tab_spanner(
      label = paste0("Placebo Group (N=", n_placebo_safety, ")"),
      columns = c("G12_Placebo", "G34_Placebo", "G5_Placebo")
    ) %>%
    gt::cols_label(
      Event = "Event",
      G12_Atezo = gt::html("Grade&nbsp;1&nbsp;or&nbsp;2"),
      G34_Atezo = gt::html("Grade&nbsp;3&nbsp;or&nbsp;4"),
      G5_Atezo  = gt::html("Grade&nbsp;5"),
      G12_Placebo = gt::html("Grade&nbsp;1&nbsp;or&nbsp;2"),
      G34_Placebo = gt::html("Grade&nbsp;3&nbsp;or&nbsp;4"),
      G5_Placebo  = gt::html("Grade&nbsp;5")
    ) %>%
    gt::tab_style(
      style = list(gt::cell_fill(color = "#f4eee2"), gt::cell_text(weight = "bold")),
      locations = gt::cells_title(groups = "title")
    ) %>%
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold")),
      locations = gt::cells_column_labels(gt::everything())
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_spanners()
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#f9f4e8"),
      locations = gt::cells_body(rows = seq(1, nrow(table3_final), 2))
    ) %>%
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(14),
      heading.align = "left"
    )
}
