#' Simulate IMpower133 complete ITT population
#'
#' Simulate 403 virtual patients' baseline characteristics, efficacy outcomes,
#' and safety data based on published summary statistics from Tables 1-3.
#'
#' @param seed Random seed, default 21
#' @param os_path Path to reconstructed OS data CSV (optional)
#' @param pfs_path Path to reconstructed PFS data CSV (optional)
#' @param method Matching method: "random" or "risk_score" (default)
#' @param output_path Path to save output CSV
#' @return A data.frame with 403 rows and all variables (invisibly)
#' @export
#'
#' @examples
#' \donttest{
#' df <- simulate_impower133()
#' head(df)
#' }
simulate_impower133 <- function(seed = 21,
                                os_path = NULL,
                                pfs_path = NULL,
                                method = c("risk_score", "random"),
                                output_path = tempfile(fileext = ".csv")) {
  method <- match.arg(method)
  set.seed(seed)

  n_atezo <- 201
  n_placebo <- 202
  n_total <- n_atezo + n_placebo

  treatment <- c(rep("Atezolizumab", n_atezo), rep("Placebo", n_placebo))

  # ========== Table 1 variables ==========

  age <- c(
    round(rtriang(n_atezo, 28, 64, 90), 0),
    round(rtriang(n_placebo, 26, 64, 87), 0)
  )
  age_group <- ifelse(age < 65, "<65 yr", ">=65 yr")

  sex <- c(
    sample(c("Male", "Female"), n_atezo, replace = TRUE, prob = c(0.642, 0.358)),
    sample(c("Male", "Female"), n_placebo, replace = TRUE, prob = c(0.653, 0.347))
  )

  ecog <- c(
    sample(c(0, 1), n_atezo, replace = TRUE, prob = c(0.363, 0.637)),
    sample(c(0, 1), n_placebo, replace = TRUE, prob = c(0.332, 0.668))
  )

  smoking <- c(
    sample(c("Never smoked", "Current smoker", "Former smoker"), n_atezo,
           replace = TRUE, prob = c(0.045, 0.368, 0.587)),
    sample(c("Never smoked", "Current smoker", "Former smoker"), n_placebo,
           replace = TRUE, prob = c(0.015, 0.371, 0.614))
  )

  brain_mets <- c(
    sample(c("Yes", "No"), n_atezo, replace = TRUE, prob = c(0.085, 0.915)),
    sample(c("Yes", "No"), n_placebo, replace = TRUE, prob = c(0.089, 0.911))
  )

  tmb_avail <- c(
    sample(c(TRUE, FALSE), n_atezo, replace = TRUE, prob = c(173/201, 28/201)),
    sample(c(TRUE, FALSE), n_placebo, replace = TRUE, prob = c(178/202, 24/202))
  )

  tmb_10 <- rep(NA_character_, n_total)
  tmb_10[treatment == "Atezolizumab" & tmb_avail] <-
    sample(c("<10 mut/Mb", ">=10 mut/Mb"), sum(tmb_avail & treatment == "Atezolizumab"),
           replace = TRUE, prob = c(0.410, 0.590))
  tmb_10[treatment == "Placebo" & tmb_avail] <-
    sample(c("<10 mut/Mb", ">=10 mut/Mb"), sum(tmb_avail & treatment == "Placebo"),
           replace = TRUE, prob = c(0.382, 0.618))

  tmb_16 <- rep(NA_character_, n_total)
  tmb_16[treatment == "Atezolizumab" & tmb_avail] <-
    sample(c("<16 mut/Mb", ">=16 mut/Mb"), sum(tmb_avail & treatment == "Atezolizumab"),
           replace = TRUE, prob = c(0.769, 0.231))
  tmb_16[treatment == "Placebo" & tmb_avail] <-
    sample(c("<16 mut/Mb", ">=16 mut/Mb"), sum(tmb_avail & treatment == "Placebo"),
           replace = TRUE, prob = c(0.775, 0.225))

  target_lesion <- c(
    round(rtriang(n_atezo, 12, 113, 325), 1),
    round(rtriang(n_placebo, 15, 105.5, 353), 1)
  )

  prior_chemo <- c(
    sample(c("Yes", "No"), n_atezo, replace = TRUE, prob = c(0.04, 0.96)),
    sample(c("Yes", "No"), n_placebo, replace = TRUE, prob = c(0.059, 0.941))
  )
  prior_rt <- c(
    sample(c("Yes", "No"), n_atezo, replace = TRUE, prob = c(0.124, 0.876)),
    sample(c("Yes", "No"), n_placebo, replace = TRUE, prob = c(0.139, 0.861))
  )
  prior_surgery <- c(
    sample(c("Yes", "No"), n_atezo, replace = TRUE, prob = c(0.164, 0.836)),
    sample(c("Yes", "No"), n_placebo, replace = TRUE, prob = c(0.124, 0.876))
  )

  # ========== Table 2 variables ==========

  measurable <- stats::rbinom(n_total, 1, 0.95)

  best_response <- rep(NA_character_, n_total)
  n_atezo_meas <- sum(measurable == 1 & treatment == "Atezolizumab")
  n_placebo_meas <- sum(measurable == 1 & treatment == "Placebo")

  best_response[measurable == 1 & treatment == "Atezolizumab"] <- sample(
    c("CR", "PR", "SD", "PD", "NE"), n_atezo_meas, replace = TRUE,
    prob = c(5/201, 116/201, 42/201, 22/201, 16/201))
  best_response[measurable == 1 & treatment == "Placebo"] <- sample(
    c("CR", "PR", "SD", "PD", "NE"), n_placebo_meas, replace = TRUE,
    prob = c(2/202, 128/202, 43/202, 14/202, 15/202))

  responder <- best_response %in% c("CR", "PR")

  dor <- rep(NA_real_, n_total)
  n_resp_atezo <- sum(responder & treatment == "Atezolizumab", na.rm = TRUE)
  n_resp_placebo <- sum(responder & treatment == "Placebo", na.rm = TRUE)
  dor[responder & treatment == "Atezolizumab"] <- round(rtriang(n_resp_atezo, 1.4, 4.2, 19.5), 1)
  dor[responder & treatment == "Placebo"] <- round(rtriang(n_resp_placebo, 2.0, 3.9, 16.1), 1)

  ongoing <- rep(NA, n_total)
  ongoing[responder & treatment == "Atezolizumab"] <- sample(
    c(TRUE, FALSE), n_resp_atezo, replace = TRUE, prob = c(18/121, 1-18/121))
  ongoing[responder & treatment == "Placebo"] <- sample(
    c(TRUE, FALSE), n_resp_placebo, replace = TRUE, prob = c(7/130, 1-7/130))

  # ========== Table 3 variables ==========

  received_tx <- c(
    sample(c(TRUE, FALSE), n_atezo, replace = TRUE, prob = c(198/201, 3/201)),
    sample(c(TRUE, FALSE), n_placebo, replace = TRUE, prob = c(196/202, 6/202))
  )

  n_tx_atezo <- sum(received_tx & treatment == "Atezolizumab")
  n_tx_placebo <- sum(received_tx & treatment == "Placebo")

  any_ae <- c(
    sample(c(1, 0), n_tx_atezo, replace = TRUE, prob = c(0.949, 0.051)),
    sample(c(1, 0), n_tx_placebo, replace = TRUE, prob = c(0.923, 0.077))
  )

  grade <- character(n_tx_atezo + n_tx_placebo)
  tx_grp <- c(rep("Atezolizumab", n_tx_atezo), rep("Placebo", n_tx_placebo))
  for(i in seq_along(grade)) {
    if(any_ae[i] == 1) {
      if(tx_grp[i] == "Atezolizumab") {
        grade[i] <- sample(c("G1-2", "G3-4", "G5"), 1,
                           prob = c(0.369, 0.566, 0.015))
      } else {
        grade[i] <- sample(c("G1-2", "G3-4", "G5"), 1,
                           prob = c(0.347, 0.561, 0.015))
      }
    } else {
      grade[i] <- "None"
    }
  }

  ae_names <- c("Neutropenia","Anemia","Alopecia","Nausea","Fatigue",
                "Decreased neutrophil count","Decreased appetite","Thrombocytopenia",
                "Decreased platelet count","Vomiting","Constipation","Leukopenia",
                "Decreased white-cell count","Diarrhea","Febrile neutropenia",
                "Infusion-related reaction")

  ae_probs <- list(
    atezo_any  = c(0.359,0.389,0.348,0.318,0.212,0.176,0.207,0.162,0.121,0.136,0.101,0.126,0.081,0.096,0.030,0.051),
    atezo_g34  = c(0.227,0.141,0,0.005,0.015,0.141,0.010,0.101,0.035,0.010,0.005,0.051,0.030,0.020,0.030,0.020),
    atezo_g5   = c(0.005,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    placebo_any= c(0.347,0.332,0.337,0.301,0.194,0.230,0.133,0.148,0.143,0.112,0.128,0.092,0.128,0.097,0.061,0.051),
    placebo_g34= c(0.245,0.122,0,0.005,0.005,0.168,0,0.077,0.036,0.015,0,0.041,0.046,0.005,0.061,0.005),
    placebo_g5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  )

  ae_cols <- list()
  for(j in seq_along(ae_names)) {
    ae_col <- paste0("AE_", ae_names[j])
    ae_vals <- character(n_tx_atezo + n_tx_placebo)

    for(i in seq_along(ae_vals)) {
      if(any_ae[i] == 1) {
        grp <- tx_grp[i]
        if(grp == "Atezolizumab") {
          if(stats::rbinom(1, 1, ae_probs$atezo_any[j]) == 1) {
            g12_prob <- ae_probs$atezo_any[j] - ae_probs$atezo_g34[j] - ae_probs$atezo_g5[j]
            ae_vals[i] <- sample(c("G1-2", "G3-4", "G5"), 1,
                                 prob = c(g12_prob, ae_probs$atezo_g34[j], ae_probs$atezo_g5[j]))
          }
        } else {
          if(stats::rbinom(1, 1, ae_probs$placebo_any[j]) == 1) {
            g12_prob <- ae_probs$placebo_any[j] - ae_probs$placebo_g34[j] - ae_probs$placebo_g5[j]
            ae_vals[i] <- sample(c("G1-2", "G3-4", "G5"), 1,
                                 prob = c(g12_prob, ae_probs$placebo_g34[j], ae_probs$placebo_g5[j]))
          }
        }
      }
    }
    ae_cols[[ae_col]] <- ae_vals
  }

  # ========== Build base data frame ==========
  df <- data.frame(
    treatment = factor(treatment, levels = c("Atezolizumab", "Placebo")),
    age = age,
    age_group = factor(age_group, levels = c("<65 yr", ">=65 yr")),
    sex = factor(sex, levels = c("Male", "Female")),
    ecog = factor(ecog, levels = c(0, 1)),
    smoking = factor(smoking, levels = c("Never smoked", "Current smoker", "Former smoker")),
    brain_mets = factor(brain_mets, levels = c("Yes", "No")),
    tmb_10 = factor(tmb_10, levels = c("<10 mut/Mb", ">=10 mut/Mb")),
    tmb_16 = factor(tmb_16, levels = c("<16 mut/Mb", ">=16 mut/Mb")),
    target_lesion = target_lesion,
    prior_chemo = factor(prior_chemo, levels = c("Yes", "No")),
    prior_rt = factor(prior_rt, levels = c("Yes", "No")),
    prior_surgery = factor(prior_surgery, levels = c("Yes", "No")),
    measurable = measurable,
    best_response = factor(best_response, levels = c("CR", "PR", "SD", "PD", "NE")),
    responder = responder,
    dor = dor,
    ongoing = ongoing,
    received_tx = received_tx,
    stringsAsFactors = FALSE
  )

  df$grade <- "None"
  tx_rows <- which(received_tx)
  df$grade[tx_rows] <- grade

  for(name in ae_names) {
    ae_col <- paste0("AE_", name)
    df[[ae_col]] <- "0"
    df[[ae_col]][tx_rows] <- ae_cols[[ae_col]]
  }

  # ========== Merge OS data ==========
  if (!is.null(os_path)) {
    os_data <- utils::read.csv(os_path, header = TRUE)

    df$treatment <- factor(df$treatment, levels = c("Placebo", "Atezolizumab"))

    atezo_os <- os_data[os_data$Treatment_group == "Atezolizumab", ]
    placebo_os <- os_data[os_data$Treatment_group == "Placebo", ]

    n_atezo_tx <- sum(df$treatment == "Atezolizumab")
    n_placebo_tx <- sum(df$treatment == "Placebo")

    if (method == "risk_score") {
      risk_score <- rep(0, nrow(df))
      risk_score <- risk_score + (ecog == 1) * 1.0
      risk_score <- risk_score + (brain_mets == "Yes") * 1.0
      risk_score <- risk_score + (age > 65) * 0.5
      risk_score <- risk_score + (tmb_10 == "<10 mut/Mb" & !is.na(tmb_10)) * 0.5

      placebo_os <- placebo_os[order(placebo_os$Survival_time), ]
      atezo_os <- atezo_os[order(atezo_os$Survival_time), ]

      placebo_rows <- which(df$treatment == "Placebo")
      placebo_rank <- order(risk_score[placebo_rows], decreasing = TRUE)

      atezo_rows <- which(df$treatment == "Atezolizumab")
      atezo_rank <- order(risk_score[atezo_rows], decreasing = TRUE)

      df$os_time <- NA
      df$os_status <- NA
      df$os_time[placebo_rows][placebo_rank] <- placebo_os$Survival_time
      df$os_status[placebo_rows][placebo_rank] <- placebo_os$Status
      df$os_time[atezo_rows][atezo_rank] <- atezo_os$Survival_time
      df$os_status[atezo_rows][atezo_rank] <- atezo_os$Status
    } else {
      atezo_os_idx <- sample(1:nrow(atezo_os), n_atezo_tx, replace = FALSE)
      placebo_os_idx <- sample(1:nrow(placebo_os), n_placebo_tx, replace = FALSE)

      df$os_time <- NA
      df$os_status <- NA
      df$os_time[df$treatment == "Placebo"] <- placebo_os$Survival_time[placebo_os_idx]
      df$os_status[df$treatment == "Placebo"] <- placebo_os$Status[placebo_os_idx]
      df$os_time[df$treatment == "Atezolizumab"] <- atezo_os$Survival_time[atezo_os_idx]
      df$os_status[df$treatment == "Atezolizumab"] <- atezo_os$Status[atezo_os_idx]
    }

    message("Merged OS data (", nrow(os_data), " rows)")
  }

  # ========== Merge PFS data ==========
  if (!is.null(pfs_path)) {
    pfs_data <- utils::read.csv(pfs_path, header = TRUE)

    atezo_pfs <- pfs_data[pfs_data$Treatment_group == "Atezolizumab", ]
    placebo_pfs <- pfs_data[pfs_data$Treatment_group == "Placebo", ]

    if (method == "risk_score") {
      risk_score <- rep(0, nrow(df))
      risk_score <- risk_score + (ecog == 1) * 1.0
      risk_score <- risk_score + (brain_mets == "Yes") * 1.0
      risk_score <- risk_score + (age > 65) * 0.5
      risk_score <- risk_score + (tmb_10 == "<10 mut/Mb" & !is.na(tmb_10)) * 0.5

      placebo_pfs <- placebo_pfs[order(placebo_pfs$Survival_time), ]
      atezo_pfs <- atezo_pfs[order(atezo_pfs$Survival_time), ]

      placebo_rows <- which(df$treatment == "Placebo")
      placebo_rank <- order(risk_score[placebo_rows], decreasing = TRUE)

      atezo_rows <- which(df$treatment == "Atezolizumab")
      atezo_rank <- order(risk_score[atezo_rows], decreasing = TRUE)

      df$pfs_time <- NA
      df$pfs_status <- NA
      df$pfs_time[placebo_rows][placebo_rank] <- placebo_pfs$Survival_time
      df$pfs_status[placebo_rows][placebo_rank] <- placebo_pfs$Status
      df$pfs_time[atezo_rows][atezo_rank] <- atezo_pfs$Survival_time
      df$pfs_status[atezo_rows][atezo_rank] <- atezo_pfs$Status
    } else {
      atezo_pfs_idx <- sample(1:nrow(atezo_pfs), sum(df$treatment == "Atezolizumab"), replace = FALSE)
      placebo_pfs_idx <- sample(1:nrow(placebo_pfs), sum(df$treatment == "Placebo"), replace = FALSE)

      df$pfs_time <- NA
      df$pfs_status <- NA
      df$pfs_time[df$treatment == "Placebo"] <- placebo_pfs$Survival_time[placebo_pfs_idx]
      df$pfs_status[df$treatment == "Placebo"] <- placebo_pfs$Status[placebo_pfs_idx]
      df$pfs_time[df$treatment == "Atezolizumab"] <- atezo_pfs$Survival_time[atezo_pfs_idx]
      df$pfs_status[df$treatment == "Atezolizumab"] <- atezo_pfs$Status[atezo_pfs_idx]
    }

    message("Merged PFS data (", nrow(pfs_data), " rows)")
  }

  df$treatment <- factor(df$treatment, levels = c("Atezolizumab", "Placebo"))

  # ========== Save ==========
  utils::write.csv(df, output_path, row.names = FALSE, fileEncoding = "UTF-8")
  message("Complete simulated data saved to ", output_path)

  return(invisible(df))
}
