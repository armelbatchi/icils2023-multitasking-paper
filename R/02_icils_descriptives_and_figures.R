# ICILS descriptives and figures

if (!output_exists("table_desc_multitask.csv")) {
  desc_multitask <- bsg %>%
    filter(!is.na(S_ACMULT), !is.na(TOTWGTS)) %>%
    group_by(CNTRYID, country_label, region) %>%
    summarise(n = n(),
              wt_mean = weighted.mean(S_ACMULT, TOTWGTS, na.rm = TRUE),
              wt_sd = sqrt(Hmisc::wtd.var(S_ACMULT, TOTWGTS, na.rm = TRUE)),
              .groups = "drop") %>%
    arrange(wt_mean)
  write_csv(desc_multitask, file.path(OUT_DIR, "table_desc_multitask.csv"))
} else {
  desc_multitask <- read_csv(file.path(OUT_DIR, "table_desc_multitask.csv"),
                             show_col_types = FALSE)
  cat("Loaded cached table_desc_multitask.csv\n")
}

if (!output_exists("table_desc_cil.csv")) {
  desc_cil <- bsg %>%
    filter(!is.na(PV1CIL), !is.na(TOTWGTS)) %>%
    group_by(CNTRYID, country_label, region) %>%
    summarise(cil_mean = weighted.mean(PV1CIL, TOTWGTS, na.rm = TRUE),
              cil_sd = sqrt(Hmisc::wtd.var(PV1CIL, TOTWGTS, na.rm = TRUE)),
              .groups = "drop")
  write_csv(desc_cil, file.path(OUT_DIR, "table_desc_cil.csv"))
} else {
  desc_cil <- read_csv(file.path(OUT_DIR, "table_desc_cil.csv"), show_col_types = FALSE)
  cat("Loaded cached table_desc_cil.csv\n")
}

if (!output_exists("fig1_ridgeplot.png")) {
  fig1_data <- bsg %>%
    filter(!is.na(S_ACMULT)) %>%
    left_join(desc_multitask %>% select(CNTRYID, wt_mean), by = "CNTRYID") %>%
    mutate(country_label = fct_reorder(country_label, wt_mean))

  fig1 <- ggplot(fig1_data, aes(x = S_ACMULT, y = country_label, fill = region)) +
    geom_density_ridges(scale = 1.35, rel_min_height = .005, alpha = .8,
                        colour = "white", linewidth = .3,
                        quantile_lines = TRUE, quantiles = 2) +
    scale_fill_manual(values = pal_region[names(pal_region) != "Pooled"], name = "Region") +
    scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 10)) +
    labs(title = "Cross-national variation in academic-media multitasking during homework",
         subtitle = "Distribution of S_ACMULT (IRT WLE, M=50, SD=10) by country, ICILS 2023.",
         x = "S_ACMULT", y = NULL,
         caption = "Source: ICILS 2023. Weighted estimates. Vertical lines = weighted median.")
  save_both(fig1, "fig1_ridgeplot", 10, 9)
} else {
  cat("fig1_ridgeplot already exists, skipping.\n")
}

# ── Checks ────────────────────────────────────────────────────────────────────
needed_vars <- c("S_ACMULT", "TOTWGTS", "country_label")
missing_vars <- setdiff(needed_vars, names(bsg))

if (length(missing_vars) > 0) {
  stop(
    paste0(
      "The following variables are missing from bsg: ",
      paste(missing_vars, collapse = ", ")
    )
  )
}

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ── Weighted helpers (base R only) ────────────────────────────────────────────
weighted_quantile <- function(x, w, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(keep)) return(rep(NA_real_, length(probs)))

  x <- x[keep]
  w <- w[keep]

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  cw <- cumsum(w) / sum(w)

  sapply(probs, function(p) {
    idx <- which(cw >= p)[1]
    x[idx]
  })
}

weighted_sd <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  if (!any(keep)) return(NA_real_)

  x <- x[keep]
  w <- w[keep]

  mu <- sum(w * x) / sum(w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}

fmt1 <- function(x) ifelse(is.na(x), "", sprintf("%.1f", x))
fmt0 <- function(x) ifelse(is.na(x), "", format(round(x), big.mark = ",", scientific = FALSE))

# ── Build full descriptive table ──────────────────────────────────────────────
dat <- bsg[is.finite(bsg$S_ACMULT) & is.finite(bsg$TOTWGTS) & bsg$TOTWGTS > 0, , drop = FALSE]

split_list <- split(dat, dat$country_label, drop = TRUE)

table_full <- do.call(
  rbind,
  lapply(names(split_list), function(cty) {
    d <- split_list[[cty]]

    qs <- weighted_quantile(d$S_ACMULT, d$TOTWGTS, probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
    reg <- if ("region" %in% names(d)) {
      vals <- unique(as.character(d$region[!is.na(d$region)]))
      if (length(vals) == 0) "" else vals[1]
    } else {
      ""
    }

    wt_mean <- weighted.mean(d$S_ACMULT, d$TOTWGTS, na.rm = TRUE)
    wt_sd_v <- weighted_sd(d$S_ACMULT, d$TOTWGTS)
    pct_high60 <- 100 * weighted.mean(as.numeric(d$S_ACMULT >= 60), d$TOTWGTS, na.rm = TRUE)
    pct_low40  <- 100 * weighted.mean(as.numeric(d$S_ACMULT <= 40), d$TOTWGTS, na.rm = TRUE)

    data.frame(
      Country = cty,
      Region = reg,
      N_students = nrow(d),
      Weighted_mean = wt_mean,
      Weighted_SD = wt_sd_v,
      Weighted_median = qs[3],
      P10 = qs[1],
      P25 = qs[2],
      P75 = qs[4],
      P90 = qs[5],
      IQR = qs[4] - qs[2],
      Min = min(d$S_ACMULT, na.rm = TRUE),
      Max = max(d$S_ACMULT, na.rm = TRUE),
      Pct_high_60 = pct_high60,
      Pct_low_40 = pct_low40,
      Interpretation = ifelse(
        wt_mean >= 55,
        "Higher average multitasking",
        ifelse(wt_mean < 45, "Lower average multitasking", "Around ICILS average")
      ),
      stringsAsFactors = FALSE
    )
  })
)

row.names(table_full) <- NULL

table_full <- table_full[order(-table_full$Weighted_mean, table_full$Country), ]
table_full$Rank <- seq_len(nrow(table_full))

table_full <- table_full[, c(
  "Rank", "Country", "Region", "N_students", "Weighted_mean", "Weighted_SD",
  "Weighted_median", "P10", "P25", "P75", "P90", "IQR", "Min", "Max",
  "Pct_high_60", "Pct_low_40", "Interpretation"
)]

# Save full CSV
write.csv(
  table_full,
  file = file.path(OUT_DIR, "table_multitask_word_descriptive_full.csv"),
  row.names = FALSE,
  na = ""
)

table_word <- data.frame(
  Rank = table_full$Rank,
  Country = table_full$Country,
  Region = table_full$Region,
  N = fmt0(table_full$N_students),
  Mean = fmt1(table_full$Weighted_mean),
  SD = fmt1(table_full$Weighted_SD),
  Median = fmt1(table_full$Weighted_median),
  `IQR (P25-P75)` = paste0(fmt1(table_full$P25), "–", fmt1(table_full$P75)),
  `P10-P90` = paste0(fmt1(table_full$P10), "–", fmt1(table_full$P90)),
  `% >=60` = fmt1(table_full$Pct_high_60),
  `% <=40` = fmt1(table_full$Pct_low_40),
  Interpretation = table_full$Interpretation,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  table_word,
  file = file.path(OUT_DIR, "table_multitask_word_descriptive_compact.csv"),
  row.names = FALSE,
  na = ""
)

rtf_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x
}

make_rtf_row <- function(values, widths, header = FALSE, font_size = 16) {
  vals <- rtf_escape(values)
  cellx <- cumsum(widths)

  row_def <- paste0(
    "\\trowd\\trgaph70\\trleft0",
    paste0("\\clvertalc\\cellx", cellx, collapse = "")
  )

  if (header) {
    cells <- paste0("\\intbl\\qc\\b\\fs", font_size, " ", vals, "\\b0\\cell")
  } else {
    cells <- paste0("\\intbl\\ql\\fs", font_size, " ", vals, "\\cell")
  }

  paste(c(row_def, cells, "\\row"), collapse = "\n")
}

col_widths <- c(500, 2100, 1100, 700, 700, 700, 800, 1200, 1100, 800, 800, 1900)

rtf_header <- paste0(
  "{\\rtf1\\ansi\\deff0",
  "{\\fonttbl{\\f0 Times New Roman;}}",
  "\\paperw15840\\paperh12240\\margl720\\margr720\\margt720\\margb720\\landscape",
  "\\fs20"
)

title_block <- paste0(
  "\\pard\\qc\\b\\fs24 ",
  rtf_escape("Table X. Cross-national descriptive distribution of academic-media multitasking during homework (S_ACMULT), ICILS 2023"),
  "\\b0\\par\\par ",
  "\\pard\\ql\\fs18 ",
  rtf_escape("This table complements the ridge plot by showing country-specific weighted location, spread, tail values, and the proportions of students with comparatively high or low multitasking."),
  "\\par\\par "
)

header_vals <- names(table_word)
body_rows <- apply(table_word, 1, function(x) make_rtf_row(x, col_widths, header = FALSE, font_size = 16))
header_row <- make_rtf_row(header_vals, col_widths, header = TRUE, font_size = 16)

footer_block <- paste0(
  "\\pard\\ql\\fs16\\par ",
  rtf_escape("Notes: Values are weighted using TOTWGTS. Higher S_ACMULT indicates more frequent academic-media multitasking during homework. Thresholds for % >=60 and % <=40 are descriptive only. A full version with additional columns (including Min-Max) is also saved as CSV."),
  "\\par}",
  collapse = ""
)

rtf_lines <- c(
  rtf_header,
  title_block,
  header_row,
  body_rows,
  footer_block
)

rtf_path <- file.path(OUT_DIR, "Table_multitasking_descriptive_ICILS2023.rtf")
writeLines(rtf_lines, con = rtf_path, useBytes = TRUE)

cat("Saved files:\n")
cat(" - ", file.path(OUT_DIR, "table_multitask_word_descriptive_full.csv"), "\n", sep = "")
cat(" - ", file.path(OUT_DIR, "table_multitask_word_descriptive_compact.csv"), "\n", sep = "")
cat(" - ", rtf_path, "\n", sep = "")

if (!output_exists("fig2_scatter.png")) {
  fig2_data <- inner_join(desc_multitask, desc_cil,
                          by = c("CNTRYID", "country_label", "region"))
  fig2 <- ggplot(fig2_data, aes(wt_mean, cil_mean)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "#457B9D", fill = "#457B9D", alpha = .15, linetype = "dashed") +
    geom_point(aes(colour = region, size = n), alpha = .85) +
    geom_text_repel(aes(label = country_label, colour = region),
                    size = 3.1, fontface = "bold", max.overlaps = 25, show.legend = FALSE) +
    scale_colour_manual(values = pal_region[names(pal_region) != "Pooled"], name = "Region") +
    scale_size_continuous(range = c(2, 8), labels = comma, name = "Sample size") +
    labs(title = "National multitasking versus CIL",
         subtitle = "Country-level weighted means: S_ACMULT vs CIL, ICILS 2023.",
         x = "Mean S_ACMULT", y = "Mean CIL",
         caption = "Source: ICILS 2023. CIL based on PV1. Dashed line = OLS fit.")
  save_both(fig2, "fig2_scatter", 9.5, 7.5)
} else {
  cat("fig2_scatter already exists, skipping.\n")
}
