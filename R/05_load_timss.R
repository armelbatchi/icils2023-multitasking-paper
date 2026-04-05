# Load and prepare TIMSS 2023 data

all_timss <- require_data_files(
  path = TIMSS_DIR,
  pattern = "\\.sav$",
  study_label = "TIMSS 2023",
  repository_url = "https://www.iea.nl/data-tools/repository/timss"
)

# File groups
bsg_files <- all_timss[grepl("^bsg", basename(all_timss), ignore.case = TRUE)]
btm_files <- all_timss[grepl("^btm", basename(all_timss), ignore.case = TRUE)]
bts_files <- all_timss[grepl("^bts", basename(all_timss), ignore.case = TRUE)]

cat(sprintf("TIMSS files found: %d BSG (student), %d BTM (math teacher), %d BTS (sci teacher)\n",
            length(bsg_files), length(btm_files), length(bts_files)))

if (length(bsg_files) == 0) {
  stop(
    paste0(
      "No TIMSS student background files were found.\n",
      "Download the Grade 8 SPSS files from the official TIMSS repository and place them in the local TIMSS data folder."
    ),
    call. = FALSE
  )
}

cat(sprintf("Reading %d student background files...\n", length(bsg_files)))
tsg <- bind_rows(lapply(bsg_files, function(f) {
  d <- tryCatch(read_spss_file(f), error = function(e) {
    warning("Failed: ", basename(f), " - ", e$message); NULL
  })
  if (is.null(d)) return(NULL)
  names(d) <- toupper(names(d)); d
}))

cat(sprintf("Student data loaded: %s rows, %d columns\n",
            format(nrow(tsg), big.mark = ","), ncol(tsg)))

tsg$CNTRYID <- as.character(tsg$IDCNTRY)

# Prong 1: learning games frequency
if ("BSBG12F" %in% names(tsg)) {
    # Reverse coding so higher values mean more frequent use
  tsg$GAME_FREQ    <- 5 - as.numeric(tsg$BSBG12F)
  tsg$GAME_FREQ_z  <- as.numeric(scale(tsg$GAME_FREQ))
  tsg$GAME_MONTHLY <- as.integer(as.numeric(tsg$BSBG12F) <= 3)  # at least monthly
  cat("  Prong 1 (BSBG12F): OK -", sum(!is.na(tsg$GAME_FREQ_z)), "non-missing\n")
} else {
  cat("  WARNING: BSBG12F not found in student data!\n")
}

# Prong 2: internet use for schoolwork composite
inet_items <- paste0("BSBG12", LETTERS[1:6])
inet_avail <- intersect(inet_items, names(tsg))
cat(sprintf("  Prong 2 items found: %d of 6 (%s)\n",
            length(inet_avail), paste(inet_avail, collapse = ", ")))

if (length(inet_avail) >= 3) {
    # Reverse coding so higher values mean more frequent use
  for (v in inet_avail) tsg[[paste0(v, "_r")]] <- 5 - as.numeric(tsg[[v]])
  inet_r <- paste0(inet_avail, "_r")
  tsg$INET_SCHOOL <- rowMeans(tsg[, inet_r], na.rm = FALSE)
    # Allow one missing item
  nv <- rowSums(!is.na(tsg[, inet_r]))
  tsg$INET_SCHOOL[nv < length(inet_r) - 1] <- NA
  tsg$INET_SCHOOL_z <- as.numeric(scale(tsg$INET_SCHOOL))
  tsg$INET_HIGH <- as.integer(
    tsg$INET_SCHOOL >= quantile(tsg$INET_SCHOOL, .67, na.rm = TRUE))
  cat("  Prong 2 (INET_SCHOOL): OK -", sum(!is.na(tsg$INET_SCHOOL_z)), "non-missing\n")
}

# Prong 3: teacher distraction climate
if (length(btm_files) > 0 || length(bts_files) > 0) {
  cat("  Loading teacher files for Prong 3...\n")

  load_teacher <- function(files, prefix) {
    if (length(files) == 0) return(NULL)
    d <- bind_rows(lapply(files, function(f) {
      x <- tryCatch(read_spss_file(f), error = function(e) NULL)
      if (is.null(x)) return(NULL)
      names(x) <- toupper(names(x)); x
    }))
    cat(sprintf("    %s: %s rows\n", prefix, format(nrow(d), big.mark = ",")))
    d
  }

  btm_df <- load_teacher(btm_files, "BTM (math)")
  bts_df <- load_teacher(bts_files, "BTS (science)")

    # Aggregate teacher indicators to the school level
  teacher_vars <- list()

  if (!is.null(btm_df) && "BTBG13G" %in% names(btm_df)) {
    teacher_vars$btbg13g_math <- btm_df %>%
      mutate(BTBG13G = as.numeric(BTBG13G)) %>%
      group_by(IDCNTRY, IDSCHOOL) %>%
      summarise(BTBG13G_math_sch = mean(BTBG13G, na.rm = TRUE), .groups = "drop")
  }

  if (!is.null(btm_df) && "BTBM18C" %in% names(btm_df)) {
    teacher_vars$btbm18c <- btm_df %>%
      mutate(BTBM18C = as.numeric(BTBM18C)) %>%
      group_by(IDCNTRY, IDSCHOOL) %>%
      summarise(BTBM18C_sch = mean(BTBM18C, na.rm = TRUE), .groups = "drop")
  }

  if (!is.null(bts_df) && "BTBS21C" %in% names(bts_df)) {
    teacher_vars$btbs21c <- bts_df %>%
      mutate(BTBS21C = as.numeric(BTBS21C)) %>%
      group_by(IDCNTRY, IDSCHOOL) %>%
      summarise(BTBS21C_sch = mean(BTBS21C, na.rm = TRUE), .groups = "drop")
  }

  if (!is.null(bts_df) && "BTBG13G" %in% names(bts_df)) {
    teacher_vars$btbg13g_sci <- bts_df %>%
      mutate(BTBG13G = as.numeric(BTBG13G)) %>%
      group_by(IDCNTRY, IDSCHOOL) %>%
      summarise(BTBG13G_sci_sch = mean(BTBG13G, na.rm = TRUE), .groups = "drop")
  }

  # Merge all teacher school-level variables onto student data
  for (nm in names(teacher_vars)) {
    tsg <- left_join(tsg, teacher_vars[[nm]], by = c("IDCNTRY", "IDSCHOOL"))
  }

  # Build distraction climate composite
  dist_cols <- intersect(c("BTBG13G_math_sch", "BTBM18C_sch",
                           "BTBS21C_sch", "BTBG13G_sci_sch"), names(tsg))
  if (length(dist_cols) >= 1) {
    for (dc in dist_cols) tsg[[paste0(dc, "_z")]] <- as.numeric(scale(tsg[[dc]]))
    tsg$DISTRACT_CLIMATE_z <- rowMeans(
      tsg[, paste0(dist_cols, "_z"), drop = FALSE], na.rm = TRUE)
    tsg$DISTRACT_CLIMATE_z[rowSums(!is.na(tsg[, paste0(dist_cols, "_z")])) == 0] <- NA
    cat(sprintf("  Prong 3 (DISTRACT_CLIMATE): OK - %d non-missing (from %d items)\n",
                sum(!is.na(tsg$DISTRACT_CLIMATE_z)), length(dist_cols)))
  }

  # Clean up large teacher objects
  rm(btm_df, bts_df, teacher_vars); gc()

} else {
  cat("  No teacher files found. Prong 3 will be skipped.\n")
  cat("  Teacher files were not found. Prong 3 will be skipped unless the TIMSS teacher files are added locally.\n")
}

# Covariates
if ("BSBGSEC" %in% names(tsg)) tsg$DIGSE_z <- as.numeric(scale(as.numeric(tsg$BSBGSEC)))
if ("ITSEX" %in% names(tsg))   tsg$FEMALE  <- ifelse(as.numeric(tsg$ITSEX) == 2, 1, 0)
if ("BSBGHER" %in% names(tsg)) tsg$SES_z   <- as.numeric(scale(as.numeric(tsg$BSBGHER)))

timss_covs <- c()
for (v in c("DIGSE_z", "FEMALE", "SES_z")) {
  if (v %in% names(tsg) && sum(!is.na(tsg[[v]])) > 1000)
    timss_covs <- c(timss_covs, v)
}

timss_countries <- sort(unique(tsg$CNTRYID))
cat(sprintf("\nTIMSS SUMMARY: %s students, %d countries, covariates: %s\n",
            format(nrow(tsg), big.mark = ","), length(timss_countries),
            paste(timss_covs, collapse = ", ")))

cat(sprintf("  PVs: BSMMAT01=%d non-missing, BSSSCI01=%d non-missing\n",
            sum(!is.na(tsg$BSMMAT01)), sum(!is.na(tsg$BSSSCI01))))
cat(sprintf("  Weights: TOTWGT=%d non-missing, JKZONE=%d, JKREP=%d\n",
            sum(!is.na(tsg$TOTWGT)), sum(!is.na(tsg$JKZONE)), sum(!is.na(tsg$JKREP))))

# TIMSS country map
timss_country_map <- tribble(
  ~CNTRYID, ~country_label,
  "36",   "Australia",
  "40",   "Austria",
  "31",   "Azerbaijan",
  "48",   "Bahrain",
  "76",   "Brazil",
  "152",  "Chile",
  "158",  "Chinese Taipei",
  "196",  "Cyprus",
  "203",  "Czech Republic",
  "246",  "Finland",
  "250",  "France",
  "268",  "Georgia",
  "275",  "Palestine",
  "344",  "Hong Kong SAR",
  "348",  "Hungary",
  "364",  "Iran",
  "372",  "Ireland",
  "376",  "Israel",
  "380",  "Italy",
  "384",  "Côte d'Ivoire",
  "392",  "Japan",
  "398",  "Kazakhstan",
  "400",  "Jordan",
  "410",  "Korea",
  "414",  "Kuwait",
  "440",  "Lithuania",
  "458",  "Malaysia",
  "470",  "Malta",
  "504",  "Morocco",
  "512",  "Oman",
  "554",  "New Zealand",
  "620",  "Portugal",
  "634",  "Qatar",
  "642",  "Romania",
  "682",  "Saudi Arabia",
  "702",  "Singapore",
  "710",  "South Africa",
  "752",  "Sweden",
  "784",  "United Arab Emirates",
  "7841", "Abu Dhabi (UAE)",
  "7842", "Dubai (UAE)",
  "7843", "Sharjah (UAE)",
  "792",  "Türkiye",
  "840",  "United States",
  "860",  "Uzbekistan",
  "926",  "England"
)

label_timss <- function(df) {
  df$CNTRYID <- as.character(df$CNTRYID)
  df <- left_join(df, timss_country_map, by = "CNTRYID")
  df$country_label <- coalesce(df$country_label, df$CNTRYID)
  df
}
