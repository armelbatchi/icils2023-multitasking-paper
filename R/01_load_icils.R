# Load and prepare ICILS 2023 data

sav_files <- require_data_files(
  path = ICILS_DIR,
  pattern = "^BSG.*\\.sav$",
  study_label = "ICILS 2023",
  repository_url = "https://www.iea.nl/data-tools/repository/icils"
)
cat(sprintf("Reading %d ICILS SPSS files...\n", length(sav_files)))

bsg_raw <- bind_rows(lapply(sav_files, function(f) {
  d <- tryCatch(read_spss_file(f), error = function(e) {
    warning("Failed to read: ", f, " - ", e$message); NULL
  })
  if (is.null(d)) return(NULL)
  names(d) <- toupper(names(d)); d
}))

bsg <- bsg_raw
bsg$CNTRYID <- as.character(bsg$IDCNTRY)
bsg <- left_join(bsg, region_map, by = "CNTRYID") %>%
  mutate(country_label = coalesce(country_label, paste0("Country_", CNTRYID)),
         region = coalesce(region, "Other"))

# Canonical covariates
SES_VAR    <- pick_first(c("S_NISB", "S_IISB", "NISB"), names(bsg))
GENDER_VAR <- pick_first(c("S_GENDER", "SGENDER"), names(bsg))
EXCOMP_VAR <- pick_first(c("S_EXCOMP"), names(bsg))

bsg$SES_CAN    <- if (!is.null(SES_VAR)) as.numeric(bsg[[SES_VAR]]) else NA_real_
bsg$GENDER_CAN <- if (!is.null(GENDER_VAR)) as.numeric(bsg[[GENDER_VAR]]) else NA_real_
bsg$EXCOMP_CAN <- if (!is.null(EXCOMP_VAR)) as.numeric(bsg[[EXCOMP_VAR]]) else NA_real_

available_covs <- c("SES_CAN", "EXCOMP_CAN", "GENDER_CAN")
available_covs <- available_covs[sapply(available_covs, function(v)
  v %in% names(bsg) && sum(!is.na(bsg[[v]])) > 100)]

countries <- sort(unique(bsg$CNTRYID))
cat(sprintf("Loaded: %s students, %d countries, covariates: %s\n",
            format(nrow(bsg), big.mark = ","), length(countries),
            paste(available_covs, collapse = ", ")))
