# Setup, configuration, and shared helpers

# Package setup
pkgs <- c(
  "dplyr", "tidyr", "purrr", "stringr", "forcats", "readr", "tibble",
  "survey", "mitools", "splines", "mgcv", "metafor",
  "ggplot2", "ggrepel", "patchwork", "scales", "ggridges",
  "RColorBrewer", "viridis", "ggtext", "ggh4x", "Hmisc",
  "gt", "lme4", "lmerTest", "broom.mixed", "quantreg"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

# Fallback SPSS reader setup
HAVE_HAVEN <- tryCatch({
  suppressPackageStartupMessages(library(haven))
  TRUE
}, error = function(e) {
  message("haven binary failed: ", e$message)
  message("Attempting source install...")
  ok <- tryCatch({
    install.packages("haven", type = "source", repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(haven))
    TRUE
  }, error = function(e2) {
    message("Source install also failed. Using foreign::read.spss() instead.")
    if (!requireNamespace("foreign", quietly = TRUE))
      install.packages("foreign", repos = "https://cloud.r-project.org")
    library(foreign)
    FALSE
  })
  ok
})

HAVE_LME4     <- requireNamespace("lme4", quietly = TRUE)
HAVE_QUANTREG <- requireNamespace("quantreg", quietly = TRUE)

cat("\n--- Package status ---\n")
cat("  R version :", R.version.string, "\n")
cat("  haven     :", HAVE_HAVEN, "\n")
cat("  lme4      :", HAVE_LME4, "\n")
cat("  quantreg  :", HAVE_QUANTREG, "\n\n")

# Shared SPSS reader
read_spss_file <- function(f) {
  if (HAVE_HAVEN) {
    haven::read_sav(f)
  } else {
    foreign::read.spss(f, to.data.frame = TRUE, use.value.labels = FALSE)
  }
}

DATA_DIR  <- Sys.getenv("IEA_DATA_DIR", unset = ".local_data")
ICILS_DIR <- Sys.getenv("ICILS_DIR", unset = file.path(DATA_DIR, "icils"))
TIMSS_DIR <- Sys.getenv("TIMSS_DIR", unset = file.path(DATA_DIR, "timss"))
OUT_DIR   <- Sys.getenv("OUT_DIR", unset = ".local_output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

require_data_files <- function(path, pattern, study_label, repository_url) {
  files <- list.files(path, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop(
      paste0(
        "No matching source files were found for ", study_label, ".\n",
        "Download the data from the official IEA repository and place the files in a local folder of your choice.\n",
        "Repository: ", repository_url, "\n",
        "Expected folder: ", path
      ),
      call. = FALSE
    )
  }
  files
}

# Constants
PV_CIL  <- paste0("PV", 1:5, "CIL")
PV_CT   <- paste0("PV", 1:5, "CT")
PV_MATH <- paste0("BSMMAT0", 1:5)
PV_SCI  <- paste0("BSSSCI0", 1:5)
N_REP   <- 75
JRR_FAC <- 0.5

# Cache helper
output_exists <- function(...) {
  all(file.exists(file.path(OUT_DIR, c(...))))
}

# ICILS country map
region_map <- tribble(
  ~CNTRYID, ~country_label,           ~region,
  "31",     "Azerbaijan",             "Central Asia",
  "40",     "Austria",                "Europe",
  "70",     "Bosnia & Herzegovina",   "Southeast Europe",
  "152",    "Chile",                  "Latin America",
  "158",    "Chinese Taipei",         "East Asia",
  "191",    "Croatia",                "Europe",
  "196",    "Cyprus",                 "Europe",
  "203",    "Czechia",                "Europe",
  "208",    "Denmark",                "Europe",
  "246",    "Finland",                "Europe",
  "250",    "France",                 "Europe",
  "276",    "Germany",                "Europe",
  "276001", "Germany (NRW)",          "Europe",
  "300",    "Greece",                 "Europe",
  "348",    "Hungary",                "Europe",
  "380",    "Italy",                  "Europe",
  "398",    "Kazakhstan",             "Central Asia",
  "410",    "Republic of Korea",      "East Asia",
  "411",    "Korea (benchmarking)",   "East Asia",
  "428",    "Latvia",                 "Europe",
  "442",    "Luxembourg",             "Europe",
  "470",    "Malta",                  "Europe",
  "512",    "Oman",                   "Middle East",
  "528",    "Netherlands",            "Europe",
  "578",    "Norway",                 "Europe",
  "620",    "Portugal",               "Europe",
  "688",    "Serbia",                 "Southeast Europe",
  "703",    "Slovak Republic",        "Europe",
  "705",    "Slovenia",               "Europe",
  "724",    "Spain",                  "Europe",
  "752",    "Sweden",                 "Europe",
  "840",    "United States",          "North America",
  "858",    "Uruguay",                "Latin America",
  "956",    "Kosovo",                 "Southeast Europe",
  "9642",   "Belgium (Flemish)",      "Europe"
)

# Plot theme and palette
pal_region <- c(
  "East Asia" = "#E63946", "Europe" = "#457B9D", "Latin America" = "#2A9D8F",
  "North America" = "#E9C46A", "Central Asia" = "#F4A261",
  "Middle East" = "#9B5DE5", "Southeast Europe" = "#00BBF9", "Pooled" = "black"
)

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "sans", colour = "grey20"),
      plot.title       = ggtext::element_textbox_simple(
        size = base_size + 3, face = "bold", lineheight = 1.15, margin = margin(b = 8)),
      plot.subtitle    = ggtext::element_textbox_simple(
        size = base_size, colour = "grey40", lineheight = 1.10, margin = margin(b = 10)),
      plot.caption     = element_text(size = base_size - 2, colour = "grey45",
                                      hjust = 0, margin = margin(t = 10)),
      axis.title       = element_text(size = base_size, face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      plot.margin      = margin(12, 12, 12, 12)
    )
}
theme_set(theme_paper())

# Helper functions
save_both <- function(p, stub, w, h, dpi = 300) {
  ggsave(file.path(OUT_DIR, paste0(stub, ".png")), p, width = w, height = h, dpi = dpi)
  ggsave(file.path(OUT_DIR, paste0(stub, ".pdf")), p, width = w, height = h, dpi = dpi)
}

pick_first <- function(cands, nm) {
  hit <- intersect(cands, nm)
  if (length(hit) == 0) NULL else hit[1]
}

z_rowmean <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (length(vars) == 0) return(rep(NA_real_, nrow(df)))
  X <- scale(as.matrix(df[, vars, drop = FALSE]))
  out <- rowMeans(X, na.rm = TRUE)
  out[rowSums(!is.na(X)) == 0] <- NA_real_
  as.numeric(out)
}

# Plausible-value and jackknife engine for ICILS
pv_jrr <- function(data, pv_stem, model_fn, n_rep = N_REP, jrr_fac = JRR_FAC,
                   wt_base = "TOTWGTS", wt_prefix = "SRWGT") {
  full_ests <- lapply(pv_stem, function(pv) model_fn(data, pv, wt_base))
  param_names <- names(full_ests[[1]])
  M <- length(pv_stem)

  Q_bar <- setNames(
    sapply(param_names, function(p) mean(sapply(full_ests, `[`, p), na.rm = TRUE)),
    param_names)

  sampling_vars <- lapply(pv_stem, function(pv) {
    theta_full <- model_fn(data, pv, wt_base)
    rep_ests <- lapply(seq_len(n_rep), function(r) {
      wc <- paste0(wt_prefix, r)
      if (!wc %in% names(data)) return(theta_full * NA_real_)
      model_fn(data, pv, wc)
    })
    rep_mat <- do.call(rbind, rep_ests)
    sapply(param_names, function(p)
      jrr_fac * sum((rep_mat[, p] - theta_full[p])^2, na.rm = TRUE))
  })
  U_bar <- Reduce(`+`, sampling_vars) / M
  B_m   <- sapply(param_names, function(p) var(sapply(full_ests, `[`, p), na.rm = TRUE))
  se    <- sqrt(U_bar + (1 + 1/M) * B_m)

  tibble(parameter = param_names, estimate = Q_bar, se = se,
         t_value = Q_bar / se,
         p_value = 2 * pt(-abs(Q_bar / se), df = max(n_rep - 1, 1)),
         ci_lo = Q_bar - 1.96 * se, ci_hi = Q_bar + 1.96 * se)
}

# Plausible-value and jackknife engine for TIMSS
make_jk2_weights <- function(data) {
  zones <- sort(unique(data$JKZONE[!is.na(data$JKZONE)]))
  rw <- matrix(data$TOTWGT, nrow = nrow(data), ncol = length(zones))
  for (z in seq_along(zones)) {
    iz <- data$JKZONE == zones[z]; iz[is.na(iz)] <- FALSE
    r0 <- data$JKREP == 0 & iz; r1 <- data$JKREP == 1 & iz
    r0[is.na(r0)] <- FALSE; r1[is.na(r1)] <- FALSE
    rw[r0, z] <- 0; rw[r1, z] <- 2 * data$TOTWGT[r1]
  }
  as.data.frame(rw)
}

pv_jk2 <- function(data, pv_names, model_fn, coef_names = NULL) {
  pv_ests <- lapply(pv_names, function(pv) model_fn(data, pv, "TOTWGT"))
  if (is.null(coef_names)) coef_names <- names(pv_ests[[1]])
  M <- length(pv_names)
  theta_bar <- rowMeans(do.call(cbind, pv_ests)); names(theta_bar) <- coef_names
  B <- rowMeans(do.call(cbind, lapply(pv_ests, function(e) (e - theta_bar)^2)))
  rep_wts <- make_jk2_weights(data); n_reps <- ncol(rep_wts)
  samp_vars <- lapply(pv_names, function(pv) {
    tf <- model_fn(data, pv, "TOTWGT")
    re <- lapply(seq_len(n_reps), function(r) {
      data$REP_WT <- rep_wts[[r]]
      tryCatch(model_fn(data, pv, "REP_WT"),
               error = function(e) rep(NA_real_, length(coef_names)))
    })
    rowSums((do.call(cbind, re) - tf)^2, na.rm = TRUE)
  })
  U_bar <- rowMeans(do.call(cbind, samp_vars), na.rm = TRUE)
  se <- sqrt(U_bar + (1 + 1/M) * B)
  tibble(parameter = coef_names, estimate = theta_bar, se = se,
         t_val = theta_bar / se,
         p_val = 2 * pt(-abs(theta_bar / se), df = max(n_reps - 1, 1)))
}

# Forest-plot helper
make_forest <- function(df, meta_obj, title, xlab, scale_mult = 1, caption = "") {
  pooled <- tibble(country_label = "Pooled (RE)", estimate = as.numeric(meta_obj$beta),
    ci_lo = meta_obj$ci.lb, ci_hi = meta_obj$ci.ub, region = "Pooled", se = meta_obj$se)
  fd <- bind_rows(df, pooled) %>%
    mutate(est_s = estimate * scale_mult, lo_s = ci_lo * scale_mult, hi_s = ci_hi * scale_mult,
           is_pooled = country_label == "Pooled (RE)",
           lab = sprintf("%+.2f", est_s),
           country_label = fct_reorder(country_label, est_s),
           country_label = fct_relevel(country_label, "Pooled (RE)"))
  ggplot(fd, aes(x = est_s, y = country_label, colour = ifelse(is_pooled, "Pooled", region))) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey60", linetype = "dashed") +
    geom_errorbarh(aes(xmin = lo_s, xmax = hi_s), height = 0.25, linewidth = 0.55) +
    geom_point(aes(size = ifelse(is_pooled, 4, 2.6)), show.legend = FALSE) +
    geom_text(aes(label = lab), hjust = ifelse(fd$est_s >= 0, -0.15, 1.15),
              size = 3, colour = "grey20", show.legend = FALSE) +
    scale_colour_manual(values = pal_region, name = "Region") +
    scale_size_identity() +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.16))) +
    labs(title = title, x = xlab, y = NULL, caption = caption)
}

cat("Setup complete.\n")
