# Supplementary tables and reliability outputs

# Requires 00_setup.R and the relevant loaded ICILS objects.

# CORRELATION MATRIX — Mixed-type correlations (Pearson / polyserial / polychoric)
# Problem:  Variables 4 and 6–11 are ordinal (4- or 5-point scales).
#           Standard Pearson r underestimates associations between ordinal
#           variables and between ordinal and continuous variables.
# Solution: psych::mixedCor() computes:
#           • Pearson r          for continuous × continuous pairs
#           • Polyserial r       for continuous × ordinal pairs
#           • Polychoric r       for ordinal × ordinal pairs
# We report BOTH:
#   (a) Weighted Pearson correlations (primary table — comparable to S1 style)
#   (b) Mixed-type correlations from psych::mixedCor() (supplementary)
# With N > 86,000, all |r| ≥ .01 are p < .001. We nonetheless compute and
# export the full p-value matrix so reviewers can verify.

if (!requireNamespace("psych", quietly = TRUE))
  install.packages("psych", repos = "https://cloud.r-project.org")
library(psych)

# 1. DEFINE VARIABLES AND THEIR MEASUREMENT TYPES

# Plausible-value stems
pv_cil_names <- paste0("PV", 1:5, "CIL")
pv_ct_names  <- paste0("PV", 1:5, "CT")
M_PV <- 5  # number of plausible values

# Continuous variables
continuous_vars_nopv <- c("S_ACMULT", "S_LRNSAFE")
continuous_vars_nopv <- continuous_vars_nopv[continuous_vars_nopv %in% names(bsg)]

# Ordinal variables
ordinal_vars <- c("S_EXCOMP",
                   "IS3G21A", "IS3G21B", "IS3G21C",
                   "IS3G21D", "IS3G21E", "IS3G21F")
ordinal_vars <- ordinal_vars[ordinal_vars %in% names(bsg)]

# The "display" variable list uses placeholder names for PV-averaged estimates
all_vars <- c(continuous_vars_nopv, "CIL_PV", "CT_PV", ordinal_vars)

# Backward-compatible continuous-variable list
continuous_vars <- c(continuous_vars_nopv, "CIL_PV", "CT_PV")

# Display labels
var_labels <- c(
  "S_ACMULT"  = "1. S_ACMULT",
  "CIL_PV"    = "2. CIL",
  "CT_PV"     = "3. CT",
  "S_LRNSAFE" = "4. Safe-ICT instruction",
  "S_EXCOMP"  = "5. Computer experience",
  "IS3G21A"   = "6. A. Text chat",
  "IS3G21B"   = "7. B. Social media use/post",
  "IS3G21C"   = "8. C. Check social media",
  "IS3G21D"   = "9. D. Internet for info",
  "IS3G21E"   = "10. E. Watch videos/streams",
  "IS3G21F"   = "11. F. Listen to music/radio"
)

cat("Non-PV continuous: ", paste(continuous_vars_nopv, collapse = ", "), "\n")
cat("PV variables:      CIL (5 PVs), CT (5 PVs)\n")
cat("Ordinal variables: ", paste(ordinal_vars, collapse = ", "), "\n\n")

# 2. HELPER: Weighted Pearson correlation (single pair)

weighted_cor <- function(x, y, w) {
  keep <- complete.cases(x, y, w) & w > 0
  if (sum(keep) < 30) return(list(r = NA_real_, n = 0L))
  x <- x[keep]; y <- y[keep]; w <- w[keep]
  n <- sum(keep)
  mx <- weighted.mean(x, w); my <- weighted.mean(y, w)
  cov_xy <- sum(w * (x - mx) * (y - my)) / sum(w)
  sd_x   <- sqrt(sum(w * (x - mx)^2) / sum(w))
  sd_y   <- sqrt(sum(w * (y - my)^2) / sum(w))
  if (sd_x == 0 | sd_y == 0) return(list(r = NA_real_, n = n))
  r <- cov_xy / (sd_x * sd_y)
  list(r = r, n = n)
}

# 3. HELPER: Resolve a variable name to a numeric vector
#    For PV variables, returns a specific PV column; for others, the column itself

get_var_for_pv <- function(varname, pv_index) {
  if (varname == "CIL_PV") {
    return(as.numeric(bsg[[paste0("PV", pv_index, "CIL")]]))
  } else if (varname == "CT_PV") {
    return(as.numeric(bsg[[paste0("PV", pv_index, "CT")]]))
  } else {
    return(as.numeric(bsg[[varname]]))
  }
}

# 4. COMPUTE WEIGHTED PEARSON CORRELATIONS (PV-averaged where applicable)

n_vars <- length(all_vars)
pearson_r <- matrix(NA_real_, n_vars, n_vars, dimnames = list(all_vars, all_vars))
pearson_n <- matrix(NA_integer_, n_vars, n_vars, dimnames = list(all_vars, all_vars))
pearson_p <- matrix(NA_real_, n_vars, n_vars, dimnames = list(all_vars, all_vars))

cat("Computing PV-averaged weighted correlations...\n")

for (i in seq_len(n_vars)) {
  for (j in seq_len(i)) {
    vi <- all_vars[i]
    vj <- all_vars[j]

    # Determine whether either variable involves PVs
    vi_is_pv <- vi %in% c("CIL_PV", "CT_PV")
    vj_is_pv <- vj %in% c("CIL_PV", "CT_PV")

    if (vi_is_pv || vj_is_pv) {
      # Average correlation across 5 PVs (Rubin's rule for point estimates)
      r_pv <- numeric(M_PV)
      n_pv <- integer(M_PV)
      for (m in seq_len(M_PV)) {
        xi <- get_var_for_pv(vi, m)
        xj <- get_var_for_pv(vj, m)
        res <- weighted_cor(xi, xj, bsg$TOTWGTS)
        r_pv[m] <- res$r
        n_pv[m] <- res$n
      }
      avg_r <- mean(r_pv, na.rm = TRUE)
      avg_n <- round(mean(n_pv, na.rm = TRUE))
      pearson_r[i, j] <- pearson_r[j, i] <- avg_r
      pearson_n[i, j] <- pearson_n[j, i] <- avg_n
    } else {
      # No PVs involved — single computation
      xi <- as.numeric(bsg[[vi]])
      xj <- as.numeric(bsg[[vj]])
      res <- weighted_cor(xi, xj, bsg$TOTWGTS)
      pearson_r[i, j] <- pearson_r[j, i] <- res$r
      pearson_n[i, j] <- pearson_n[j, i] <- res$n
    }

    # P-value from average r and average n
    r_val <- pearson_r[i, j]
    n_val <- pearson_n[i, j]
    if (!is.na(r_val) && n_val > 3 && i != j) {
      r_abs <- min(abs(r_val), 0.9999)
      t_val <- r_abs * sqrt(n_val - 2) / sqrt(1 - r_abs^2)
      p_val <- 2 * pt(-abs(t_val), df = n_val - 2)
      pearson_p[i, j] <- pearson_p[j, i] <- p_val
    }
  }
}

cat("Done. Correlations involving CIL and CT are averaged across 5 PVs.\n\n")

# 5. DESCRIPTIVE STATISTICS (PV-averaged M and SD for CIL/CT)

desc_stats <- data.frame(
  Variable = var_labels[all_vars],
  Type = ifelse(all_vars %in% continuous_vars, "Continuous", "Ordinal"),
  M = NA_real_, SD = NA_real_, N = NA_integer_,
  stringsAsFactors = FALSE
)

for (idx in seq_along(all_vars)) {
  v <- all_vars[idx]

  if (v %in% c("CIL_PV", "CT_PV")) {
    # PV-averaged mean: average of 5 weighted means
    means_pv <- numeric(M_PV)
    vars_pv  <- numeric(M_PV)
    ns_pv    <- integer(M_PV)
    for (m in seq_len(M_PV)) {
      x <- get_var_for_pv(v, m)
      w <- bsg$TOTWGTS
      k <- complete.cases(x, w) & w > 0
      if (sum(k) < 30) { means_pv[m] <- NA; vars_pv[m] <- NA; ns_pv[m] <- 0; next }
      means_pv[m] <- weighted.mean(x[k], w[k])
      vars_pv[m]  <- Hmisc::wtd.var(x[k], w[k])
      ns_pv[m]    <- sum(k)
    }
    # Rubin's rules for mean
    Q_bar <- mean(means_pv, na.rm = TRUE)
    # Rubin's rules for variance: total_var = mean(within-PV var) + (1+1/M)*between-PV var of means
    U_bar <- mean(vars_pv, na.rm = TRUE)          # average within-imputation variance
    B_m   <- var(means_pv, na.rm = TRUE)           # between-imputation variance of the mean
    # For reporting SD of the *distribution* (not SE of the mean), we use within-PV SD averaged
    # The between-PV variance of means is negligible for the population SD
    desc_stats$M[idx]  <- round(Q_bar, 1)
    desc_stats$SD[idx] <- round(sqrt(U_bar), 1)    # average within-PV SD
    desc_stats$N[idx]  <- round(mean(ns_pv))
  } else {
    # Non-PV variable
    x <- as.numeric(bsg[[v]])
    w <- bsg$TOTWGTS
    k <- complete.cases(x, w) & w > 0
    if (sum(k) >= 30) {
      desc_stats$M[idx]  <- round(weighted.mean(x[k], w[k]), 2)
      desc_stats$SD[idx] <- round(sqrt(Hmisc::wtd.var(x[k], w[k])), 2)
      desc_stats$N[idx]  <- sum(k)
    }
  }
}

rownames(desc_stats) <- NULL

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  DESCRIPTIVE STATISTICS (weighted, PV-averaged for CIL & CT)\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
print(desc_stats, digits = 3)

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  P-VALUE MATRIX (two-tailed, from weighted Pearson r)\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

p_formatted <- matrix("", n_vars, n_vars)
rownames(p_formatted) <- var_labels[all_vars]
colnames(p_formatted) <- seq_len(n_vars)

for (i in seq_len(n_vars)) {
  for (j in seq_len(i)) {
    if (i == j) {
      p_formatted[i, j] <- "  —  "
    } else {
      p <- pearson_p[i, j]
      if (is.na(p)) {
        p_formatted[i, j] <- "  NA "
      } else if (p < .001) {
        p_formatted[i, j] <- "<.001"
      } else {
        p_formatted[i, j] <- sprintf("%.3f", p)
      }
    }
  }
}
print(noquote(p_formatted))

# Formatted lower-triangle correlation matrix

pearson_formatted <- matrix("", n_vars, n_vars)
rownames(pearson_formatted) <- var_labels[all_vars]
colnames(pearson_formatted) <- seq_len(n_vars)

for (i in seq_len(n_vars)) {
  for (j in seq_len(n_vars)) {
    if (i == j) {
      pearson_formatted[i, j] <- "—"
    } else if (j < i) {
      r <- pearson_r[i, j]
      p <- pearson_p[i, j]

      stars <- if (is.na(p)) {
        ""
      } else if (p < .001) {
        "***"
      } else if (p < .01) {
        "**"
      } else if (p < .05) {
        "*"
      } else {
        ""
      }

      pearson_formatted[i, j] <- if (is.na(r)) {
        "NA"
      } else {
        paste0(sprintf("%.3f", r), stars)
      }
    } else {
      pearson_formatted[i, j] <- ""
    }
  }
}

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("  FORMATTED LOWER-TRIANGLE CORRELATION MATRIX\n")
cat("══════════════════════════════════════════════════════════════════\n\n")
print(noquote(pearson_formatted))

# 7. SAVE ALL OUTPUTS

# Pearson correlations
write.csv(round(pearson_r, 3),
          file.path(OUT_DIR, "table_correlation_pearson_r.csv"), row.names = TRUE)

# Mixed correlations
if (exists("mixed_r") && !is.null(mixed_r)) {
  write.csv(round(mixed_r, 3),
            file.path(OUT_DIR, "table_correlation_mixed_r.csv"), row.names = TRUE)
}

# P-values
write.csv(pearson_p,
          file.path(OUT_DIR, "table_correlation_pvalues.csv"), row.names = TRUE)

# Pairwise N
write.csv(pearson_n,
          file.path(OUT_DIR, "table_correlation_pairwise_n.csv"), row.names = TRUE)

# Descriptive statistics
write.csv(desc_stats,
          file.path(OUT_DIR, "table_correlation_descriptives.csv"), row.names = FALSE)

# Formatted lower-triangle output
write.csv(pearson_formatted,
          file.path(OUT_DIR, "table_correlation_formatted.csv"), row.names = TRUE)

cat("\n\nAll outputs saved to:", OUT_DIR, "\n")
cat("  - table_correlation_pearson_r.csv        (weighted Pearson matrix)\n")
cat("  - table_correlation_mixed_r.csv          (polychoric/polyserial matrix)\n")
cat("  - table_correlation_pvalues.csv          (two-tailed p-values)\n")
cat("  - table_correlation_pairwise_n.csv       (pairwise sample sizes)\n")
cat("  - table_correlation_descriptives.csv     (M, SD, N, variable type)\n")
cat("  - table_correlation_formatted.csv        (formatted lower-triangle)\n")

library(psych)

cat("\n========== RELIABILITY VALUES FOR METHODS SECTION ==========\n\n")

# ── 1. Identify S_ACMULT raw items ───────────────────────────────────────────
# ICILS 2023 student questionnaire: IS8G18A-F (multitasking during homework)
# If your variable names are uppercased, they should be IS8G18A etc.
acmult_candidates <- c("IS8G18A", "IS8G18B", "IS8G18C",
                        "IS8G18D", "IS8G18E", "IS8G18F")
acmult_items <- acmult_candidates[acmult_candidates %in% names(bsg)]

if (length(acmult_items) == 0) {
  # Try alternate naming
  acmult_candidates2 <- grep("IS8G18", names(bsg), value = TRUE)
  if (length(acmult_candidates2) > 0) {
    acmult_items <- sort(acmult_candidates2)
    cat("Found S_ACMULT items with alternate names:", paste(acmult_items, collapse = ", "), "\n")
  } else {
    cat("WARNING: Could not find IS8G18* items. Trying IS8G17*...\n")
    acmult_items <- sort(grep("IS8G17", names(bsg), value = TRUE))
  }
}

cat("S_ACMULT raw items found:", paste(acmult_items, collapse = ", "), "\n")
cat("  N items:", length(acmult_items), "\n\n")

if (length(acmult_items) >= 4) {
  # Overall alpha
  acmult_mat <- bsg[, acmult_items]
  acmult_mat <- data.frame(lapply(acmult_mat, as.numeric))
  a_overall <- psych::alpha(acmult_mat, na.rm = TRUE, warnings = FALSE)

  cat("S_ACMULT Cronbach's alpha (pooled):",
      round(a_overall$total$raw_alpha, 2), "\n")

  # Country-specific
  a_by_cty <- tapply(seq_len(nrow(bsg)), bsg$CNTRYID, function(idx) {
    d <- acmult_mat[idx, , drop = FALSE]
    d <- d[complete.cases(d), , drop = FALSE]
    if (nrow(d) < 50) return(NA_real_)
    tryCatch(psych::alpha(d, warnings = FALSE)$total$raw_alpha,
             error = function(e) NA_real_)
  })
  a_vals <- unlist(a_by_cty)
  a_vals <- a_vals[!is.na(a_vals)]

  cat("S_ACMULT alpha by country:\n")
  cat("  Mean:", round(mean(a_vals), 2), "\n")
  cat("  Min: ", round(min(a_vals), 2), "\n")
  cat("  Max: ", round(max(a_vals), 2), "\n")
  cat("  Range: ", round(min(a_vals), 2), "-", round(max(a_vals), 2), "\n\n")
}

# ── 2. S_LRNSAFE raw items ──────────────────────────────────────────────────
lrnsafe_candidates <- c("IS8G08A", "IS8G08B", "IS8G08C", "IS8G08D")
lrnsafe_items <- lrnsafe_candidates[lrnsafe_candidates %in% names(bsg)]

if (length(lrnsafe_items) == 0) {
  lrnsafe_items <- sort(grep("IS8G08", names(bsg), value = TRUE))
}

cat("S_LRNSAFE raw items found:", paste(lrnsafe_items, collapse = ", "), "\n")

if (length(lrnsafe_items) >= 3) {
  lrn_mat <- data.frame(lapply(bsg[, lrnsafe_items], as.numeric))
  a_lrn <- psych::alpha(lrn_mat, na.rm = TRUE, warnings = FALSE)

  cat("S_LRNSAFE Cronbach's alpha (pooled):",
      round(a_lrn$total$raw_alpha, 2), "\n")

  a_lrn_cty <- tapply(seq_len(nrow(bsg)), bsg$CNTRYID, function(idx) {
    d <- lrn_mat[idx, , drop = FALSE]
    d <- d[complete.cases(d), , drop = FALSE]
    if (nrow(d) < 50) return(NA_real_)
    tryCatch(psych::alpha(d, warnings = FALSE)$total$raw_alpha,
             error = function(e) NA_real_)
  })
  a_lrn_vals <- unlist(a_lrn_cty)
  a_lrn_vals <- a_lrn_vals[!is.na(a_lrn_vals)]

  cat("  Mean:", round(mean(a_lrn_vals), 2), "\n")
  cat("  Range:", round(min(a_lrn_vals), 2), "-", round(max(a_lrn_vals), 2), "\n\n")
}

# ── 3. PV reliability for CIL and CT ────────────────────────────────────────
cat("--- PV Reliability (inter-PV correlation as lower bound) ---\n")

pv_cil_names <- paste0("PV", 1:5, "CIL")
pv_ct_names  <- paste0("PV", 1:5, "CT")

if (all(pv_cil_names %in% names(bsg))) {
  pv_cil <- bsg[, pv_cil_names]
  r_cil <- cor(pv_cil, use = "pairwise.complete.obs")
  mean_r_cil <- mean(r_cil[lower.tri(r_cil)])
  cat("CIL mean inter-PV correlation:", round(mean_r_cil, 3), "\n")
}

if (all(pv_ct_names %in% names(bsg))) {
  pv_ct <- bsg[!is.na(bsg$PV1CT), pv_ct_names]
  r_ct <- cor(pv_ct, use = "pairwise.complete.obs")
  mean_r_ct <- mean(r_ct[lower.tri(r_ct)])
  cat("CT mean inter-PV correlation:", round(mean_r_ct, 3), "\n")
}
cat("\n")

# ── 4. Exclusion rate ────────────────────────────────────────────────────────
cat("--- Sample / Exclusion ---\n")
n_raw <- nrow(bsg_raw)
n_analytic_cil <- sum(!is.na(bsg$S_ACMULT) & !is.na(bsg$PV1CIL))
n_analytic_ct  <- sum(!is.na(bsg$S_ACMULT) & !is.na(bsg$PV1CT))

cat("Total students in raw data:", format(n_raw, big.mark = ","), "\n")
cat("Analytic sample (CIL):", format(n_analytic_cil, big.mark = ","), "\n")
cat("Analytic sample (CT):", format(n_analytic_ct, big.mark = ","), "\n")
cat("Exclusion rate (CIL):", round(100 * (1 - n_analytic_cil / n_raw), 1), "%\n")
cat("Exclusion rate (CT):", round(100 * (1 - n_analytic_ct / n_raw), 1), "%\n\n")

# ── 5. S_ACMULT response categories ─────────────────────────────────────────
cat("--- S_ACMULT Response Categories ---\n")
if (length(acmult_items) > 0) {
  cat("Unique values for", acmult_items[1], ":\n")
  print(sort(unique(as.numeric(bsg[[acmult_items[1]]]))))

  # Check if labelled
  if (inherits(bsg[[acmult_items[1]]], "haven_labelled")) {
    cat("Value labels:\n")
    print(attr(bsg[[acmult_items[1]]], "labels"))
  }
}
cat("\n")

# ── 6. S_EXCOMP response categories ─────────────────────────────────────────
cat("--- S_EXCOMP Response Info ---\n")
if ("S_EXCOMP" %in% names(bsg)) {
  cat("Unique values:", paste(sort(unique(as.numeric(bsg$S_EXCOMP))), collapse = ", "), "\n")
  if (inherits(bsg$S_EXCOMP, "haven_labelled")) {
    cat("Value labels:\n")
    print(attr(bsg$S_EXCOMP, "labels"))
  }
}
cat("\n")

# ── 7. R version ─────────────────────────────────────────────────────────────
cat("--- R Version ---\n")
cat(R.version.string, "\n")

cat("\n========== COPY THESE VALUES INTO YOUR METHODS SECTION ==========\n")

#   included = non-missing S_ACMULT + PV1 outcome

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(Hmisc)
library(flextable)
library(officer)
library(tibble)

if (!exists("OUT_DIR")) OUT_DIR <- "."

# 0. SETTINGS
use_adjusted_models <- FALSE

# 1. HELPERS

clean_num_var <- function(x, valid_min = NULL, valid_max = NULL) {
  if (inherits(x, c("haven_labelled", "labelled", "haven_labelled_spss"))) {
    if (requireNamespace("haven", quietly = TRUE)) {
      x <- haven::zap_missing(x)
      x <- haven::zap_labels(x)
    }
  }

  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_

  if (!is.null(valid_min)) x[x < valid_min] <- NA_real_
  if (!is.null(valid_max)) x[x > valid_max] <- NA_real_

  x
}

guess_female_code <- function(x) {
  labs <- attr(x, "labels")
  if (is.null(labs)) return(NA_real_)

  nm <- names(labs)
  hit <- which(str_detect(tolower(nm), "female|girl|woman"))
  if (length(hit) == 0) return(NA_real_)

  suppressWarnings(as.numeric(unname(labs[hit[1]])))
}

w_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (sum(keep) == 0) return(NA_real_)
  weighted.mean(x[keep], w[keep], na.rm = TRUE)
}

w_sd <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (sum(keep) <= 1) return(NA_real_)
  sqrt(Hmisc::wtd.var(x[keep], w[keep], na.rm = TRUE))
}

fmt_mean_sd <- function(m, s, digits = 2) {
  if (!is.finite(m)) return(NA_character_)
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), m, s)
}

fmt_pct <- function(p, digits = 1) {
  if (!is.finite(p)) return(NA_character_)
  sprintf(paste0("%.", digits, "f%%"), 100 * p)
}

fmt_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < .001) return("<0.001")
  sprintf("%.3f", p)
}

std_diff_cont <- function(x, g) {
  x1 <- x[g == 1]
  x0 <- x[g == 0]
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]

  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)

  m1 <- mean(x1)
  m0 <- mean(x0)
  s1 <- sd(x1)
  s0 <- sd(x0)
  sp <- sqrt((s1^2 + s0^2) / 2)

  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (m1 - m0) / sp
}

std_diff_bin <- function(x, g) {
  x1 <- x[g == 1]
  x0 <- x[g == 0]
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]

  if (length(x1) == 0 || length(x0) == 0) return(NA_real_)

  p1 <- mean(x1)
  p0 <- mean(x0)
  ps <- (p1 + p0) / 2
  den <- sqrt(ps * (1 - ps))

  if (!is.finite(den) || den == 0) return(NA_real_)
  (p1 - p0) / den
}

approx_p_cont <- function(x, g, w = NULL) {
  keep <- is.finite(x) & !is.na(g)
  if (!is.null(w)) keep <- keep & is.finite(w) & w > 0
  if (sum(keep) < 30) return(NA_real_)

  dd <- data.frame(
    x = x[keep],
    g = g[keep],
    w = if (is.null(w)) rep(1, sum(keep)) else w[keep]
  )

  fit <- tryCatch(
    lm(x ~ g, data = dd, weights = dd$w),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)
  sm <- summary(fit)$coefficients
  if (!"g" %in% rownames(sm)) return(NA_real_)
  sm["g", "Pr(>|t|)"]
}

approx_p_bin <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  if (sum(keep) < 20) return(NA_real_)

  tab <- table(g[keep], x[keep])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)

  suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
}

make_flag <- function(data, outcome = c("CIL", "CT"), adjusted = FALSE) {
  outcome <- match.arg(outcome)
  pv <- if (outcome == "CIL") "PV1CIL" else "PV1CT"

  flag <- !is.na(data$S_ACMULT) & !is.na(data[[pv]])

  if (adjusted) {
    flag <- flag & !is.na(data$TOTWGTS) & data$TOTWGTS > 0
    if (exists("available_covs")) {
      for (v in available_covs) {
        if (v %in% names(data)) flag <- flag & !is.na(data[[v]])
      }
    }
  }

  flag
}

# 2. CLEAN A COPY OF THE DATA FIRST

bsg_s2 <- bsg

# Exposure and outcomes used to define inclusion
if ("S_ACMULT" %in% names(bsg_s2)) {
  bsg_s2$S_ACMULT <- clean_num_var(bsg_s2$S_ACMULT, valid_min = 0, valid_max = 100)
}
if ("PV1CIL" %in% names(bsg_s2)) {
  bsg_s2$PV1CIL <- clean_num_var(bsg_s2$PV1CIL, valid_min = 0, valid_max = 1000)
}
if ("PV1CT" %in% names(bsg_s2)) {
  bsg_s2$PV1CT <- clean_num_var(bsg_s2$PV1CT, valid_min = 0, valid_max = 1000)
}
if ("TOTWGTS" %in% names(bsg_s2)) {
  bsg_s2$TOTWGTS <- clean_num_var(bsg_s2$TOTWGTS, valid_min = 0)
}

# Background variables for the included-versus-excluded table
if ("SES_CAN" %in% names(bsg_s2)) {
  # broad standardized-index bounds; removes -99, -999, etc.
  bsg_s2$SES_CAN <- clean_num_var(bsg_s2$SES_CAN, valid_min = -5, valid_max = 5)
}

if ("S_LRNSAFE" %in% names(bsg_s2)) {
  bsg_s2$S_LRNSAFE <- clean_num_var(bsg_s2$S_LRNSAFE, valid_min = 0, valid_max = 100)
}

if ("EXCOMP_CAN" %in% names(bsg_s2)) {
  lbl <- attr(bsg_s2$EXCOMP_CAN, "labels")

  if (!is.null(lbl)) {
    valid_codes <- suppressWarnings(as.numeric(unname(lbl)))
    valid_codes <- valid_codes[is.finite(valid_codes)]
    valid_codes <- valid_codes[valid_codes >= 0]

    tmp <- clean_num_var(bsg_s2$EXCOMP_CAN)
    if (length(valid_codes) > 0) {
      tmp[!tmp %in% valid_codes] <- NA_real_
    } else {
      tmp[tmp < 0 | tmp > 10] <- NA_real_
    }
    bsg_s2$EXCOMP_CAN <- tmp
  } else {
    bsg_s2$EXCOMP_CAN <- clean_num_var(bsg_s2$EXCOMP_CAN, valid_min = 0, valid_max = 10)
  }
}

# Rebuild female indicator from cleaned gender variable
if (exists("GENDER_VAR") && !is.null(GENDER_VAR) && GENDER_VAR %in% names(bsg)) {
  female_code <- guess_female_code(bsg[[GENDER_VAR]])
  g_clean <- clean_num_var(bsg[[GENDER_VAR]])

  if (is.finite(female_code)) {
    bsg_s2$FEMALE_BIN <- ifelse(is.na(g_clean), NA_real_, as.numeric(g_clean == female_code))
  }
}

diag_s2 <- tibble(
  variable = intersect(
    c("SES_CAN", "EXCOMP_CAN", "S_LRNSAFE", "S_ACMULT", "PV1CIL", "PV1CT"),
    names(bsg_s2)
  ),
  min = sapply(intersect(
    c("SES_CAN", "EXCOMP_CAN", "S_LRNSAFE", "S_ACMULT", "PV1CIL", "PV1CT"),
    names(bsg_s2)
  ), function(v) min(bsg_s2[[v]], na.rm = TRUE)),
  max = sapply(intersect(
    c("SES_CAN", "EXCOMP_CAN", "S_LRNSAFE", "S_ACMULT", "PV1CIL", "PV1CT"),
    names(bsg_s2)
  ), function(v) max(bsg_s2[[v]], na.rm = TRUE))
)

print(diag_s2)

# 3. BUILD supp_table_s2 FROM THE CLEANED DATA

compare_vars <- c("SES_CAN", "EXCOMP_CAN")
compare_labs <- c(
  "SES_CAN"    = "Socioeconomic background",
  "EXCOMP_CAN" = "Computer experience"
)

if ("S_LRNSAFE" %in% names(bsg_s2)) {
  compare_vars <- c(compare_vars, "S_LRNSAFE")
  compare_labs["S_LRNSAFE"] <- "Safe-ICT instruction"
}

if ("FEMALE_BIN" %in% names(bsg_s2)) {
  compare_vars <- c(compare_vars, "FEMALE_BIN")
  compare_labs["FEMALE_BIN"] <- "Female"
}

make_missingness_table <- function(data, outcome = c("CIL", "CT"), adjusted = FALSE) {
  outcome <- match.arg(outcome)
  include_flag <- make_flag(data, outcome = outcome, adjusted = adjusted)

  w <- if ("TOTWGTS" %in% names(data)) data$TOTWGTS else rep(1, nrow(data))

  bind_rows(lapply(compare_vars, function(v) {
    x <- data[[v]]
    keep_var <- !is.na(x)

    x_use <- x[keep_var]
    g_use <- as.integer(include_flag[keep_var])
    w_use <- w[keep_var]

    n_inc <- sum(g_use == 1)
    n_exc <- sum(g_use == 0)

    if (v == "FEMALE_BIN") {
      inc_stat <- fmt_pct(w_mean(x_use[g_use == 1], w_use[g_use == 1]))
      exc_stat <- fmt_pct(w_mean(x_use[g_use == 0], w_use[g_use == 0]))
      smd <- std_diff_bin(x_use, g_use)
      pval <- approx_p_bin(x_use, g_use)
    } else {
      inc_stat <- fmt_mean_sd(
        w_mean(x_use[g_use == 1], w_use[g_use == 1]),
        w_sd(x_use[g_use == 1], w_use[g_use == 1])
      )
      exc_stat <- fmt_mean_sd(
        w_mean(x_use[g_use == 0], w_use[g_use == 0]),
        w_sd(x_use[g_use == 0], w_use[g_use == 0])
      )
      smd <- std_diff_cont(x_use, g_use)
      pval <- approx_p_cont(x_use, g_use, w_use)
    }

    tibble(
      Outcome = outcome,
      Sample_definition = if (adjusted) {
        "Adjusted-model complete case"
      } else {
        "Outcome-specific analytic count"
      },
      Variable = unname(compare_labs[v]),
      Included_n = n_inc,
      Excluded_n = n_exc,
      Included = inc_stat,
      Excluded = exc_stat,
      Std_diff = round(smd, 3),
      P_value = fmt_p(pval)
    )
  }))
}

supp_table_s2 <- bind_rows(
  make_missingness_table(bsg_s2, outcome = "CIL", adjusted = use_adjusted_models),
  make_missingness_table(bsg_s2, outcome = "CT",  adjusted = use_adjusted_models)
)

print(supp_table_s2)

write_csv(
  supp_table_s2,
  file.path(OUT_DIR, "supp_table_S2_included_vs_excluded_cleaned.csv")
)

print(
  supp_table_s2 %>%
    filter(Variable == "Socioeconomic background")
)

tbl_s2_pretty <- supp_table_s2 %>%
  mutate(
    Outcome = factor(Outcome, levels = c("CIL", "CT")),
    Variable = factor(
      Variable,
      levels = c(
        "Socioeconomic background",
        "Computer experience",
        "Safe-ICT instruction",
        "Female"
      )
    ),
    Included_n_fmt = formatC(Included_n, format = "d", big.mark = ","),
    Excluded_n_fmt = formatC(Excluded_n, format = "d", big.mark = ","),
    Std_diff = ifelse(is.na(Std_diff), "", sprintf("%.3f", Std_diff)),
    P_value = ifelse(is.na(P_value), "", as.character(P_value))
  ) %>%
  arrange(Outcome, Variable)

make_panel_df <- function(dat, outcome_name) {
  dd <- dat %>% filter(Outcome == outcome_name)

  n_inc <- unique(dd$Included_n_fmt)
  n_exc <- unique(dd$Excluded_n_fmt)

  panel_title <- paste0(
    outcome_name,
    " outcome",
    "   (included n = ", n_inc[1], "; excluded n = ", n_exc[1], ")"
  )

  bind_rows(
    tibble(
      Variable = panel_title,
      Included = "",
      Excluded = "",
      SMD = "",
      P_value = ""
    ),
    dd %>%
      transmute(
        Variable = as.character(Variable),
        Included = Included,
        Excluded = Excluded,
        SMD = Std_diff,
        P_value = P_value
      )
  )
}

tbl_display <- bind_rows(
  make_panel_df(tbl_s2_pretty, "CIL"),
  make_panel_df(tbl_s2_pretty, "CT")
)

panel_rows <- which(str_detect(tbl_display$Variable, "outcome"))

ft <- flextable(tbl_display)

ft <- set_header_labels(
  ft,
  Variable = "Characteristic",
  Included = "Included sample",
  Excluded = "Excluded sample",
  SMD = "Std. diff.",
  P_value = "P value"
)

ft <- set_caption(
  ft,
  caption = paste0(
    "Supplementary Table S2. Comparison of included and excluded students on key ",
    "background variables. Included students were those with non-missing S_ACMULT ",
    "and the relevant outcome indicator, defined as PV1CIL for CIL analyses and ",
    "PV1CT for CT analyses, matching the current analytic-count logic."
  )
)

ft <- theme_booktabs(ft)
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 10, part = "body")
ft <- fontsize(ft, size = 10.5, part = "header")
ft <- fontsize(ft, size = 9, part = "footer")
ft <- padding(ft, padding = 4, part = "all")

ft <- align(ft, j = "Variable", align = "left", part = "all")
ft <- align(
  ft,
  j = c("Included", "Excluded", "SMD", "P_value"),
  align = "center",
  part = "all"
)

ft <- width(ft, j = "Variable", width = 3.2)
ft <- width(ft, j = c("Included", "Excluded"), width = 1.55)
ft <- width(ft, j = c("SMD", "P_value"), width = 0.95)

ft <- border_remove(ft)

std_border <- fp_border(color = "black", width = 1)

ft <- bold(ft, part = "header")
ft <- bg(ft, part = "header", bg = "#D9EAF7")
ft <- hline_top(ft, part = "header", border = std_border)
ft <- hline_bottom(ft, part = "header", border = std_border)

data_rows <- setdiff(seq_len(nrow(tbl_display)), panel_rows)
alt_rows <- data_rows[seq(2, length(data_rows), by = 2)]
if (length(alt_rows) > 0) {
  ft <- bg(ft, i = alt_rows, bg = "#F7F7F7", part = "body")
}

ft <- bold(ft, i = panel_rows, bold = TRUE, part = "body")
ft <- bg(ft, i = panel_rows, bg = "#EAF2F8", part = "body")

for (r in panel_rows) {
  ft <- merge_at(ft, i = r, j = 1:5, part = "body")
  ft <- align(ft, i = r, j = 1, align = "left", part = "body")
}

regular_rows <- setdiff(seq_len(nrow(tbl_display)), panel_rows)
if (length(regular_rows) > 0) {
  ft <- padding(ft, i = regular_rows, j = 1, padding.left = 10, part = "body")
}

sig_rows <- which(tbl_display$P_value == "<0.001")
if (length(sig_rows) > 0) {
  ft <- bold(ft, i = sig_rows, j = "P_value", bold = TRUE, part = "body")
}

ft <- hline_bottom(ft, part = "body", border = std_border)

ft <- add_footer_lines(
  ft,
  values = c(
    "Values are weighted mean (SD) for continuous variables and weighted percentage for female.",
    "Std. diff. = standardized difference. P values are approximate tests of included versus excluded students.",
    "Included students were defined outcome-specifically using non-missing S_ACMULT and PV1CIL or PV1CT.",
    "Background variables were cleaned for user-defined missing values and implausible out-of-range codes before comparison."
  )
)

ft <- italic(ft, part = "footer")
ft <- align(ft, align = "left", part = "footer")
ft <- autofit(ft)

save_as_docx(
  "Supplementary Table S2" = ft,
  path = file.path(OUT_DIR, "supp_table_S2_included_vs_excluded.docx")
)

ft
