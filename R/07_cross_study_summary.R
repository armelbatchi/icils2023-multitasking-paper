# Cross-study comparison tables and meta summaries

# Requires prior ICILS and TIMSS analysis scripts.

comparison <- tribble(
  ~Study, ~Prong, ~Predictor, ~Outcome,
  "ICILS 2023", "---", "S_ACMULT", "CIL",
  "ICILS 2023", "---", "S_ACMULT", "CT"
)

# ICILS values
comparison$Pooled_beta <- c(
  round(as.numeric(meta_cil$beta), 3),
  if (exists("meta_ct")) round(as.numeric(meta_ct$beta), 3) else NA)
comparison$CI_95 <- c(
  sprintf("[%.3f, %.3f]", meta_cil$ci.lb, meta_cil$ci.ub),
  if (exists("meta_ct")) sprintf("[%.3f, %.3f]", meta_ct$ci.lb, meta_ct$ci.ub) else NA)
comparison$I2 <- c(
  round(meta_cil$I2, 1),
  if (exists("meta_ct")) round(meta_ct$I2, 1) else NA)
comparison$k <- c(
  meta_cil$k,
  if (exists("meta_ct")) meta_ct$k else NA)

# TIMSS values
add_timss_row <- function(prong, pred, outcome, res) {
  if (is.null(res) || is.null(res$meta)) return(NULL)
  m <- res$meta
  tibble(Study = "TIMSS 2023", Prong = prong, Predictor = pred, Outcome = outcome,
         Pooled_beta = round(as.numeric(coef(m)), 3),
         CI_95 = sprintf("[%.3f, %.3f]", m$ci.lb, m$ci.ub),
         I2 = round(m$I2, 1),
         k = as.integer(m$k))
}

# Prong 3 values from CSV
p3m_res <- read_csv(file.path(OUT_DIR, "P3_interaction_math.csv"), show_col_types = FALSE)
p3s_res <- read_csv(file.path(OUT_DIR, "P3_interaction_sci.csv"), show_col_types = FALSE)

p3_int_math <- p3m_res %>% filter(parameter == "INET_SCHOOL_z:DISTRACT_CLIMATE_z")
p3_int_sci  <- p3s_res %>% filter(parameter == "INET_SCHOOL_z:DISTRACT_CLIMATE_z")

comparison <- bind_rows(comparison,
  add_timss_row("P1", "GAME_FREQ",   "Math",    r_p1m),
  add_timss_row("P1", "GAME_FREQ",   "Science", r_p1s),
  add_timss_row("P2", "INET_SCHOOL", "Math",    r_p2m),
  add_timss_row("P2", "INET_SCHOOL", "Science", r_p2s),
  tibble(Study = "TIMSS 2023", Prong = "P3", Predictor = "INET × DISTRACT",
         Outcome = "Math",
         Pooled_beta = round(p3_int_math$estimate, 3),
         CI_95 = sprintf("[%.3f, %.3f]",
                         p3_int_math$estimate - 1.96 * p3_int_math$se,
                         p3_int_math$estimate + 1.96 * p3_int_math$se),
         I2 = NA_real_, k = NA_integer_),
  tibble(Study = "TIMSS 2023", Prong = "P3", Predictor = "INET × DISTRACT",
         Outcome = "Science",
         Pooled_beta = round(p3_int_sci$estimate, 3),
         CI_95 = sprintf("[%.3f, %.3f]",
                         p3_int_sci$estimate - 1.96 * p3_int_sci$se,
                         p3_int_sci$estimate + 1.96 * p3_int_sci$se),
         I2 = NA_real_, k = NA_integer_)
)

write_csv(comparison, file.path(OUT_DIR, "table_cross_study_comparison.csv"))
cat("\n===== CROSS-STUDY COMPARISON =====\n")
print(comparison)

meta_summary <- tibble()

if (exists("meta_cil")) {
  meta_summary <- bind_rows(meta_summary, tibble(
    outcome = "CIL", pooled = as.numeric(meta_cil$beta), se = meta_cil$se,
    ci_lo = meta_cil$ci.lb, ci_hi = meta_cil$ci.ub,
    I2 = meta_cil$I2, tau2 = meta_cil$tau2, k = meta_cil$k))
}
if (exists("meta_ct")) {
  meta_summary <- bind_rows(meta_summary, tibble(
    outcome = "CT", pooled = as.numeric(meta_ct$beta), se = meta_ct$se,
    ci_lo = meta_ct$ci.lb, ci_hi = meta_ct$ci.ub,
    I2 = meta_ct$I2, tau2 = meta_ct$tau2, k = meta_ct$k))
}

write_csv(meta_summary, file.path(OUT_DIR, "meta_summaries.csv"))

cat("\n=============================================================\n")
cat("  ALL ANALYSES COMPLETE\n")
cat("  Output directory:", OUT_DIR, "\n")
cat("=============================================================\n")
cat("Files created:\n")
list.files(OUT_DIR) %>% walk(~ cat("  ", .x, "\n"))
cat("=============================================================\n")
