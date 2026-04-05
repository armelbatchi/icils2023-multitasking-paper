# ICILS secondary analyses

if (!output_exists("table_moderation_cil.csv") && "S_LRNSAFE" %in% names(bsg)) {
  mod_fn <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$S_ACMULT) & !is.na(d$S_LRNSAFE) & !is.na(d[[wc]])
    for (cv in available_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 500) return(c(`S_ACMULT` = NA, `S_LRNSAFE` = NA, `S_ACMULT:S_LRNSAFE` = NA))
    rhs <- paste(c("S_ACMULT * S_LRNSAFE", available_covs, "factor(CNTRYID)"), collapse = " + ")
    m <- lm(as.formula(paste(oc, "~", rhs)), data = d, weights = d[[wc]])
    cc <- coef(m)
    needed <- c("S_ACMULT", "S_LRNSAFE", "S_ACMULT:S_LRNSAFE")
    cc[intersect(needed, names(cc))]
  }
  mod_res <- pv_jrr(bsg, PV_CIL, mod_fn)
  write_csv(mod_res, file.path(OUT_DIR, "table_moderation_cil.csv"))
  cat("Moderation results:\n"); print(mod_res)

  # Moderation figure
  int_data <- bsg %>%
    filter(!is.na(PV1CIL), !is.na(S_ACMULT), !is.na(S_LRNSAFE), !is.na(TOTWGTS)) %>%
    mutate(lrnsafe_grp = case_when(
      S_LRNSAFE <= quantile(S_LRNSAFE, .25, na.rm = TRUE) ~ "Low (bottom 25%)",
      S_LRNSAFE >= quantile(S_LRNSAFE, .75, na.rm = TRUE) ~ "High (top 25%)",
      TRUE ~ "Middle 50%") %>%
        factor(levels = c("Low (bottom 25%)", "Middle 50%", "High (top 25%)")))

  fig6 <- ggplot(int_data, aes(S_ACMULT, PV1CIL, colour = lrnsafe_grp, fill = lrnsafe_grp)) +
    geom_smooth(method = "lm", formula = y ~ x, alpha = .12, linewidth = 1) +
    scale_colour_manual(values = c("#E76F51", "#E9C46A", "#2A9D8F"),
                        name = "Safe-ICT instruction") +
    scale_fill_manual(values = c("#E76F51", "#E9C46A", "#2A9D8F"),
                      name = "Safe-ICT instruction") +
    coord_cartesian(xlim = c(30, 70)) +
    labs(title = "Safe-ICT instruction moderates the multitasking-CIL link",
         subtitle = "Simple slopes at low, middle, and high levels of S_LRNSAFE.",
         x = "S_ACMULT", y = "CIL (PV1)",
         caption = "Source: ICILS 2023. Pooled sample, weighted OLS.")
  save_both(fig6, "fig6_moderation", 9, 6.5)
} else {
  cat("Moderation results already exist or S_LRNSAFE not available.\n")
}

if (HAVE_LME4 && !output_exists("multilevel_summaries.txt")) {
  ml_data <- bsg %>%
    filter(!is.na(PV1CIL), !is.na(S_ACMULT), !is.na(TOTWGTS), !is.na(IDSCHOOL)) %>%
    mutate(school_id = paste(CNTRYID, IDSCHOOL, sep = "_"),
           wt_norm = TOTWGTS / mean(TOTWGTS))
  for (cv in available_covs) ml_data <- ml_data %>% filter(!is.na(.data[[cv]]))

  ml_rhs <- paste(c("S_ACMULT", available_covs), collapse = " + ")

  cat("Fitting M1 (random intercept)...\n")
  m1 <- lmer(as.formula(paste("PV1CIL ~", ml_rhs,
                               "+ factor(CNTRYID) + (1|school_id)")),
             data = ml_data, weights = wt_norm, REML = TRUE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 30000)))

  cat("Fitting M2 (random slope)...\n")
  m2 <- lmer(as.formula(paste("PV1CIL ~", ml_rhs,
                               "+ factor(CNTRYID) + (1 + S_ACMULT|school_id)")),
             data = ml_data, weights = wt_norm, REML = TRUE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))

  sink(file.path(OUT_DIR, "multilevel_summaries.txt"))
  cat("=== M1 (random intercept) ===\n"); print(summary(m1), correlation = FALSE)
  cat("\n=== M2 (random slope) ===\n"); print(summary(m2), correlation = FALSE)
  cat("\n=== LR test ===\n"); print(anova(m1, m2))
  sink()

  # Caterpillar plot
  re <- ranef(m2)$school_id %>%
    as_tibble(rownames = "sid") %>%
    rename(int = `(Intercept)`, slope = S_ACMULT) %>%
    mutate(CNTRYID = str_extract(sid, "^[^_]+")) %>%
    left_join(region_map, by = "CNTRYID") %>%
    arrange(slope) %>%
    mutate(rank = row_number(), total = slope + fixef(m2)["S_ACMULT"])

  fig8 <- ggplot(re, aes(rank, total)) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_hline(yintercept = fixef(m2)["S_ACMULT"], colour = "#E63946", linewidth = .6) +
    geom_point(aes(colour = region), alpha = .4, size = .8) +
    scale_colour_manual(values = pal_region[names(pal_region) != "Pooled"], name = "Region") +
    annotate("text", x = max(re$rank) * .95, y = fixef(m2)["S_ACMULT"] + .15,
             label = sprintf("Fixed = %.2f", fixef(m2)["S_ACMULT"]),
             colour = "#E63946", size = 3.2, hjust = 1, fontface = "italic") +
    labs(title = "School-level variation in the multitasking slope",
         x = "School rank", y = "School-specific beta for S_ACMULT on CIL",
         caption = "Source: ICILS 2023. Two-level model (PV1CIL).")
  save_both(fig8, "fig8_caterpillar", 10, 6)
} else {
  cat("Multilevel summaries already exist or lme4 unavailable.\n")
}

if (!output_exists("quantile_results.csv")) {
  qd <- bsg %>%
    filter(!is.na(PV1CIL), !is.na(S_ACMULT), !is.na(TOTWGTS)) %>%
    mutate(wt_norm = TOTWGTS / mean(TOTWGTS))
  for (cv in available_covs) qd <- qd %>% filter(!is.na(.data[[cv]]))

  taus <- c(.10, .25, .50, .75, .90)
  qr_covs <- available_covs[available_covs %in% names(qd)]

  qr_res <- if (HAVE_QUANTREG && nrow(qd) > 1000) {
    lapply(taus, function(tau) {
      cat(sprintf("  tau = %.2f ...\n", tau))
      m <- rq(as.formula(paste("PV1CIL ~ S_ACMULT +",
                                paste(qr_covs, collapse = " + "))),
              data = qd, tau = tau, weights = wt_norm)
      s <- summary(m, se = "boot", R = 200)
      cf <- coef(s)
      tibble(tau = tau,
             estimate = cf["S_ACMULT", "Value"],
             se = cf["S_ACMULT", "Std. Error"])
    }) %>% bind_rows()
  } else {
    cat("Using binned-OLS fallback...\n")
    qd <- qd %>% mutate(q5 = ntile(PV1CIL, 5))
    lapply(1:5, function(q) {
      d <- qd %>% filter(q5 == q)
      m <- lm(as.formula(paste("PV1CIL ~ S_ACMULT +",
                                paste(qr_covs, collapse = " + "))),
              data = d, weights = wt_norm)
      cf <- coef(summary(m))
      tibble(tau = taus[q],
             estimate = cf["S_ACMULT", "Estimate"],
             se = cf["S_ACMULT", "Std. Error"])
    }) %>% bind_rows()
  }

  qr_res <- qr_res %>% mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)
  write_csv(qr_res, file.path(OUT_DIR, "quantile_results.csv"))
} else {
  qr_res <- read_csv(file.path(OUT_DIR, "quantile_results.csv"), show_col_types = FALSE)
  cat("Loaded cached quantile_results.csv\n")
}

fig11 <- ggplot(qr_res, aes(tau, estimate)) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "#457B9D", alpha = .18) +
  geom_line(colour = "#264653", linewidth = 1) +
  geom_point(colour = "#264653", size = 3) +
  scale_x_continuous(breaks = c(.1, .25, .5, .75, .9),
                     labels = paste0("Q", c(10, 25, 50, 75, 90))) +
  labs(title = "Multitasking-CIL association across the performance distribution",
       subtitle = "Coefficient of S_ACMULT at each quantile of CIL.",
       x = "CIL quantile", y = "S_ACMULT coefficient",
       caption = "Source: ICILS 2023 (PV1CIL).")
save_both(fig11, "fig11_quantile", 8, 6)

if (!output_exists("table_robustness_binary.csv")) {
  bsg$freq_sm <- as.integer(bsg$IS3G21C >= 3)

  bin_fn <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$freq_sm) & !is.na(d[[wc]])
    for (cv in available_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 200) return(c(freq_sm = NA_real_))
    m <- lm(as.formula(paste(oc, "~ freq_sm +",
                              paste(available_covs, collapse = " + "),
                              "+ factor(CNTRYID)")),
            data = d, weights = d[[wc]])
    coef(m)[intersect(c("freq_sm", "SES_CAN"), names(coef(m)))]
  }

  rob <- pv_jrr(bsg, PV_CIL, bin_fn)
  write_csv(rob, file.path(OUT_DIR, "table_robustness_binary.csv"))
  cat("Binary exposure robustness:\n"); print(rob)
} else {
  cat("Binary robustness already exists.\n")
}
