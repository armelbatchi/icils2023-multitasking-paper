# TIMSS analyses and figures

if (!output_exists("P1_game_math.csv")) {
  rhs_p1 <- paste(c("GAME_FREQ_z", timss_covs), collapse = " + ")
  cn_p1  <- c("(Intercept)", "GAME_FREQ_z", timss_covs)

  run_p1 <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$GAME_FREQ_z) & !is.na(d[[wc]])
    for (cv in timss_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 50) return(setNames(rep(NA_real_, length(cn_p1)), cn_p1))
    m <- tryCatch(lm(as.formula(paste(oc, "~", rhs_p1)),
                      data = d, weights = d[[wc]]),
                  error = function(e) NULL)
    if (is.null(m)) return(setNames(rep(NA_real_, length(cn_p1)), cn_p1))
    cc <- coef(m)
    out <- setNames(rep(NA_real_, length(cn_p1)), cn_p1)
    out[names(cc)] <- cc; out
  }

  p1m <- list(); p1s <- list()
  for (cnt in timss_countries) {
    cd <- tsg[tsg$CNTRYID == cnt, , drop = FALSE]
    if (sum(!is.na(cd$BSMMAT01) & !is.na(cd$GAME_FREQ_z)) >= 100) {
      cat(sprintf("  P1 Math: %s ...", cnt))
      r <- tryCatch(pv_jk2(cd, PV_MATH, run_p1, cn_p1), error = function(e) NULL)
      if (!is.null(r)) { p1m[[cnt]] <- r %>% mutate(CNTRYID = cnt); cat(" ok\n") }
      else cat(" err\n")
    }
    if (sum(!is.na(cd$BSSSCI01) & !is.na(cd$GAME_FREQ_z)) >= 100) {
      r <- tryCatch(pv_jk2(cd, PV_SCI, run_p1, cn_p1), error = function(e) NULL)
      if (!is.null(r)) p1s[[cnt]] <- r %>% mutate(CNTRYID = cnt)
    }
  }
  write_csv(bind_rows(p1m), file.path(OUT_DIR, "P1_game_math.csv"))
  write_csv(bind_rows(p1s), file.path(OUT_DIR, "P1_game_sci.csv"))
  cat("Prong 1 results saved.\n")
} else {
  cat("Prong 1 results already exist.\n")
}

p1_math_df <- read_csv(file.path(OUT_DIR, "P1_game_math.csv"), show_col_types = FALSE)
p1_sci_df  <- read_csv(file.path(OUT_DIR, "P1_game_sci.csv"), show_col_types = FALSE)

if (!output_exists("P2_inet_math.csv")) {
  rhs_p2 <- paste(c("INET_SCHOOL_z", timss_covs), collapse = " + ")
  cn_p2  <- c("(Intercept)", "INET_SCHOOL_z", timss_covs)

  run_p2 <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$INET_SCHOOL_z) & !is.na(d[[wc]])
    for (cv in timss_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 50) return(setNames(rep(NA_real_, length(cn_p2)), cn_p2))
    m <- tryCatch(lm(as.formula(paste(oc, "~", rhs_p2)),
                      data = d, weights = d[[wc]]),
                  error = function(e) NULL)
    if (is.null(m)) return(setNames(rep(NA_real_, length(cn_p2)), cn_p2))
    cc <- coef(m)
    out <- setNames(rep(NA_real_, length(cn_p2)), cn_p2)
    out[names(cc)] <- cc; out
  }

  p2m <- list(); p2s <- list()
  for (cnt in timss_countries) {
    cd <- tsg[tsg$CNTRYID == cnt, , drop = FALSE]
    if (sum(!is.na(cd$BSMMAT01) & !is.na(cd$INET_SCHOOL_z)) >= 100) {
      cat(sprintf("  P2 Math: %s ...", cnt))
      r <- tryCatch(pv_jk2(cd, PV_MATH, run_p2, cn_p2), error = function(e) NULL)
      if (!is.null(r)) { p2m[[cnt]] <- r %>% mutate(CNTRYID = cnt); cat(" ok\n") }
      else cat(" err\n")
    }
    if (sum(!is.na(cd$BSSSCI01) & !is.na(cd$INET_SCHOOL_z)) >= 100) {
      r <- tryCatch(pv_jk2(cd, PV_SCI, run_p2, cn_p2), error = function(e) NULL)
      if (!is.null(r)) p2s[[cnt]] <- r %>% mutate(CNTRYID = cnt)
    }
  }
  write_csv(bind_rows(p2m), file.path(OUT_DIR, "P2_inet_math.csv"))
  write_csv(bind_rows(p2s), file.path(OUT_DIR, "P2_inet_sci.csv"))
  cat("Prong 2 results saved.\n")
} else {
  cat("Prong 2 results already exist.\n")
}

p2_math_df <- read_csv(file.path(OUT_DIR, "P2_inet_math.csv"), show_col_types = FALSE)
p2_sci_df  <- read_csv(file.path(OUT_DIR, "P2_inet_sci.csv"), show_col_types = FALSE)

if (!output_exists("P3_interaction_math.csv") && "DISTRACT_CLIMATE_z" %in% names(tsg)) {
  int_cn <- c("INET_SCHOOL_z", "DISTRACT_CLIMATE_z",
              "INET_SCHOOL_z:DISTRACT_CLIMATE_z", "(Intercept)")

  run_int <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$INET_SCHOOL_z) &
            !is.na(d$DISTRACT_CLIMATE_z) & !is.na(d[[wc]])
    for (cv in timss_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 200) return(setNames(rep(NA_real_, 4), int_cn))
    rhs <- paste(c("INET_SCHOOL_z * DISTRACT_CLIMATE_z", timss_covs,
                    "factor(CNTRYID)"), collapse = " + ")
    m <- tryCatch(lm(as.formula(paste(oc, "~", rhs)),
                      data = d, weights = d[[wc]]),
                  error = function(e) NULL)
    if (is.null(m)) return(setNames(rep(NA_real_, 4), int_cn))
    cc <- coef(m)
    out <- setNames(rep(NA_real_, 4), int_cn)
    out[intersect(int_cn, names(cc))] <- cc[intersect(int_cn, names(cc))]
    out
  }

  pd <- tsg %>% filter(!is.na(INET_SCHOOL_z), !is.na(DISTRACT_CLIMATE_z))
  cat("Running Prong 3 interaction (Math)...\n")
  p3m <- tryCatch(pv_jk2(pd, PV_MATH, run_int, int_cn), error = function(e) NULL)
  cat("Running Prong 3 interaction (Science)...\n")
  p3s <- tryCatch(pv_jk2(pd, PV_SCI,  run_int, int_cn), error = function(e) NULL)

  if (!is.null(p3m)) { write_csv(p3m, file.path(OUT_DIR, "P3_interaction_math.csv")); print(p3m) }
  if (!is.null(p3s)) { write_csv(p3s, file.path(OUT_DIR, "P3_interaction_sci.csv")); print(p3s) }
} else {
  cat("Prong 3 results already exist or DISTRACT_CLIMATE_z not available.\n")
}

make_timss_panel_forest <- function(df_list, param_list, panel_labels,
                                     main_title, xlab, stub,
                                     width = 14, height = 16) {
  panels <- list()
  meta_results <- list()

  for (i in seq_along(df_list)) {
    fd <- df_list[[i]] %>%
      label_timss() %>%
      filter(parameter == param_list[i], !is.na(estimate), !is.na(se), se > 0) %>%
      mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)

    if (nrow(fd) < 3) next

    meta <- metafor::rma(yi = estimate, sei = se, data = fd, method = "REML")
    meta_results[[panel_labels[i]]] <- meta

    pooled <- tibble(
      country_label = "Pooled (RE)",
      estimate = as.numeric(meta$beta),
      ci_lo = meta$ci.lb, ci_hi = meta$ci.ub,
      se = meta$se, CNTRYID = "Pooled"
    )

    fd <- bind_rows(fd, pooled) %>%
      mutate(
        is_pooled = country_label == "Pooled (RE)",
        country_label = fct_reorder(country_label, estimate),
        country_label = fct_relevel(country_label, "Pooled (RE)"),
        panel = panel_labels[i]
      )

    panels[[i]] <- fd
  }

  all_data <- bind_rows(panels) %>%
    mutate(panel = factor(panel, levels = panel_labels))

  fig <- ggplot(all_data, aes(x = estimate, y = country_label)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey60", linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   height = 0.3, linewidth = 0.45, colour = "grey55") +
    geom_point(aes(colour = is_pooled, shape = is_pooled, size = is_pooled)) +
    scale_colour_manual(values = c("FALSE" = "#457B9D", "TRUE" = "#E63946"), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2.2, "TRUE" = 4.5), guide = "none") +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = main_title,
      subtitle = "Country-specific adjusted coefficients with 95% CIs. Red diamond = pooled random-effects estimate.",
      x = xlab,
      y = NULL,
      caption = "Source: TIMSS 2023. Adjusted for digital self-efficacy, gender, and home educational resources."
    ) +
    theme(
      strip.text = element_text(size = 11, face = "bold", hjust = 0),
      axis.text.y = element_text(size = 8),
      panel.spacing = unit(1.5, "lines")
    )

  save_both(fig, stub, width, height)

  for (nm in names(meta_results)) {
    m <- meta_results[[nm]]
    cat(sprintf("  %s: pooled=%.2f [%.2f, %.2f], I²=%.1f%%, k=%d\n",
                nm, coef(m), m$ci.lb, m$ci.ub, m$I2, m$k))
  }

  meta_results
}

# ── Reload cached results ──
p1_math_df <- read_csv(file.path(OUT_DIR, "P1_game_math.csv"), show_col_types = FALSE)
p1_sci_df  <- read_csv(file.path(OUT_DIR, "P1_game_sci.csv"), show_col_types = FALSE)
p2_math_df <- read_csv(file.path(OUT_DIR, "P2_inet_math.csv"), show_col_types = FALSE)
p2_sci_df  <- read_csv(file.path(OUT_DIR, "P2_inet_sci.csv"), show_col_types = FALSE)

# ── Prong 1 paneled figure ──
cat("\n===== PRONG 1: Learning games frequency =====\n")
meta_p1 <- make_timss_panel_forest(
  df_list      = list(p1_math_df, p1_sci_df),
  param_list   = c("GAME_FREQ_z", "GAME_FREQ_z"),
  panel_labels = c("(a) Mathematics", "(b) Science"),
  main_title   = "Prong 1: Learning games frequency (BSBG12F) and achievement",
  xlab         = "Achievement points per 1 SD increase in game frequency",
  stub         = "fig_prong1_panels",
  width = 14, height = 16
)

# ── Prong 2 paneled figure ──
cat("\n===== PRONG 2: Internet use for schoolwork =====\n")
meta_p2 <- make_timss_panel_forest(
  df_list      = list(p2_math_df, p2_sci_df),
  param_list   = c("INET_SCHOOL_z", "INET_SCHOOL_z"),
  panel_labels = c("(a) Mathematics", "(b) Science"),
  main_title   = "Prong 2: Internet use for schoolwork (BSBG12A\u2013F) and achievement",
  xlab         = "Achievement points per 1 SD increase in internet use",
  stub         = "fig_prong2_panels",
  width = 14, height = 16
)

# ── Prong 3 figure ──
cat("\n===== PRONG 3: Distraction climate interaction =====\n")

p3m <- read_csv(file.path(OUT_DIR, "P3_interaction_math.csv"), show_col_types = FALSE)
p3s <- read_csv(file.path(OUT_DIR, "P3_interaction_sci.csv"), show_col_types = FALSE)

p3_plot <- bind_rows(
  p3m %>% mutate(outcome = "Mathematics"),
  p3s %>% mutate(outcome = "Science")
) %>%
  filter(parameter != "(Intercept)") %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    parameter = case_when(
      parameter == "INET_SCHOOL_z" ~ "Internet use for schoolwork\n(student-level main effect)",
      parameter == "DISTRACT_CLIMATE_z" ~ "Teacher-reported distraction\n(school-level main effect)",
      parameter == "INET_SCHOOL_z:DISTRACT_CLIMATE_z" ~ "Internet use × Distraction climate\n(cross-level interaction)",
      TRUE ~ parameter
    ),
    parameter = factor(parameter, levels = c(
      "Internet use for schoolwork\n(student-level main effect)",
      "Teacher-reported distraction\n(school-level main effect)",
      "Internet use × Distraction climate\n(cross-level interaction)"
    )),
    signif = ifelse(abs(estimate / se) > 1.96, "p < .05", "n.s."),
    outcome = factor(outcome, levels = c("Mathematics", "Science"))
  )

fig_p3 <- ggplot(p3_plot, aes(x = estimate, y = parameter, colour = outcome, shape = signif)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey60", linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.25, linewidth = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 4.5, position = position_dodge(width = 0.6)) +
  scale_colour_manual(
    values = c("Mathematics" = "#1D3557", "Science" = "#E76F51"),
    name = "Outcome"
  ) +
  scale_shape_manual(
    values = c("p < .05" = 16, "n.s." = 1),
    name = "Significance"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    title = "Prong 3: Teacher-reported distraction climate as school-level moderator",
    subtitle = paste0(
      "Pooled regression with country fixed effects. Filled = significant (p < .05); open = not significant.\n",
      "The cross-level interaction is negative but not statistically significant for either outcome."
    ),
    x = "Coefficient (achievement score points)",
    y = NULL,
    caption = paste0(
      "Source: TIMSS 2023. Adjusted for digital self-efficacy, gender, home educational resources, and country FE.\n",
      "Distraction climate = school mean of BTBG13G (distracted students), BTBM18C / BTBS21C (keeping students on task)."
    )
  ) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.subtitle = ggtext::element_textbox_simple(
      size = 10, colour = "grey40", lineheight = 1.2, margin = margin(b = 10))
  )

save_both(fig_p3, "fig_prong3_interaction", 10, 6)
cat("Prong 3 figure saved.\n")
cat("\nProng 3 Math:\n"); print(p3m)
cat("\nProng 3 Science:\n"); print(p3s)

# Store meta-analysis objects for later use
r_p1m <- list(meta = meta_p1[["(a) Mathematics"]])
r_p1s <- list(meta = meta_p1[["(b) Science"]])
r_p2m <- list(meta = meta_p2[["(a) Mathematics"]])
r_p2s <- list(meta = meta_p2[["(b) Science"]])

timss_region_lookup <- tribble(
  ~country_label,        ~region,
  "Chinese Taipei",      "East Asia",
  "Hong Kong SAR",       "East Asia",
  "Japan",               "East Asia",
  "Korea",               "East Asia",
  "Singapore",           "East Asia",

  "Uzbekistan",          "Central Asia",
  "Kazakhstan",          "Central Asia",
  "Azerbaijan",          "Central Asia",
  "Georgia",             "Central Asia",

  "Czech Republic",      "Europe",
  "England",             "Europe",
  "Hungary",             "Europe",
  "Portugal",            "Europe",
  "France",              "Europe",
  "Cyprus",              "Europe",
  "Italy",               "Europe",
  "Sweden",              "Europe",
  "Ireland",             "Europe",
  "Lithuania",           "Europe",
  "Austria",             "Europe",
  "Finland",             "Europe",
  "Romania",             "Europe",
  "Malta",               "Europe",

  "United States",       "North America",

  "Chile",               "Latin America",
  "Brazil",              "Latin America",

  "Australia",           "Oceania",
  "New Zealand",         "Oceania",

  "Morocco",             "Africa",
  "Côte d'Ivoire",       "Africa",
  "South Africa",        "Africa",

  "Palestine",           "Middle East",
  "Malaysia",            "Middle East",
  "Oman",                "Middle East",
  "Iran",                "Middle East",
  "Jordan",              "Middle East",
  "Saudi Arabia",        "Middle East",
  "Abu Dhabi (UAE)",     "Middle East",
  "Kuwait",              "Middle East",
  "Qatar",               "Middle East",
  "Bahrain",             "Middle East",
  "Israel",              "Middle East",
  "Sharjah (UAE)",       "Middle East",
  "United Arab Emirates","Middle East",
  "Dubai (UAE)",         "Middle East",
  "Türkiye",             "Middle East"
)

# extend palette only if needed
timss_region_cols <- pal_region
if (!"Africa" %in% names(timss_region_cols))  timss_region_cols["Africa"]  <- "#6D597A"
if (!"Oceania" %in% names(timss_region_cols)) timss_region_cols["Oceania"] <- "#43AA8B"
if (!"Other" %in% names(timss_region_cols))   timss_region_cols["Other"]   <- "grey60"
if (!"Pooled" %in% names(timss_region_cols))  timss_region_cols["Pooled"]  <- "black"

make_timss_panel_forest <- function(df_list, param_list, panel_labels,
                                    main_title, xlab, stub,
                                    width = 12, height = 9) {

  panels <- list()
  meta_results <- list()

  for (i in seq_along(df_list)) {
    fd <- df_list[[i]] %>%
      label_timss() %>%
      filter(parameter == param_list[i], !is.na(estimate), !is.na(se), se > 0) %>%
      left_join(timss_region_lookup, by = "country_label") %>%
      mutate(
        region = ifelse(is.na(region) | region == "", "Other", region),
        ci_lo = estimate - 1.96 * se,
        ci_hi = estimate + 1.96 * se,
        sig   = (ci_lo > 0 | ci_hi < 0),
        panel = panel_labels[i]
      )

    if (nrow(fd) < 3) next

    meta <- metafor::rma(
      yi     = estimate,
      sei    = se,
      data   = fd,
      method = "REML",
      test   = "knha"
    )
    meta_results[[panel_labels[i]]] <- meta

    pooled <- tibble(
      country_label = "Pooled",
      estimate = as.numeric(coef(meta)),
      ci_lo = meta$ci.lb,
      ci_hi = meta$ci.ub,
      se = meta$se,
      CNTRYID = "Pooled",
      sig = (meta$ci.lb > 0 | meta$ci.ub < 0),
      panel = panel_labels[i],
      region = "Pooled"
    )

    fd <- bind_rows(
      fd %>% select(country_label, estimate, ci_lo, ci_hi, se, CNTRYID, sig, panel, region),
      pooled
    )

    panels[[i]] <- fd
  }

  all_data <- bind_rows(panels) %>%
    mutate(panel = factor(panel, levels = panel_labels))

  # Shared ordering across facets: pooled bottom, higher estimates top
  country_order_nonpooled <- all_data %>%
    filter(country_label != "Pooled") %>%
    group_by(country_label) %>%
    summarise(order_value = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
    arrange(order_value) %>%
    pull(country_label)

  country_levels <- c("Pooled", country_order_nonpooled)

  all_data <- all_data %>%
    mutate(country_label = factor(country_label, levels = country_levels))

  guide_df <- tibble(country_label = factor(country_levels, levels = country_levels))

  # Put pooled estimate into facet titles
  panel_title_map <- vapply(panel_labels, function(lbl) {
    m <- meta_results[[lbl]]
    sprintf("%s\nPE = %.2f [%.2f, %.2f]", lbl, coef(m), m$ci.lb, m$ci.ub)
  }, character(1))

  # keep only colours actually used
  region_cols_use <- timss_region_cols[names(timss_region_cols) %in% unique(all_data$region)]

  x_min <- min(all_data$ci_lo, na.rm = TRUE)
  x_max <- max(all_data$ci_hi, na.rm = TRUE)
  x_pad <- 0.08 * (x_max - x_min)

  fig <- ggplot(all_data, aes(x = estimate, y = country_label)) +
    geom_hline(
      data = guide_df,
      aes(yintercept = country_label),
      colour = "grey75",
      linewidth = 0.25,
      alpha = 0.20
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey55"
    ) +
    geom_segment(
      aes(x = ci_lo, xend = ci_hi, yend = country_label, colour = region),
      linewidth = 0.9,
      alpha = 0.9,
      lineend = "round"
    ) +
    # open dots for all estimates
    geom_point(
      aes(colour = region),
      shape = 21,
      fill = "white",
      size = 2.8,
      stroke = 0.8
    ) +
    # filled dots only for significant estimates
    geom_point(
      data = dplyr::filter(all_data, sig),
      aes(fill = region, colour = region),
      shape = 21,
      size = 2.8,
      stroke = 0.8
    ) +
    facet_grid(
      . ~ panel,
      labeller = labeller(panel = as_labeller(panel_title_map))
    ) +
    scale_colour_manual(values = region_cols_use, name = "Region") +
    scale_fill_manual(values = region_cols_use, guide = "none") +
    coord_cartesian(
      xlim = c(x_min - x_pad, x_max + x_pad),
      clip = "off"
    ) +
    labs(
      title = main_title,
      subtitle = "Filled dots indicate 95% confidence intervals excluding zero; open dots indicate non-significant estimates.",
      x = xlab,
      y = NULL,
      caption = "Source: TIMSS 2023. Adjusted for digital self-efficacy, gender, and home educational resources."
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.7),
      axis.line = element_line(colour = "grey35", linewidth = 0.4),
      axis.ticks = element_line(colour = "grey35", linewidth = 0.4),
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.7),
      strip.text = element_text(face = "bold", size = 11),
      panel.spacing.x = unit(1.8, "lines"),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9, colour = "grey35"),
      plot.margin = margin(10, 20, 10, 10)
    )

  save_both(fig, stub, width, height)
  print(fig)

  cat("\n", main_title, "\n", sep = "")
  for (nm in names(meta_results)) {
    m <- meta_results[[nm]]
    cat(sprintf(
      "  %s: pooled = %.2f [%.2f, %.2f], I² = %.1f%%, k = %d\n",
      nm, coef(m), m$ci.lb, m$ci.ub, m$I2, m$k
    ))
  }

  invisible(meta_results)
}

# ── Prong 1 ───────────────────────────────────────────────────────────────────
meta_p1 <- make_timss_panel_forest(
  df_list      = list(p1_math_df, p1_sci_df),
  param_list   = c("GAME_FREQ_z", "GAME_FREQ_z"),
  panel_labels = c("(a) Mathematics", "(b) Science"),
  main_title   = "Prong 1: Learning games frequency (BSBG12F) and achievement",
  xlab         = "Achievement points per 1 SD increase in game frequency",
  stub         = "fig_prong1_panels",
  width = 12, height = 9
)

# ── Prong 2 ───────────────────────────────────────────────────────────────────
meta_p2 <- make_timss_panel_forest(
  df_list      = list(p2_math_df, p2_sci_df),
  param_list   = c("INET_SCHOOL_z", "INET_SCHOOL_z"),
  panel_labels = c("(a) Mathematics", "(b) Science"),
  main_title   = "Prong 2: Internet use for schoolwork (BSBG12A–F) and achievement",
  xlab         = "Achievement points per 1 SD increase in internet use",
  stub         = "fig_prong2_panels",
  width = 12, height = 9
)

if (!output_exists("robustness_binary_game_math.csv")) {
  # Binary: monthly game access
  run_bin_game <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$GAME_MONTHLY) & !is.na(d[[wc]])
    for (cv in timss_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 200) return(c(GAME_MONTHLY = NA_real_))
    m <- tryCatch(lm(as.formula(paste(oc, "~ GAME_MONTHLY +",
                                       paste(timss_covs, collapse = " + "),
                                       "+ factor(CNTRYID)")),
                      data = d, weights = d[[wc]]),
                  error = function(e) NULL)
    if (is.null(m)) return(c(GAME_MONTHLY = NA_real_))
    coef(m)["GAME_MONTHLY"]
  }
  rb_gm <- tryCatch(pv_jk2(tsg, PV_MATH, run_bin_game, "GAME_MONTHLY"),
                     error = function(e) NULL)
  if (!is.null(rb_gm)) write_csv(rb_gm, file.path(OUT_DIR, "robustness_binary_game_math.csv"))

  # Binary: top-tertile internet use
  run_bin_inet <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$INET_HIGH) & !is.na(d[[wc]])
    for (cv in timss_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 200) return(c(INET_HIGH = NA_real_))
    m <- tryCatch(lm(as.formula(paste(oc, "~ INET_HIGH +",
                                       paste(timss_covs, collapse = " + "),
                                       "+ factor(CNTRYID)")),
                      data = d, weights = d[[wc]]),
                  error = function(e) NULL)
    if (is.null(m)) return(c(INET_HIGH = NA_real_))
    coef(m)["INET_HIGH"]
  }
  rb_im <- tryCatch(pv_jk2(tsg, PV_MATH, run_bin_inet, "INET_HIGH"),
                     error = function(e) NULL)
  if (!is.null(rb_im)) write_csv(rb_im, file.path(OUT_DIR, "robustness_binary_inet_math.csv"))
  cat("TIMSS robustness checks saved.\n")
} else {
  cat("TIMSS robustness results already exist.\n")
}
