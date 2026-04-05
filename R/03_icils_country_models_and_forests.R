# ICILS country models and figures

# Requires 00_setup.R, 01_load_icils.R, and 02_icils_descriptives_and_figures.R.

if (!output_exists("table_country_coefs.csv")) {
  rhs <- paste(c("S_ACMULT", available_covs), collapse = " + ")
  cnames <- c("(Intercept)", "S_ACMULT", available_covs)

  run_model <- function(d, oc, wc) {
    keep <- !is.na(d[[oc]]) & !is.na(d$S_ACMULT) & !is.na(d[[wc]])
    for (cv in available_covs) keep <- keep & !is.na(d[[cv]])
    d <- d[keep, , drop = FALSE]
    if (nrow(d) < 50) return(setNames(rep(NA_real_, length(cnames)), cnames))
    m <- tryCatch(lm(as.formula(paste(oc, "~", rhs)), data = d, weights = d[[wc]]),
                  error = function(e) NULL)
    if (is.null(m)) return(setNames(rep(NA_real_, length(cnames)), cnames))
    cc <- coef(m)
    out <- setNames(rep(NA_real_, length(cnames)), cnames)
    out[names(cc)] <- cc
    out
  }

  cil_results <- list(); ct_results <- list()
  for (cnt in countries) {
    cd <- bsg[bsg$CNTRYID == cnt, , drop = FALSE]
    # CIL
    if (sum(!is.na(cd$PV1CIL)) >= 100) {
      cat(sprintf("  CIL: %s ...", cnt))
      res <- tryCatch(pv_jrr(cd, PV_CIL, run_model), error = function(e) NULL)
      if (!is.null(res)) {
        cil_results[[cnt]] <- res %>% mutate(CNTRYID = cnt, outcome = "CIL")
        cat(" OK\n")
      } else cat(" ERROR\n")
    }
    # CT
    if (sum(!is.na(cd$PV1CT)) >= 100) {
      cat(sprintf("  CT:  %s ...", cnt))
      res <- tryCatch(pv_jrr(cd, PV_CT, run_model), error = function(e) NULL)
      if (!is.null(res)) {
        ct_results[[cnt]] <- res %>% mutate(CNTRYID = cnt, outcome = "CT")
        cat(" OK\n")
      } else cat(" ERROR\n")
    }
  }

  all_results <- bind_rows(bind_rows(cil_results), bind_rows(ct_results)) %>%
    left_join(region_map, by = "CNTRYID")
  write_csv(all_results, file.path(OUT_DIR, "table_country_coefs.csv"))
  cat("Country regression results saved.\n")
} else {
  all_results <- read_csv(file.path(OUT_DIR, "table_country_coefs.csv"),
                          show_col_types = FALSE)
  cat("Loaded cached table_country_coefs.csv\n")
}

cil_df <- all_results %>% filter(outcome == "CIL")
ct_df  <- all_results %>% filter(outcome == "CT")

pick_first_present <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

code_var_cil <- pick_first_present(cil_df, c("CNTRYID", "country_id", "idcountry", "code", "CNT"))
code_var_ct  <- pick_first_present(ct_df,  c("CNTRYID", "country_id", "idcountry", "code", "CNT"))

if (is.na(code_var_cil) || is.na(code_var_ct)) {
  stop("Could not identify the country code column in cil_df and/or ct_df.")
}

# Region labels
lookup_df <- region_map %>%
  transmute(
    code = as.character(CNTRYID),
    country_map = as.character(country_label),
    region_map2 = as.character(region)
  ) %>%
  distinct()

prep_pair_df <- function(df, outcome, code_var, lookup_df) {
  df %>%
    filter(parameter == "S_ACMULT", !is.na(estimate), !is.na(se), se > 0) %>%
    mutate(
      code = as.character(.data[[code_var]]),
      outcome = outcome,
      ci_lb = estimate - 1.96 * se,
      ci_ub = estimate + 1.96 * se
    ) %>%
    select(-any_of(c("country", "country_label", "region"))) %>%
    left_join(lookup_df, by = "code") %>%
    mutate(
      country = dplyr::coalesce(country_map, code),
      region  = dplyr::coalesce(region_map2, "Other")
    ) %>%
    select(code, country, region, outcome, estimate, se, ci_lb, ci_ub)
}

forest_cil_pair <- prep_pair_df(cil_df, "CIL", code_var_cil, lookup_df)
forest_ct_pair  <- prep_pair_df(ct_df,  "CT",  code_var_ct,  lookup_df)

# Keep countries present in both outcomes
pair_df <- bind_rows(forest_cil_pair, forest_ct_pair) %>%
  group_by(code, country, region) %>%
  filter(n_distinct(outcome) == 2) %>%
  ungroup()

# REML and Knapp-Hartung pooled estimates
meta_cil_pair <- rma(
  yi = estimate,
  sei = se,
  data = forest_cil_pair,
  method = "REML",
  test = "knha"
)

meta_ct_pair <- rma(
  yi = estimate,
  sei = se,
  data = forest_ct_pair,
  method = "REML",
  test = "knha"
)

pooled_rows <- bind_rows(
  tibble(
    code = "Pooled",
    country = "Pooled",
    region = "Pooled",
    outcome = "CIL",
    estimate = as.numeric(coef(meta_cil_pair)),
    se = meta_cil_pair$se,
    ci_lb = meta_cil_pair$ci.lb,
    ci_ub = meta_cil_pair$ci.ub
  ),
  tibble(
    code = "Pooled",
    country = "Pooled",
    region = "Pooled",
    outcome = "CT",
    estimate = as.numeric(coef(meta_ct_pair)),
    se = meta_ct_pair$se,
    ci_lb = meta_ct_pair$ci.lb,
    ci_ub = meta_ct_pair$ci.ub
  )
)

plot_df <- bind_rows(pair_df, pooled_rows) %>%
  mutate(
    sig = (ci_lb > 0 | ci_ub < 0)
  )

# Order countries for plotting
country_order_nonpooled <- pair_df %>%
  group_by(country, region) %>%
  summarise(order_value = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  arrange(order_value) %>%
  pull(country)

country_levels <- c("Pooled", country_order_nonpooled)

plot_df <- plot_df %>%
  mutate(
    country = factor(country, levels = country_levels),
    outcome = factor(outcome, levels = c("CT", "CIL"))
  )

# Facet titles with pooled estimates
facet_title_map <- plot_df %>%
  filter(country == "Pooled") %>%
  mutate(
    facet_title = dplyr::case_when(
      outcome == "CT"  ~ sprintf("A. CT\nPE = %.2f [%.2f, %.2f]", estimate, ci_lb, ci_ub),
      outcome == "CIL" ~ sprintf("B. CIL\nPE = %.2f [%.2f, %.2f]", estimate, ci_lb, ci_ub)
    )
  ) %>%
  select(outcome, facet_title) %>%
  tibble::deframe()

# Region colors
if (exists("pal_region")) {
  region_cols <- pal_region
} else {
  reg_names <- unique(plot_df$region)
  region_cols <- setNames(scales::hue_pal()(length(reg_names)), reg_names)
}
if (!"Pooled" %in% names(region_cols)) region_cols <- c(region_cols, Pooled = "black")
if (!"Other"  %in% names(region_cols)) region_cols  <- c(region_cols, Other = "grey60")

# Plot limits
x_min <- min(plot_df$ci_lb, na.rm = TRUE)
x_max <- max(plot_df$ci_ub, na.rm = TRUE)
x_pad <- 0.08 * (x_max - x_min)
right_lim <- x_max + x_pad

fig2ab_facets <- ggplot(plot_df, aes(x = estimate, y = country)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey55"
  ) +
  geom_segment(
    aes(x = ci_lb, xend = ci_ub, yend = country, colour = region),
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
    data = dplyr::filter(plot_df, sig),
    aes(fill = region, colour = region),
    shape = 21,
    size = 2.8,
    stroke = 0.8
  ) +
  facet_grid(
    . ~ outcome,
    labeller = labeller(outcome = as_labeller(facet_title_map))
  ) +
  scale_colour_manual(values = region_cols, name = "Region") +
  scale_fill_manual(values = region_cols, guide = "none") +
  coord_cartesian(
    xlim = c(x_min - x_pad, right_lim),
    clip = "off"
  ) +
  labs(
    title = "Country-specific associations of academic-media multitasking with CIL and CT",
    subtitle = "Filled dots indicate estimates with 95% confidence intervals excluding zero; open dots indicate non-significant estimates",
    x = "Regression coefficient for S_ACMULT",
    y = NULL,
    caption = "Horizontal lines show 95% confidence intervals. Pooled estimates use REML random-effects meta-analysis with Knapp-Hartung adjustment."
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.7),
    axis.line = element_line(colour = "grey35", linewidth = 0.4),
    axis.ticks = element_line(colour = "grey35", linewidth = 0.4),
    axis.text.y = element_text(size = 8.5),
    strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing.x = unit(1.8, "lines"),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 9, colour = "grey35"),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(file.path(OUT_DIR, "fig4_forest_ct_cil.png"), fig2ab_facets, width = 12, height = 9, dpi = 300)
ggsave(file.path(OUT_DIR, "fig4_forest_ct_cil.pdf"), fig2ab_facets, width = 12, height = 9)
print(fig2ab_facets)

# Plot order for item and grouped estimates

if (!output_exists("fig_pooled_item_grouped_forest.png")) {

  # ── Define exposures ──
  pooled_terms <- tribble(
    ~var,                    ~label,                                            ~type,      ~order,
    "ITEM_A_Z",              "A. Text chat with others",                        "Item",     1,
    "ITEM_B_Z",              "B. Use social media to view/post content",        "Item",     2,
    "ITEM_C_Z",              "C. Check social media",                           "Item",     3,
    "ITEM_D_Z",              "D. Use the internet to find information",         "Item",     4,
    "ITEM_E_Z",              "E. Watch videos / live streams / television",     "Item",     5,
    "ITEM_F_Z",              "F. Listen to music / radio / podcasts",           "Item",     6,
    "ACMULT_SOCIAL_ABC_Z",   "A/B/C. Social interruption",                     "Grouped",  7,
    "ACMULT_INFO_D_Z",       "D. Information seeking",                         "Grouped",  8,
    "ACMULT_ENT_EF_Z",       "E/F. Entertainment / background media",          "Grouped",  9,
    "ACMULT_OFFTASK_Z",      "A+B+C+E+F. Off-task multitasking",              "Grouped",  10,
    "S_ACMULT_Z",            "Overall S_ACMULT scale",                         "Overall",  11
  )

  # ── Create item-level z-scores if not already present ──
  item_vars <- c(A = "IS3G21A", B = "IS3G21B", C = "IS3G21C",
                 D = "IS3G21D", E = "IS3G21E", F = "IS3G21F")
  item_vars <- item_vars[item_vars %in% names(bsg)]

  for (nm in names(item_vars)) {
    zname <- paste0("ITEM_", nm, "_Z")
    if (!zname %in% names(bsg))
      bsg[[zname]] <- as.numeric(scale(bsg[[item_vars[[nm]]]]))
  }

  # Grouped composites
  social_items <- unname(item_vars[c("A","B","C")[c("A","B","C") %in% names(item_vars)]])
  ent_items    <- unname(item_vars[c("E","F")[c("E","F") %in% names(item_vars)]])

  if (!"ACMULT_SOCIAL_ABC_Z" %in% names(bsg))
    bsg$ACMULT_SOCIAL_ABC_Z <- z_rowmean(bsg, social_items)
  if (!"ACMULT_INFO_D_Z" %in% names(bsg) && "IS3G21D" %in% names(bsg))
    bsg$ACMULT_INFO_D_Z <- as.numeric(scale(bsg$IS3G21D))
  if (!"ACMULT_ENT_EF_Z" %in% names(bsg))
    bsg$ACMULT_ENT_EF_Z <- z_rowmean(bsg, ent_items)
  if (!"ACMULT_OFFTASK_Z" %in% names(bsg))
    bsg$ACMULT_OFFTASK_Z <- z_rowmean(bsg, c(social_items, ent_items))
  if (!"S_ACMULT_Z" %in% names(bsg))
    bsg$S_ACMULT_Z <- as.numeric(scale(bsg$S_ACMULT))

  # Keep only terms whose variables exist
  pooled_terms <- pooled_terms %>% filter(var %in% names(bsg))

  # ── Covariates for pooled model ──
  pooled_covs <- unique(c(available_covs))
  pooled_covs <- pooled_covs[pooled_covs %in% names(bsg) &
    sapply(pooled_covs, function(v) sum(!is.na(bsg[[v]])) > 100)]

  # ── Pooled FE model function ──
  run_pooled_fe <- function(data, outcome_col, weight_col, exposure_var, covars) {
    keep <- !is.na(data[[outcome_col]]) & !is.na(data[[exposure_var]]) &
            !is.na(data[[weight_col]]) & !is.na(data$CNTRYID)
    for (cv in covars) keep <- keep & !is.na(data[[cv]])
    d <- data[keep, , drop = FALSE]
    param_names <- c("(Intercept)", exposure_var, covars)
    if (nrow(d) < 500)
      return(setNames(rep(NA_real_, length(param_names)), param_names))
    f <- as.formula(paste(outcome_col, "~",
                          paste(c(exposure_var, covars), collapse = " + "),
                          "+ factor(CNTRYID)"))
    m <- tryCatch(lm(f, data = d, weights = d[[weight_col]]), error = function(e) NULL)
    if (is.null(m)) return(setNames(rep(NA_real_, length(param_names)), param_names))
    cc <- coef(m)
    out <- setNames(rep(NA_real_, length(param_names)), param_names)
    out[names(cc)] <- cc; out
  }

  fit_pooled <- function(data, pv_names, exposure_var, covars) {
    res <- tryCatch(
      pv_jrr(data, pv_names, function(d, oc, wc)
        run_pooled_fe(d, oc, wc, exposure_var, covars)),
      error = function(e) NULL)
    if (is.null(res)) return(NULL)
    res %>% filter(parameter == exposure_var)
  }

  # ── Run all models ──
  pooled_item_results <- list()
  for (i in seq_len(nrow(pooled_terms))) {
    v   <- pooled_terms$var[i]
    lab <- pooled_terms$label[i]
    typ <- pooled_terms$type[i]
    ord <- pooled_terms$order[i]
    cat(sprintf("  Pooled model: %s ... ", v))

    cil_res <- fit_pooled(bsg, PV_CIL, v, pooled_covs)
    if (!is.null(cil_res))
      pooled_item_results[[paste0("CIL_", v)]] <- cil_res %>%
        mutate(outcome = "B. CIL", var = v, label = lab, type = typ, order = ord)

    if (sum(!is.na(bsg$PV1CT)) > 1000) {
      ct_res <- fit_pooled(bsg, PV_CT, v, pooled_covs)
      if (!is.null(ct_res))
        pooled_item_results[[paste0("CT_", v)]] <- ct_res %>%
          mutate(outcome = "A. CT", var = v, label = lab, type = typ, order = ord)
    }
    cat("done\n")
  }

  pooled_item_df <- bind_rows(pooled_item_results) %>%
    filter(!is.na(estimate), !is.na(se))
  write_csv(pooled_item_df, file.path(OUT_DIR, "table_pooled_item_grouped.csv"))

} else {
  pooled_item_df <- read_csv(file.path(OUT_DIR, "table_pooled_item_grouped.csv"),
                             show_col_types = FALSE)
  cat("Loaded cached pooled item/grouped results.\n")
}

ordered_labels <- c(
  "A. Text chat with others",
  "B. Use social media to view/post content",
  "C. Check social media",
  "D. Use the internet to find information",
  "E. Watch videos / live streams / television",
  "F. Listen to music / radio / podcasts",
  "A/B/C. Social interruption",
  "D. Information seeking",
  "E/F. Entertainment / background media",
  "A+B+C+E+F. Off-task multitasking",
  "Overall S_ACMULT scale"
)
ordered_labels <- ordered_labels[ordered_labels %in% pooled_item_df$label]

pooled_item_df_plot <- pooled_item_df %>%
  mutate(
    outcome = factor(outcome, levels = c("A. CT", "B. CIL")),
    type    = factor(type, levels = c("Item", "Grouped", "Overall")),
    label   = factor(label, levels = rev(ordered_labels)),
    sig     = (ci_lo > 0 | ci_hi < 0)
  )

n_item    <- sum(pooled_item_df_plot$type == "Item"    & !duplicated(paste(pooled_item_df_plot$label, pooled_item_df_plot$type)))
n_grouped <- sum(pooled_item_df_plot$type == "Grouped" & !duplicated(paste(pooled_item_df_plot$label, pooled_item_df_plot$type)))
n_overall <- sum(pooled_item_df_plot$type == "Overall" & !duplicated(paste(pooled_item_df_plot$label, pooled_item_df_plot$type)))

sep1 <- n_overall + 0.5
sep2 <- n_overall + n_grouped + 0.5

block_df <- expand.grid(
  outcome = levels(pooled_item_df_plot$outcome), stringsAsFactors = FALSE
) %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(
    block = rep(c("Overall", "Grouped", "Item"), times = 2),
    ymin  = rep(c(0.5, sep1, sep2), times = 2),
    ymax  = rep(c(sep1, sep2, n_overall + n_grouped + n_item + 0.5), times = 2),
    xmin  = -Inf, xmax = Inf
  )

sep_df <- expand.grid(
  outcome = levels(pooled_item_df_plot$outcome),
  yint = c(sep1, sep2),
  stringsAsFactors = FALSE
)

x_rng <- range(c(pooled_item_df_plot$ci_lo, pooled_item_df_plot$ci_hi), na.rm = TRUE)
x_lim <- c(floor(x_rng[1] - 0.5), ceiling(x_rng[2] + 0.5))

fig_items <- ggplot(
  pooled_item_df_plot,
  aes(x = estimate, y = label)
) +
  geom_rect(
    data = block_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = block),
    inherit.aes = FALSE, alpha = 0.06, colour = NA
  ) +
  geom_hline(
    data = sep_df,
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linewidth = 0.7,
    colour = "grey60"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    colour = "grey50",
    linetype = "dashed"
  ) +
  geom_errorbarh(
    aes(xmin = ci_lo, xmax = ci_hi, colour = type),
    height = 0.22,
    linewidth = 0.6
  ) +
  geom_point(
    aes(colour = type),
    shape = 21,
    fill = "white",
    size = 3,
    stroke = 0.9
  ) +
  geom_point(
    data = pooled_item_df_plot %>% filter(sig),
    aes(fill = type, colour = type),
    shape = 21,
    fill = "black",
    size = 3,
    stroke = 0.9
  ) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_colour_manual(
    values = c("Item" = "#457B9D", "Grouped" = "#1D3557", "Overall" = "#E76F51"),
    name = "Row type"
  ) +
  scale_fill_manual(
    values = c("Item" = "#457B9D", "Grouped" = "#1D3557", "Overall" = "#E76F51"),
    name = "Row type"
  ) +
  scale_shape_manual(
    values = c("Item" = 21, "Grouped" = 21, "Overall" = 21),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("Overall" = "grey45", "Grouped" = "grey60", "Item" = "grey75"),
    guide = "none"
  ) +
  coord_cartesian(xlim = x_lim) +
  labs(
    title = "Pooled associations between specific multitasking behaviours and digital proficiency",
    subtitle = paste0(
      "Items A–F at top, grouped families in the middle, overall S_ACMULT at bottom. ",
      "Filled points indicate 95% confidence intervals excluding zero."
    ),
    x = "Score points per 1 SD increase in exposure",
    y = NULL,
    caption = paste0(
      "Source: ICILS 2023. Country fixed-effect models. ",
      "Grouped rows are exploratory summaries; primary exposure = S_ACMULT."
    )
  ) +
  theme_paper(base_size = 12) +
  theme(
    panel.border     = element_rect(colour = "grey55", fill = NA, linewidth = 0.9),
    strip.background = element_rect(fill = "grey92", colour = "grey55", linewidth = 0.9),
    strip.text       = element_text(face = "bold", size = 12),
    panel.spacing.x  = unit(1.5, "lines"),
    axis.text.y      = element_text(size = 10),
    legend.position  = "bottom",
    legend.box       = "horizontal"
  )

save_both(fig_items, "fig_pooled_item_grouped_forest", 13, 9)
cat("Pooled item/grouped forest plot saved.\n")
print(fig_items)

if (!output_exists("fig7_heatmap.png")) {
  hd <- all_results %>%
    filter(parameter == "S_ACMULT", !is.na(estimate)) %>%
    select(country_label, outcome, estimate, p_value) %>%
    mutate(
      signif = case_when(p_value < .001 ~ "***", p_value < .01 ~ "**",
                         p_value < .05 ~ "*", TRUE ~ ""),
      label = sprintf("%.2f%s", estimate, signif))

  hc <- hd %>%
    complete(country_label, outcome) %>%
    left_join(hd %>% filter(outcome == "CIL") %>%
                select(country_label, cil_est = estimate),
              by = "country_label") %>%
    mutate(country_label = fct_reorder(country_label, cil_est, .na_rm = TRUE))

  fig7 <- ggplot(hc, aes(outcome, country_label, fill = estimate)) +
    geom_tile(data = hc %>% filter(!is.na(estimate)),
              colour = "white", linewidth = .6) +
    geom_text(data = hc %>% filter(!is.na(label)),
              aes(label = label), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#E76F51", mid = "#FEFAE0", high = "#264653",
                         midpoint = 0, name = "Beta", na.value = "transparent") +
    labs(title = "CIL vs CT coefficients across countries",
         subtitle = "*p<.05; **p<.01; ***p<.001",
         x = NULL, y = NULL,
         caption = "Source: ICILS 2023. PV x JRR estimates.") +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))
  save_both(fig7, "fig7_heatmap", 7.2, 10)
} else {
  cat("fig7_heatmap already exists, skipping.\n")
}

if (!output_exists("fig5_dose_response.png")) {
  dd <- bsg %>%
    filter(!is.na(PV1CIL), !is.na(S_ACMULT), !is.na(TOTWGTS), TOTWGTS > 0) %>%
    mutate(CNTRYID_f = factor(CNTRYID))
  gcovs <- intersect(c("SES_CAN", "EXCOMP_CAN", "GENDER_CAN"), names(dd))

  gf <- as.formula(paste("PV1CIL ~ s(S_ACMULT, bs='tp', k=4) +",
                         paste(gcovs, collapse = " + "),
                         "+ s(CNTRYID_f, bs='re')"))
  gm <- gam(gf, data = dd, weights = TOTWGTS / mean(TOTWGTS),
            method = "REML", select = TRUE)
  writeLines(capture.output(summary(gm)), file.path(OUT_DIR, "gam_summary.txt"))

  nd <- data.frame(S_ACMULT = seq(25, 75, .5))
  for (cv in gcovs) nd[[cv]] <- if (cv == "GENDER_CAN") 0 else median(dd[[cv]], na.rm = TRUE)
  nd$CNTRYID_f <- levels(dd$CNTRYID_f)[1]
  p <- predict(gm, nd, type = "response", exclude = "s(CNTRYID_f)", se.fit = TRUE)
  nd$pred <- p$fit; nd$se <- p$se.fit

  fig5 <- ggplot(nd, aes(S_ACMULT, pred)) +
    geom_ribbon(aes(ymin = pred - 1.96*se, ymax = pred + 1.96*se),
                fill = "#2A9D8F", alpha = .2) +
    geom_line(colour = "#264653", linewidth = 1) +
    geom_rug(data = dd %>% slice_sample(n = min(5000, nrow(dd))),
             aes(x = S_ACMULT), inherit.aes = FALSE, sides = "b", alpha = .03) +
    labs(title = "Dose-response: multitasking and CIL",
         subtitle = "GAM with thin-plate spline, adjusting for covariates + country RE.",
         x = "S_ACMULT", y = "Predicted CIL",
         caption = "Source: ICILS 2023 (PV1CIL). Shaded area = 95% CI.")
  save_both(fig5, "fig5_dose_response", 9, 6.5)
} else {
  cat("fig5_dose_response already exists, skipping.\n")
}
