# RQ4 - overall time to treatment

# Meta-analysis of time-to-treatment
# Based on McGrath 'median of medians' approach
library(metafor)
library(tidyverse)
library(broom)
library(gt)
library(metamedian)


# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data files
trt_q_perc <- readRDS("Rdata/trt_q_perc.RDS")
trt_unadj_n <- readRDS("Rdata/trt_unadj_n.RDS")

# Useful functions
source("chemo_analysis0_repeated_functions.R")


# Table of percentiles for appendix  ------------------------------------------
table <- trt_q_perc |>
  filter(
    trt == "Chemo",
    stage == "All stages",
    variable == "all" | variable == "tumour_topography_group",
    percentile == 25 | percentile == 50 | percentile == 75 | percentile == 90
  ) |>
  rename(cancer = variable_value) |>
  mutate(cancer = ifelse(cancer == "all", "All 8 cancers" , cancer)) |>
  mutate(site_order = factor(cancer, levels = site)) |>
  mutate(j_order = factor(jurisdiction, levels = country_order)) |>
  select(site_order, j_order, percentile, days) |>
  pivot_wider(
    names_from = percentile,
    values_from = days,
    names_prefix = "pct_"
  ) |>
  arrange(site_order, j_order) |>
  group_by(site_order) |>
  chemo_exclusions_all(jurisdiction = j_order, cancer = site_order) |>
  select(site_order, j_order, pct_25, pct_50, pct_75, pct_90)

appdx_analysis3 <- table |> gt() |>
  cols_label (
    site_order = "Cancer site",
    j_order = "Jurisdiction",
    pct_25 = md("25<sup>th</sup> centile"),
    pct_50 = md("Median (50<sup>th</sup> centile)"),
    pct_75 = md("75<sup>th</sup> centile"),
    pct_90 = md("90<sup>th</sup> centile")
    
  ) |>
  fmt_number(
    columns = c(pct_25, pct_50, pct_75, pct_90),
    decimals = 1
  ) |>
  cols_align(align = "center") |>
  cols_align(align = "left", columns = j_order) |>
  tab_style(
    style = cell_text(
      font = "Arial",
      size = 8
    ),
    locations = list(
      cells_body(), 
      cells_column_labels()
    )
  )

appdx_analysis3 |> saveRDS("results/appdx_analysis3.RDS")

# Create dataset for meta-analysis ----------------------------------------
dat <- trt_q_perc |> 
  filter(
    stage == "All stages",
    trt == "Chemo",
    variable == "all" | 
      variable == "tumour_topography_group",
    percentile == 0 | 
      percentile == 25 | 
      percentile == 50 | 
      percentile == 75 | 
      percentile == 100
  ) |>
  rename(cancer = variable_value) |>
  mutate(cancer = ifelse(cancer == "all", "All 8 cancers" , cancer)) |>
  select(trt, jurisdiction, cancer, ordering, percentile, days) |>
  pivot_wider(
    names_from = percentile,
    values_from = days,
    names_prefix = "pct_"
  ) |>
  left_join(
    filter(
      trt_unadj_n,
      stage == "All stages",
      variable == "all"
    ) |>
      select(trt, cancer, jurisdiction, n_trt),
    by = c("trt", "cancer", "jurisdiction")
  ) |> 
  chemo_dat_clean() |>
  mutate(site_order  = factor(cancer, levels = site)) |>
  chemo_exclusions_all() |>
  filter(n_trt > 0)

# meta-analysis summary stats
poly_metan <- 
  dat |>
  filter(!exclude) |>
  group_by(cancer, trt) |>
  nest() |>
  mutate(
    metan = map(
      data, 
      qe_function1
    )
  ) |>
  mutate(
    metan_summ = map(
      data, 
      qe_function2
    )
  ) |>
  select(trt, cancer, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  mutate(
    pred_lb = prop-qt(.975, df = nobs)*sqrt(tau2+(se^2)),
    pred_ub = prop+qt(.975, df = nobs)*sqrt(tau2+(se^2))
  ) |>
  select(trt, cancer, prop0 = lower, prop1 = prop, prop2 = upper, prop3 = prop, pred_lb, pred_ub) |>
  pivot_longer(cols = starts_with("prop"), names_to = "prop", values_to = "x") |>
  mutate(y = ifelse(prop == "prop0" | prop == "prop2", 0, ifelse(prop == "prop1", 0.3, ifelse(prop == "prop3", -0.3, NA)))) |>
  mutate(site_order  = factor(cancer, levels = site)) |>
  mutate(trt_order   = factor(trt, levels = c("Chemo", "Radio"))) |>
  as_tibble()


# Create the meta-analysis table ------------------------------------------
table <- dat |>
  filter(!exclude) |>
  group_by(site_order) |>
  mutate(minimum = min(pct_50)) |>
  mutate(maximum = max(pct_50)) |>
  group_by(site_order, minimum, maximum) |>
  nest() |>
  mutate(
    metan = map(
      data, 
      qe_function1
    ) 
  ) |>
  mutate(
    metan_summ = map(
      data, 
      qe_function2
    )
  ) |>
  select(site_order, minimum, maximum, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  group_by() |>
  mutate(
    pred_lb = prop-qt(.975, df = nobs)*sqrt(tau2+(se^2)),
    pred_ub = prop+qt(.975, df = nobs)*sqrt(tau2+(se^2))
  ) |>
  mutate(
    prop  = prop,
    lower = lower,
    upper = upper,
    tau = sqrt(tau2)
  ) |> 
  select(site_order, prop, lower, upper, pred_lb, pred_ub, tau, i2, minimum, maximum) |>
  arrange(site_order) 

chemo_analysis3_table <- table |>
  gt() |>
  gt_metan_common() |>
  fmt_number(
    columns = c(tau),
    decimals = 1
  ) |>
  fmt_number(
    columns = c(prop, lower, upper, pred_lb, pred_ub, minimum, maximum),
    decimals = 1
  ) |>
  cols_align(align = "left", columns = site_order) |>
  cols_label (
    site_order = "Cancer site",
    prop = "Average days-to-treatment",
    lower = "(95% confidence interval)",
    pred_lb = "95% prediction interval *",
    minimum = "Observed jurisdictional range",
    tau = md("*Tau*, natural scale (days)  **"),
    i2 = md("*I<sup>2</sup>* ***")
  ) 

chemo_analysis3_table |> saveRDS("results/chemo_analysis3_table.RDS")


# Create the forest plot --------------------------------------------------

dat <- dat |>
  rename(
    lower = pct_25,
    upper = pct_75
  )

# Draw the forest plot
p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(yvar = ordering)
p <- p + incl_point(yvar = ordering, xvar = pct_50)
# Excluded jurisdictions and estimates
p <- p + excl_error(yvar = ordering)
p <- p + excl_point(yvar = ordering, xvar = pct_50)
# Meta-analysis
p <- p + meta_error(y_shift = 19) 
p <- p + meta_diamond(y_shift = 19, groupvar = site_order)
# y axis
p <- p + jdn_yscale()
# x axis
p <- p + days_xscale()
# legend
p <- p + guides(
  alpha = "none",
  colour = "none"
)
# faceting
p <- p + facet_wrap(site_order~., ncol = 3, scales = "free_x")
# themes
p <- p + theme_bw() +
  theme_icbp() +
  icbp_colour_manual()
p 
ggsave("results/figure3_chemo.svg",
       plot = p,
       width = 15,
       height = 15,
       units = "cm")
ggsave("results/figure3_chemo.pdf",
       plot = p,
       width = 15,
       height = 15,
       units = "cm")


# Alternative layout
p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(yvar = as.numeric(site_order))
p <- p + incl_point(yvar = as.numeric(site_order), xvar = pct_50)
# Excluded jurisdictions and estimates
p <- p + excl_error(yvar = as.numeric(site_order))
p <- p + excl_point(yvar = as.numeric(site_order), xvar = pct_50)
# Meta-analysis
p <- p + meta_error(y_shift = as.numeric(site_order)) 
p <- p + meta_diamond(y_shift = as.numeric(site_order), groupvar = site_order)
# y axis
p <- p + site_yscale_time()
# x axis
p <- p + days_xscale()
# legend
p <- p + guides(
  alpha = "none",
  colour = "none"
)
# faceting
p <- p + facet_wrap(ctry_order~., ncol = 4, scales = "free_x")
# themes
p <- p + theme_bw() +
  theme_icbp() +
  icbp_colour_manual()
p <- p + theme(strip.text.x=element_text(size=7)) 
p 
ggsave("results/figure3_chemo_appendix.svg",
       plot = p,
       width = 15,
       height = 13,
       units = "cm")
ggsave("results/figure3_chemo_appendix.pdf",
       plot = p,
       width = 15,
       height = 13,
       units = "cm")


# Clean up ----------------------------------------------------------------
rm(list = ls())

