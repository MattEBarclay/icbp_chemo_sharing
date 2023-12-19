# Meta-analysis of quantile regression results
# Treating the quantile summaries as though they were linear regression

# RQ5 - are the differences in time to treatment similar by jurisdiction?

library(metafor)
library(tidyverse)
library(broom)
library(gt)

# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data files
trt_q <- readRDS("Rdata/trt_q.RDS")

# Useful functions
source("chemo_analysis0_repeated_functions.R")

# Create dataset for meta-analysis ----------------------------------------
dat <- trt_q |> 
  filter(
    stage == "All stages",
    model == "CHM_quantile_nostage" ,
      term == "15-64 vs 65-74"   | 
      term == "75-84 vs 65-74"   |
      term == "85-99 vs 65-74"   |
      term == "Women vs Men"     |
      term == "Oesophageal vs Lung"  |
      term == "Stomach vs Lung"  |
      term == "Colon vs Lung"    |
      term == "Rectal vs Lung"   |
      term == "Liver vs Lung"    |
      term == "Pancreatic vs Lung" | 
      term == "Ovarian vs Lung",
    jurisdiction != ""
  ) |> 
  mutate(term = factor(term, levels = term_list)) |>
  select(trt, jurisdiction, ordering, term, q50, q50_lower, q50_upper, q50_se) |> 
  arrange(trt, ordering) |>
  chemo_dat_clean() |>
  chemo_exclusions_comparisons() |>
  arrange(term, ordering)


# meta-analysis summary stats
poly_metan <- dat |>
  filter(!exclude) |>
  filter(
    term == "15-64 vs 65-74"   | 
      term == "75-84 vs 65-74"   |
      term == "85-99 vs 65-74"   |
      term == "Women vs Men"     
  ) |>
  group_by(trt, term) |>
  nest() |>
  mutate(
    metan = map(
      data, 
      ~ rma(yi = q50, sei = q50_se, data = .x, measure = "GEN", method = "REML") |>
        tidy(conf.int = TRUE) |>
        select(prop = estimate, lower = conf.low, upper = conf.high, se = std.error)
    )
  ) |>
  mutate(metan_summ = map(
    data, 
    ~ rma(yi = q50, sei = q50_se, data = .x, measure = "GEN", method = "REML") |>
      glance() |>
      select(tau2 = tau.squared, nobs = nobs)
  )) |>
  select(trt, term, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  group_by(trt) |>
  mutate(
    pred_lb = prop-qt(.975, df = nobs)*sqrt(tau2+(se^2)),
    pred_ub = prop+qt(.975, df = nobs)*sqrt(tau2+(se^2))
  ) |>
  select(trt, term, pred_lb, pred_ub, prop0 = lower, prop1 = prop, prop2 = upper, prop3 = prop) |>
  pivot_longer(cols = starts_with("prop"), names_to = "prop", values_to = "x") |>
  mutate(y = ifelse(prop == "prop0" | prop == "prop2", 0, ifelse(prop == "prop1", 0.3, ifelse(prop == "prop3", -0.3, NA)))) |>
  mutate(trt_order   = factor(trt, levels = c("Chemo", "Radio")))  |>
  mutate(term        = factor(term, levels= c(
    "Women vs Men", 
    "15-64 vs 65-74", 
    "75-84 vs 65-74", 
    "85-99 vs 65-74"
  ))) |>
  as_tibble()


# Appendix table of relevant data ------------------------------------------
appdx_analysis4 <- dat |> 
  filter(
    term == "15-64 vs 65-74"   | 
      term == "75-84 vs 65-74"   |
      term == "85-99 vs 65-74"   |
      term == "Women vs Men"     
  ) |>
  select(term, jurisdiction, q50, q50_lower, q50_upper) |> 
  group_by(term) |>
  gt() |>
  cols_merge(
    columns = c(q50_lower, q50_upper),
    pattern = "({1}, {2})"
  ) |>
  cols_label (
    jurisdiction = "Jurisdiction", 
    term = "Comparison",
    q50 = "Median difference",
    q50_lower = "(95% confidence interval)"
  ) |>
  fmt_number(
    columns = c(q50, q50_lower, q50_upper),
    decimals = 1
  ) |>
  cols_align(align = "center") |>
  cols_align(align = "left", columns = c(jurisdiction)) |>
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

appdx_analysis4 |> saveRDS("results/appdx_analysis4.RDS")

# Create the meta-analysis table ------------------------------------------

# Code to create table
table <- dat |>
  filter(!exclude) |>
  group_by(term) |>
  mutate(minimum = min(q50)) |>
  mutate(maximum = max(q50)) |>
  group_by(term, minimum, maximum) |>
  nest() |>
  mutate(
    metan = map(
      data, 
      ~ rma(yi = q50, sei = q50_se, data = .x, measure = "GEN", method = "REML") |>
        tidy(conf.int = TRUE) |>
        select(prop = estimate, lower = conf.low, upper = conf.high, se = std.error))
  ) |>
  mutate(metan_summ = map(
    data, ~ rma(yi = q50, sei = q50_se, data = .x, measure = "GEN", method = "REML") |>
      glance() |>
      select(i2 = i.squared, tau2 = tau.squared, nobs = nobs)
  )) |>
  select(term, minimum, maximum, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  group_by() |>
  arrange(term) |>
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
  select(term, prop, lower, upper, pred_lb, pred_ub, tau, i2, minimum, maximum) 

chemo_analysis4_table <- table |>
  gt() |>
  gt_metan_common() |>
  fmt_number(
    columns = c(prop, lower, upper, pred_lb, pred_ub, minimum, maximum),
    decimals = 1
  ) |>
  fmt_number(
    columns = c(tau),
    decimals = 1
  ) |>
  cols_align(align = "left", columns = term) |>
  cols_label (
    term = "",
    prop = "Difference in median days-to-treatment",
    lower = "(95% confidence interval)",
    pred_lb = "95% prediction interval *",
    minimum = "Observed jurisdictional range",
    tau = md("*Tau*, natural scale (days) **"),
    i2 = md("*I<sup>2</sup>* ***")
  ) 

chemo_analysis4_table |> saveRDS("results/chemo_analysis4_table.RDS")


# Create the forest plot --------------------------------------------------

# want to exclude highly uncertain quantiles
# so the plot is more readable
dat |> 
  filter(
    term == "15-64 vs 65-74"   | 
      term == "75-84 vs 65-74"   |
      term == "85-99 vs 65-74"   |
      term == "Women vs Men"     
  ) |>
  filter(q50_upper-q50_lower <= 30) |> 
  select(q50_lower) |> 
  min()

dat |> 
  filter(
    term == "15-64 vs 65-74"   | 
      term == "75-84 vs 65-74"   |
      term == "85-99 vs 65-74"   |
      term == "Women vs Men"     
  ) |>
  filter(q50_upper-q50_lower <= 30) |> 
  select(q50_upper) |> 
  max()

# some convenient edits to the data
dat <- dat |> 
  filter(
    term == "15-64 vs 65-74"   | 
      term == "75-84 vs 65-74"   |
      term == "85-99 vs 65-74"   |
      term == "Women vs Men"     
  ) |>
  filter(q50_upper-q50_lower <= 30) |> 
  mutate(
    term = factor(
      term, 
      levels= c(
        "Women vs Men", 
        "15-64 vs 65-74", 
        "75-84 vs 65-74", 
        "85-99 vs 65-74"
      )
    )
  ) |>
  rename(
    lower = q50_lower,
    upper = q50_upper
  )


# Draw the forest plot
p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(yvar = ordering)
p <- p + incl_point(yvar = ordering, xvar = q50)
# Excluded jurisdictions and estimates
p <- p + excl_error(yvar = ordering)
p <- p + excl_point(yvar = ordering, xvar = q50)
# Meta-analysis
p <- p + meta_error(y_shift = 19) 
p <- p + meta_diamond(y_shift = 19, groupvar = term)
# y axis
p <- p + jdn_yscale()
# x axis
p <- p + daydiff_xscale()
p <- p + geom_vline(xintercept = 0)
# legend
p <- p + guides(
  alpha = "none",
  colour = "none"
)
# faceting
p <- p + facet_wrap(term~., ncol = 2, scales = "free_x")
# themes
p <- p + theme_bw() +
  theme_icbp() +
  icbp_colour_manual()
p 
ggsave("results/figure4_chemo_appendix.svg",
       plot = p,
       width = 15,
       height = 11.25,
       units = "cm")
ggsave("results/figure4_chemo_appendix.pdf",
       plot = p,
       width = 15,
       height = 11.25,
       units = "cm")


# Clean up ----------------------------------------------------------------
rm(list = ls())

