## RQ2 - differences in treatment use by patient characteristic
library(metafor)
library(tidyverse)
library(broom)
library(gt)

# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data file
trt_yn <- readRDS("Rdata/trt_yn.RDS")

# Useful functions
source("chemo_analysis0_repeated_functions.R")

# Create dataset for meta-analysis ----------------------------------------
dat <- trt_yn |> 
  filter(
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
    model == "CHM_YN_all_nostage", 
    jurisdiction != "",
    stage == "All stages"
  ) |>
  mutate(lower = estimate-1.96*std.error) |>
  mutate(upper = estimate+1.96*std.error)  |>
  mutate(term = factor(term, levels = term_list)) |>
  chemo_dat_clean() |>
  chemo_exclusions_comparisons() |>
  arrange(term, ordering)


# Useful to have a 'stage' dataset for the countries with good
# stage data
dat_stage <- trt_yn |> 
  filter(
    term == "15-64 vs 65-74" | 
      term == "75-84 vs 65-74" |
      term == "85-99 vs 65-74" |
      term == "Women vs Men" ,
    model == "CHM_YN_all_stage", 
    jurisdiction == "Alberta" | 
      jurisdiction == "Nova Scotia" | 
      jurisdiction == "Saskatchewan" | 
      jurisdiction == "Manitoba" | 
      jurisdiction == "New South Wales",
    stage == "All stages"
  ) |>
  mutate(lower = estimate-1.96*std.error) |>
  mutate(upper = estimate+1.96*std.error)  |>
  mutate(term = factor(term, levels = term_list)) |>
  chemo_dat_clean() |>
  chemo_exclusions_comparisons() |>
  filter(std.error < 1) |>
  arrange(term, ordering)

# meta-analysis summary stats for plotting
poly_metan <- 
  dat |>
  filter(!exclude) |>
  group_by(term, trt) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(yi = estimate, vi = vi, data = .x, measure = "GEN", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop0 = conf.low, prop1 = estimate, prop2 = conf.high, prop3 = estimate, se = std.error))
  ) |>
  mutate(metan_summ = map(
    data, ~ rma(yi = estimate, vi = vi, data = .x, measure = "GEN", method = "REML") |>
      glance() |>
      select(tau2 = tau.squared, nobs = nobs)
  )) |>
  select(trt, term, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  group_by(trt) |>
  mutate(
    pred_lb = prop1-qt(.975, df = nobs)*sqrt(tau2+(se^2)),
    pred_ub = prop1+qt(.975, df = nobs)*sqrt(tau2+(se^2))
  ) |>
  select(trt, term, pred_lb, pred_ub, prop0, prop1, prop2, prop3) |>
  common_meta_format()


# Appendix table of relevant data ------------------------------------------
appdx_analysis2 <- dat |> 
  filter(
    term == "15-64 vs 65-74"   | 
      term == "75-84 vs 65-74"   |
      term == "85-99 vs 65-74"   |
      term == "Women vs Men"     
  ) |>
  select(term, jurisdiction, estimate, lower, upper) |> 
  mutate(
    estimate = exp(estimate),
    lower = exp(lower),
    upper = exp(upper)
  ) |>
  group_by(term) |>
  gt() |>
  gt_appdx_common() |>
  cols_align(align = "left", columns = c(jurisdiction)) |>
  fmt_number(
    columns = c(estimate, lower, upper),
    decimals = 2
  )  |>
  cols_label (
    jurisdiction = "Jurisdiction", 
    estimate = "Odds ratio",
    lower = "(95% confidence interval)"
  )

appdx_analysis2 |> saveRDS("results/appdx_analysis2.RDS")

# Create the meta-analysis table ------------------------------------------
table <- dat |>
  filter(!exclude) |>
  group_by(term) |>
  mutate(minimum = min(estimate)) |>
  mutate(maximum = max(estimate)) |>
  group_by(term, minimum, maximum) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(yi = estimate, vi = vi, data = .x, measure = "GEN", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop = estimate, lower = conf.low, upper = conf.high, se = std.error))
  ) |>
  mutate(metan_summ = map(
    data, ~ rma(yi = estimate, vi = vi, data = .x, measure = "GEN", method = "REML") |>
      glance() |>
      select(i2 = i.squared, tau2 = tau.squared, nobs = nobs)
  )) |>
  select(term, minimum, maximum, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  arrange(term) |>
  group_by() |>
  mutate(
    pred_lb = prop-qt(.975, df = nobs)*sqrt(tau2+(se^2)),
    pred_ub = prop+qt(.975, df = nobs)*sqrt(tau2+(se^2))
  ) |>
  mutate(
    prop  = exp(prop),
    lower = exp(lower),
    upper = exp(upper),
    pred_lb = exp(pred_lb),
    pred_ub = exp(pred_ub),
    tau = sqrt(tau2),
    minimum = exp(minimum),
    maximum = exp(maximum)
  ) |> 
  select(term, prop, lower, upper, pred_lb, pred_ub, tau, i2, minimum, maximum) 

chemo_analysis2_table <- table |>
  gt() |>
  gt_metan_common() |>
  fmt_number(
    columns = c(prop, lower, upper, pred_lb, pred_ub, minimum, maximum),
    decimals = 2
  ) |>
  cols_align(align = "left", columns = term) |>
  cols_label (
    term = "",
    prop = "Average odds ratio",
    lower = "(95% confidence interval)",
    pred_lb = "95% prediction interval *",
    minimum = "Observed jurisdictional range",
    tau = md("*Tau*, log-odds scale **"),
    i2 = md("*I<sup>2</sup>* ***")
  ) 

chemo_analysis2_table |> saveRDS("results/chemo_analysis2_table.RDS")



# Create the forest plot --------------------------------------------------

# xlabels for plot
xlabels <- c("1/64","1/32","1/16","1/8","1/4","1/2","1","2")

# limit to relevant terms and transform variables
dat <- dat |>
  filter(term == "Women vs Men" | term == "15-64 vs 65-74" | term == "75-84 vs 65-74" | term == "85-99 vs 65-74") |>
  mutate(
    lower = exp(lower),
    upper = exp(upper),
    estimate = exp(estimate)
  ) 

dat_stage <- dat_stage |>
  filter(term == "Women vs Men" | term == "15-64 vs 65-74" | term == "75-84 vs 65-74" | term == "85-99 vs 65-74") |>
  mutate(
    lower = exp(lower),
    upper = exp(upper),
    estimate = exp(estimate)
  ) 

poly_metan <- poly_metan |> filter(term == "Women vs Men" | term == "15-64 vs 65-74" | term == "75-84 vs 65-74" | term == "85-99 vs 65-74") |>
  mutate(
    x = exp(x),
    pred_lb = exp(pred_lb),
    pred_ub = exp(pred_ub)
  ) 

# All combined
p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(yvar = ordering)
p <- p + incl_point(yvar = ordering, xvar = estimate)
# Excluded jurisdictions and estimates
p <- p + excl_error(yvar = ordering)
p <- p + excl_point(yvar = ordering, xvar = estimate)
# Meta-analysis
p <- p + meta_error(y_shift = 19) 
p <- p + meta_diamond(y_shift = 19, groupvar = term)
### specific details for extreme low values (Newfoundland)
p <- p + geom_errorbarh(
  data = dat |> 
    filter(lower < (1/128)) |>
    mutate(lower = ifelse(lower <= 1/128, 1/128, lower)), 
  height=.1, 
  colour = "#dcdcdc",
  aes(y = ordering, xmin = lower, xmax = upper)
)
p <- p + geom_segment(
  data = dat |> filter(lower < (1/128)),
  aes(xend = 1/128, x = 1/128, yend = ordering, y = ordering),
  arrow = arrow(length = unit(0.1, "cm")),
  color = "#dcdcdc"
)
### Stage-specific results
p <- p +
  geom_errorbarh(
    data = dat_stage |>
      filter(!exclude),
    height = .1,
    aes(y = ordering + 0.2, xmin = lower, xmax = upper, colour = ctry_class)
  )
p <- p +
  geom_point(
    data = dat_stage |>
      filter(!exclude), 
    shape = 21,
    fill = "white",
    aes(y = ordering + 0.2, x = estimate, colour = ctry_class)
  )
# y axis
p <- p + jdn_yscale()
# x axis
p <- p + odds_xscale()
p <- p + geom_vline(xintercept = 1)
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
ggsave("results/chemo_odds.svg",
       plot = p,
       width = 15,
       height = 11.25,
       units = "cm")


# Clean up ----------------------------------------------------------------
rm(list = ls())


