## Meta-analysis - dotplot of chemo/radio by country and cancer

library(tidyverse)
library(broom)
library(gt)

# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data file
trt_unadj_n <- readRDS("Rdata/trt_unadj_n.RDS")

# Useful functions
source("chemo_analysis0_repeated_functions.R")

trt_unadj_n |> select(variable) |> unique()

# Create dataset for analysis ----------------------------------------
dat <- trt_unadj_n |>
  filter(
    stage == "All stages",
    (variable == "diagnosis_year"),
    jurisdiction != "",
    trt == "Chemo"
  ) |>
  mutate(diagnosis_year = variable_value) |>
  select(trt, jurisdiction, ordering, cancer, diagnosis_year, n, n_trt, prop, lower, upper) 

dat <- dat |>
  mutate(site_order  = factor(cancer, levels = site_chm)) |>
  chemo_dat_clean() |>
  chemo_exclusions_all() |>
  arrange(site_order, ordering)
  


# Appendix table of relevant data ------------------------------------------
appdx_analysis1 <- dat |> select(site_order, jurisdiction, diagnosis_year, n, n_trt, prop, lower, upper) |>
  group_by(site_order) |>
  filter(!is.na(n)) |>
  gt() |>
  gt_appdx_common() |>
  cols_align(align = "left", columns = c(jurisdiction)) |>
  fmt_number(
    columns = c(n, n_trt),
    decimals = 0
  ) |>
  fmt_percent(
    columns = c(prop, lower, upper),
    decimals = 1
  ) |>
  cols_label (
    jurisdiction = "Jurisdiction", 
    n = "Patients",
    n_trt = "Received chemotherapy",
    prop = "Average treatment %",
    lower = "(95% confidence interval)",
  ) 
  
appdx_analysis1 |> saveRDS("results/appdx_analysis5_trend.RDS")


# Line plot of treatment trends --------------------------------------------

trial <- dat |>
  filter(trt == "Chemo") |>
  filter(cancer == "All 8 cancers") |>
  select(jurisdiction, cancer, diagnosis_year, prop, lower, upper)
trial

j <- unique(trial$jurisdiction) 
c <- unique(trial$diagnosis_year)
extras <- expand.grid(j, c, stringsAsFactors = FALSE)
extras <- extras |> 
  rename(jurisdiction = Var1, diagnosis_year = Var2) |>
  mutate(ctry_order  = factor(jurisdiction, levels = country_order_reml2)) |>
  assign_ctry_class()

trial <- trial |>
  right_join(extras) |>
  arrange(ctry_order, diagnosis_year) 

trial |> filter(jurisdiction == "Prince Edward Island")
trial <- trial |> mutate(cancer = "All 8 cancers") |> mutate(diagnosis_year = as.numeric(diagnosis_year))
trial

p <- ggplot()
p <- p + geom_line(
  data = trial |> mutate(grouping = ctry_order) |> mutate(ctry_order = NULL), 
  colour = "gray", 
  aes(
    y = prop, 
    x = diagnosis_year,
    group = grouping
  )
)
p <- p + geom_line(
  data=trial, 
  size = 1,
  aes(
    y = prop, 
    x = diagnosis_year,
    group = ctry_order,
    colour = ctry_class
  )
)
p <- p + geom_line(
  data=trial, 
  size = .3,
  linetype = 2,
  aes(
    y = lower, 
    x = diagnosis_year,
    group = ctry_order,
    colour = ctry_class
  )
)
p <- p + geom_line(
  data=trial, 
  size = .3,
  linetype = 2,
  aes(
    y = upper, 
    x = diagnosis_year,
    group = ctry_order,
    colour = ctry_class
  )
)
p <- p + geom_point(
  data=trial, 
  size = 1,
  aes(
    y = prop, 
    x = diagnosis_year,
    group = ctry_order,
    colour = ctry_class
  )
)
p <- p + scale_x_continuous(
  name = "Diagnosis year",
  limits = c(2012, 2017),
  breaks = 2012:2017,
  minor_breaks = NULL,
  labels = c("2012","13","14","15","16","17")
)
p <- p + scale_y_continuous(
  limits = c(0, 0.5),
  breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), 
  labels = c("0%", "10%", "20%", "30%", "40%", "50%"), 
  minor_breaks = NULL, 
  expand = c(0,0), 
  name = "Receiving chemotherapy (%)"
)
p <- p +
  theme_bw() +
  theme_icbp() +
  icbp_colour_manual() +
  guides(
    alpha = "none",
    size = "none",
    fill = "none",
    colour = "none",
    shape = "none"
  ) 
p <- p + theme(axis.text.x=element_text(size=8, color="black", vjust = 0.5, hjust=1, angle = 90)) 
p <- p + facet_wrap(ctry_order~., ncol = 4)
p <- p + theme(strip.text.x=element_text(size=7)) 
p 
ggsave("results/figure6_chemo_overall_trend_appendix.svg",
       plot = p,
       width = 15,
       height = 13,
       units = "cm")
ggsave("results/figure6_chemo_overall_trend_appendix.pdf",
       plot = p,
       width = 15,
       height = 13,
       units = "cm")

rm(trial)
rm(j, c, extras)


# Clean up ----------------------------------------------------------------
rm(list = ls())

