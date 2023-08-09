## Meta-analysis - dotplot of chemo/radio by country and cancer

#
library(metafor)
library(tidyverse)
library(broom)
library(gt)

# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data files
trt_q_perc <- readRDS("Rdata/trt_q_perc.RDS")
trt_yn <- readRDS("Rdata/trt_yn.RDS")
trt_unadj_n <- readRDS("Rdata/trt_unadj_n.RDS")
trt_q <- readRDS("Rdata/trt_q.RDS")

# Useful functions
source("chemo_analysis0_repeated_functions.R")

# Step 1. What is the use of chemotherapy and radiotherapy overall across 
# jurisdictions and cancer sites? 
# Description of use of chemotherapy and radiotherapy at the jurisdiction level 
# (i.e. by jurisdiction, without split by personal-level predictors or stage), 
# for each of the 8 cancers and for all cancers. 

# Create dataset for meta-analysis ----------------------------------------
meta_use <- 
  trt_unadj_n |>
  chemo_exclusions_all() |>
  filter(!exclude) |>
  filter(
    stage == "All stages",
    variable == "all",
    jurisdiction != "",
    trt == "Chemo"
  ) |>
  select(trt, jurisdiction, ordering, cancer, n, n_trt, prop, lower, upper) |>
  group_by(cancer, trt) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop = estimate))
  ) |>
  select(trt, cancer, metan) |>
  unnest(cols = c(metan)) |>
  distinct() |>
  mutate(prop = expit(prop)) |>
  rename(use = prop)
  
dat_timely <- trt_q_perc |> 
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
  chemo_exclusions_all() |>
  filter(!exclude) |>
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
  mutate(site_order  = factor(cancer, levels = site)) |>
  mutate(ctry_order  = factor(jurisdiction, levels = country_order_chm)) |>
  mutate(trt_order   = factor(trt, levels = c("Chemo", "Radio"))) |>
  mutate(ordering    = as.numeric(ctry_order)) |>
  mutate(ordering    = ifelse(ordering >=  5, ordering+1, ordering)) |>
  mutate(ordering    = ifelse(ordering >= 14, ordering+1, ordering)) |>
  filter(n_trt > 0)

meta_timely <- 
  dat_timely |>
  group_by(cancer, trt) |>
  nest() |>
  mutate(
    metan = map(
      data, 
      qe_function1
    )
  ) |>
  select(trt, cancer, metan) |>
  unnest(cols = c(metan)) |>
  select(trt, cancer, prop) |>
  distinct() |>
  rename(time = prop)

dat <- inner_join(meta_use, meta_timely)
dat <- dat |> 
  filter(cancer != "All 8 cancers")

dat

p <- ggplot()
p <- p + geom_point(
  data = dat,
  aes(
    y=use, 
    x=time
  )
)
p <- p + geom_text(
  data = dat,
  hjust = 0,
  nudge_x = 3,
  aes(
    y=use, 
    x=time,
    label = cancer
  )
)
p <- p + scale_x_continuous(
  limits=c(0, 210), 
  breaks = c(0, 50, 100, 150, 200), 
  minor_breaks = NULL, 
  expand = c(0,0), 
  name = "Average days from diagnosis to treatment"
)
p <- p + scale_y_continuous(
  name = "Average chemotherapy usage (%)", 
  breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
  labels = c("0%", "10%", "20%", "30%", "40%", "50%", "60%"),
  minor_breaks = NULL
)
p <- p + 
  theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", hjust=0 )) +
  theme(axis.text.x=element_text(size=10, color="black", hjust=0.5)) +
  theme(axis.text.y=element_text(size=10, color="black")) +
  theme(axis.ticks.length = unit(0, "cm")) +
  theme(panel.spacing.x = unit(0.5, "lines")) + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
  theme(panel.spacing.y = unit(0.3,"line"))
p 
ggsave("results/chemo_correlation.png",
       plot = p,
       width = 15,
       height = 15,
       units = "cm")

cor.test(dat$use, dat$time)

# Clean up
rm(list = ls())

