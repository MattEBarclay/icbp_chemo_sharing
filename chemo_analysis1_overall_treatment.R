## Meta-analysis - dotplot of chemo/radio by country and cancer

library(metafor)
library(tidyverse)
library(broom)
library(gt)

# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data file
trt_unadj_n <- readRDS("Rdata/trt_unadj_n.RDS")

# Useful functions
source("chemo_analysis0_repeated_functions.R")

# Create dataset for meta-analysis ----------------------------------------
dat1 <- trt_unadj_n |>
  filter(
    stage == "All stages",
    (variable == "all") | ((variable == "age" | variable == "sex") & (cancer == "All 8 cancers")),
    jurisdiction != "",
    trt == "Chemo"
  ) |>
  mutate(cancer2 = case_when(
    variable == "age" & variable_value == "64" ~ "Age 15-64",
    variable == "age" & variable_value == "74" ~ "Age 65-74",
    variable == "age" & variable_value == "84" ~ "Age 75-84",
    variable == "age" & variable_value == "99" ~ "Age 85-99",
    variable == "sex" & variable_value == "M"  ~ "Men",
    variable == "sex" & variable_value == "F"  ~ "Women",
    variable == "sex" & variable_value == "m"  ~ "Men",
    variable == "sex" & variable_value == "f"  ~ "Women",
    TRUE ~ cancer
  )) |>
  mutate(cancer = cancer2) |>
  select(trt, jurisdiction, ordering, cancer, n, n_trt, prop, lower, upper) 

dat2 <- trt_unadj_n |>
  filter(
    (
        #variable_value == "2" | 
        variable_value == "3" | 
        #variable_value == "R" |  
        #variable_value == "II" | 
        variable_value == "III" | 
        #variable_value == "Regional spread (adjacent organs)" | 
        variable_value ==	"Regional spread (lymph nodes)"     
    ),
    (cancer == "Colon" & trt == "Chemo"),
    variable == "stage",
    stage == "All stages",
    jurisdiction != ""
  ) |>
  select(trt, jurisdiction, ordering, cancer, variable, n, n_trt) |>
  group_by(trt, jurisdiction, ordering, cancer, variable) |>
  summarise_all(sum) |>
  mutate(prop  = n_trt / n , 
         lower = (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2))),
         upper = (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2)))
  ) |>
  mutate(cancer = "Colon, stage III") |>
  select(trt, jurisdiction, ordering, cancer, n, n_trt, prop, lower, upper) 

# apply common cleaning and exclusion rules
site_chm2 <- c(site_chm, "Age 15-64", "Age 65-74", "Age 75-84", "Age 85-99", "Men", "Women")

dat <- bind_rows(dat1, dat2) |>
  mutate(site_order  = factor(cancer, levels = site_chm2)) |>
  chemo_dat_clean() |>
  chemo_exclusions_all() |>
  arrange(site_order, ordering)
  
rm(dat1, dat2, site_chm2)


# Appendix table of relevant data ------------------------------------------
appdx_analysis1 <- dat |> select(site_order, jurisdiction, n, n_trt, prop, lower, upper) |>
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
  
appdx_analysis1 |> saveRDS("results/appdx_analysis1.RDS")

# Meta-analysis data processing ------------------------------------------
poly_metan <- 
  dat |>
  filter(cancer %in% site_chm) |>
  filter(!exclude) |>
  group_by(cancer, trt) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop0 = conf.low, prop1 = estimate, prop2 = conf.high, prop3 = estimate, se = std.error))
  ) |>
  mutate(metan_summ = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      glance() |>
      select(tau2 = tau.squared, nobs = nobs)
  )) |>
  select(trt, cancer, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  group_by(trt) |>
  mutate(
    pred_lb = expit(prop1-qt(.975, df = nobs)*sqrt(tau2+(se^2))),
    pred_ub = expit(prop1+qt(.975, df = nobs)*sqrt(tau2+(se^2)))
  ) |>
  select(trt, cancer, pred_lb, pred_ub, prop0, prop1, prop2, prop3) |>
  mutate(site_order = factor(cancer, levels = site_chm)) |>
  common_meta_format() |>
  mutate(x = expit(x))


# Create the meta-analysis table --------------------------  ----------------
table <- dat |>
  filter(!exclude) |>
  group_by(site_order) |>
  mutate(minimum = min(prop)) |>
  mutate(maximum = max(prop)) |>
  group_by(site_order, minimum, maximum) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop = estimate, lower = conf.low, upper = conf.high, se = std.error))
  ) |>
  mutate(metan_summ = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      glance() |>
      select(i2 = i.squared, tau2 = tau.squared, nobs = nobs)
  )) |>
  select(site_order, minimum, maximum, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  mutate(
    pred_lb = expit(prop-qt(.975, df = nobs)*sqrt(tau2+(se^2))),
    pred_ub = expit(prop+qt(.975, df = nobs)*sqrt(tau2+(se^2)))
  ) |>
  mutate(
    prop  = expit(prop),
    lower = expit(lower),
    upper = expit(upper),
    tau = sqrt(tau2)
  ) |> 
  select(site_order, prop, lower, upper, pred_lb, pred_ub, tau, i2, minimum, maximum) |>
  group_by() |>
  arrange(site_order)

chemo_analysis1_table <- table |>
  gt() |>
  gt_metan_common() |>
  cols_align(align = "left", columns = site_order) |>
  fmt_percent(
    columns = c(prop, lower, upper, pred_lb, pred_ub, minimum, maximum),
    decimals = 1
  ) |>
  cols_label (
    site_order = "Cancer site",
    prop = "Average treatment %",
    lower = "(95% confidence interval)",
    pred_lb = "95% prediction interval *",
    minimum = "Observed jurisdictional range",
    tau = md("*Tau*, log-odds scale **"),
    i2 = md("*I<sup>2</sup>* ***")
  ) 

chemo_analysis1_table |> saveRDS("results/chemo_analysis1_table.RDS")



# Main version of meta-analysis plot --------------------------------------
dat <- dat |> filter(cancer %in% site_chm)

# Meta-analysis chemotherapy cancer-site ordering
chemo_ordering <- dat |>
  filter(!exclude) |>
  filter(cancer != "All 8 cancers") |>
  filter(trt == "Chemo") |>
  group_by(cancer) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop = estimate) |>
      mutate(prop = expit(prop))
  )
  ) |>
  select(cancer, metan) |>
  unnest(cols = c(metan)) |>
  arrange(prop)
chemo_ordering
cancer_chemo_order <- chemo_ordering$cancer

rm(chemo_ordering)

# Draw the forest plot
p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(yvar = ordering)
p <- p + incl_point(yvar = ordering, xvar = prop)
# Excluded jurisdictions and estimates
p <- p + excl_error(yvar = ordering)
p <- p + excl_point(yvar = ordering, xvar = prop)
# Meta-analysis
p <- p + meta_error(y_shift = 19) 
p <- p + meta_diamond(y_shift = 19, groupvar = site_order)
# y axis
p <- p + jdn_yscale()
# x axis
p <- p + prop_xscale()
# legend
p <- p + guides(
  alpha = "none",
  colour = "none"
)
# faceting
p <- p + facet_wrap(site_order~., ncol = 2, scales = "free_x")
# themes
p <- p + theme_bw() +
  theme_icbp() +
  icbp_colour_manual() 
p 
ggsave("results/chemo_overall_props.png",
       plot = p,
       width = 15,
       height = 22.5,
       units = "cm")


# Alternative layout of meta-analysis plot --------------------------------
dat <- dat |> mutate(
  site_order = factor(cancer, levels = c("All 8 cancers", cancer_chemo_order))
)
poly_metan <- poly_metan |> mutate(
  site_order = factor(cancer, levels = c("All 8 cancers", cancer_chemo_order))
)  

p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(yvar = as.numeric(site_order))
p <- p + incl_point(yvar = as.numeric(site_order), xvar = prop)
# Excluded jurisdictions and estimates
p <- p + excl_error(yvar = as.numeric(site_order))
p <- p + excl_point(yvar = as.numeric(site_order), xvar = prop)
# Meta-analysis
p <- p + meta_error(y_shift = as.numeric(site_order)) 
p <- p + meta_diamond(y_shift = as.numeric(site_order), groupvar = site_order)
# x axis
p <- p + prop_xscale()
# y axis
p <- p + site_yscale_use()
# faceting
p <- p + facet_wrap(ctry_order~., ncol = 4)
p <- p +  
  theme_bw() +
  theme_icbp() +
  icbp_colour_manual()
p <- p + theme(strip.text.x=element_text(size=7)) 
p 
ggsave("results/chemo_overall_props_alternative1.png",
       plot = p,
       width = 15,
       height = 13,
       units = "cm")

dat <- dat |> mutate(
  site_order = factor(cancer, levels = site)
)
poly_metan <- poly_metan |> mutate(
  site_order = factor(cancer, levels = site)
)  

# Line plot of treatment rates --------------------------------------------

trial <- dat |>
  filter(trt == "Chemo") |>
  filter(cancer != "All 8 cancers") |>
  select(jurisdiction, cancer, prop)
trial

j <- unique(trial$jurisdiction) 
c <- unique(trial$cancer)
extras <- expand.grid(j, c, stringsAsFactors = FALSE)
extras <- extras |> 
  rename(jurisdiction = Var1, cancer = Var2) |>
  mutate(ctry_order  = factor(jurisdiction, levels = country_order_reml2)) |>
  assign_ctry_class()

trial <- trial |>
  right_join(extras) |>
  mutate(site_order = factor(cancer, levels = cancer_chemo_order)) |>
  arrange(ctry_order, site_order) 

trial |> filter(jurisdiction == "Prince Edward Island")

p <- ggplot()
p <- p + geom_line(
  data = trial |> mutate(grouping = ctry_order) |> mutate(ctry_order = NULL), 
  colour = "gray", 
  aes(
    y = prop, 
    x = as.numeric(site_order),
    group = grouping
  )
)
p <- p + geom_line(
  data=trial, 
  size = 1,
  aes(
    y = prop, 
    x = as.numeric(site_order),
    group = ctry_order,
    colour = ctry_class
  )
)
p <- p + geom_point(
  data=trial, 
  size = 1,
  aes(
    y = prop, 
    x = as.numeric(site_order),
    group = ctry_order,
    colour = ctry_class
  )
)
p <- p + scale_x_continuous(
  name = "Cancer site",
  limits = c(1, 9),
  breaks = 1:9,
  minor_breaks = NULL,
  labels = cancer_chemo_order
)
p <- p + scale_y_continuous(
  limits = c(0, 1),
  breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
  labels = c("0%", "20%", "40%", "60%", "80%"), 
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
ggsave("results/chemo_overall_props_alternative2_lines.png",
       plot = p,
       width = 15,
       height = 13,
       units = "cm")

rm(trial)
rm(j, c, extras)


# Clean up ----------------------------------------------------------------
rm(list = ls())

