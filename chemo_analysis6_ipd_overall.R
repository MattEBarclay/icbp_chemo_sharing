## Meta-analysis - dotplot of chemo/radio by country and cancer

library(metafor)
library(tidyverse)
library(broom)
library(gt)
library(lme4)

# Useful vectors for sorting things out
load("Rdata/useful_vectors.Rdata")

# Relevant data file
trt_unadj_n <- readRDS("Rdata/trt_unadj_n.RDS")

# Useful functions
source("chemo_analysis0_repeated_functions.R")

# Create dataset for meta-analysis ----------------------------------------

# Pull out age split for each cancer site
# We do not have an age/sex split
dat <- trt_unadj_n |>
  filter(
    stage == "All stages",
    ((variable == "age") & (cancer != "All 8 cancers")),
    jurisdiction != "",
    trt == "Chemo"
  ) |>
  mutate(age = case_when(
    variable == "age" & variable_value == "64" ~ "15-64",
    variable == "age" & variable_value == "74" ~ "65-74",
    variable == "age" & variable_value == "84" ~ "75-84",
    variable == "age" & variable_value == "99" ~ "85-99",
    TRUE ~ ""
  )) |>
  mutate(age = factor(age, levels = c("15-64","65-74","75-84","85-99"))) |>
  select(trt, jurisdiction, ordering, cancer, age, n, n_trt, prop, lower, upper) 

# apply common cleaning and exclusion rules
dat <- dat |>
  mutate(site_order  = factor(cancer, levels = site_chm)) |>
  chemo_dat_clean() |>
  chemo_exclusions_all() |>
  arrange(site_order, ordering)
  


# IPD meta-analysis -------------------------------------------------------
# Well,  "IPD" - adjusted for cancer and age

all_unadj   <- glmer(cbind(n_trt,n-n_trt) ~ (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!(jurisdiction %in% c("Manitoba", "Victoria", "Wales", "Newfoundland & Labrador", "Prince Edward Island", "Saskatchewan"))))

all_cas     <- glmer(cbind(n_trt,n-n_trt) ~ site_order + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!(jurisdiction %in% c("Manitoba", "Victoria", "Wales", "Newfoundland & Labrador", "Prince Edward Island", "Saskatchewan"))))

all_cas_age <- glmer(cbind(n_trt,n-n_trt) ~ site_order*age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!(jurisdiction %in% c("Manitoba", "Victoria", "Wales", "Newfoundland & Labrador", "Prince Edward Island", "Saskatchewan"))))

summary(all_unadj)
summary(all_cas)
summary(all_cas_age)

oes_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Oesophageal"))
oes_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Oesophageal"))

sto_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Stomach"))
sto_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Stomach"))

col_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Colon"))
col_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Colon"))
  
rec_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Rectal"))
rec_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Rectal"))
  
liv_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Liver"))
liv_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Liver"))
  
pan_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Pancreatic"))
pan_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Pancreatic"))

lun_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Lung"))
lun_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Lung"))
  
ova_all <- glmer(cbind(n_trt,n-n_trt) ~       (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Ovarian"))
ova_age <- glmer(cbind(n_trt,n-n_trt) ~ age + (1|jurisdiction), family = binomial(link = "logit"), data = dat |> filter(!exclude) |> filter(site_order == "Ovarian"))

# list RE standard deviations
VarCorr(all_unadj)
VarCorr(all_cas)
VarCorr(all_cas_age)
VarCorr(oes_all)
VarCorr(oes_age)
VarCorr(sto_all)
VarCorr(sto_age)
VarCorr(col_all)
VarCorr(col_age)
VarCorr(rec_all)
VarCorr(rec_age)
VarCorr(liv_all)
VarCorr(liv_age)
VarCorr(pan_all)
VarCorr(pan_age)
VarCorr(lun_all)
VarCorr(lun_age)
VarCorr(ova_all)
VarCorr(ova_age)


 
# Clean up ----------------------------------------------------------------
rm(list = ls())

