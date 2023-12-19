

# Assign jurisdictions to country groups ----------------------------------
# Identifies the country for each jurisdiction
assign_ctry_class <- function(data) {
    
    data |> mutate(
      ctry_class = case_when(
        ctry_order == "Norway"           ~ 1,
        ctry_order == "England"          ~ 2,
        ctry_order == "Northern Ireland" ~ 2,
        ctry_order == "Scotland"         ~ 2,
        ctry_order == "Wales"            ~ 2,
        ctry_order == "Alberta"          ~ 3,
        ctry_order == "British Columbia" ~ 3,
        ctry_order == "Ontario"          ~ 3,
        ctry_order == "Saskatchewan"     ~ 3,
        ctry_order == "Manitoba"         ~ 3,
        ctry_order == "Prince Edward Island"  ~ 3,
        ctry_order == "New Brunswick"    ~ 3,
        ctry_order == "Newfoundland & Labrador" ~ 3,
        ctry_order == "Nova Scotia"      ~ 3,
        ctry_order == "New South Wales"  ~ 4,
        ctry_order == "Victoria"         ~ 4
      )
    ) |>
    mutate(
      ctry_class = factor(
        ctry_class, 
        levels = c(1, 2, 3, 4),
        labels = c("Norway", "United Kingdom", "Canada", "Australia")
      )
    )
  
}


# Identify exclusions (overall treatment data) ----------------------------
# Applies exclusion rules for the overall use and time-to-event data
chemo_exclusions_all <- function(data, jurisdiction = jurisdiction, cancer = cancer) {
  
  data |>
    filter(
      # Staging is different, remove entirely
      !(
        ({{jurisdiction}} == "Norway")  & 
        ({{cancer}} == "Colon, stage III")
      )
      &
      # filter out due to major data concerns
      # Note: missing in source data, this is just a confirmation step.
      !(
          ({{jurisdiction}} == "Prince Edward Island")  & 
          (
            {{cancer}} == "Ovarian" |
            {{cancer}} == "Liver"
          )
      )
    ) |>
    mutate(
      # exclude (but still present) where data may not be comparable
      # or specific sites with oral chemo issues
      exclude = case_when(
        {{jurisdiction}} == "Newfoundland & Labrador" ~ 1,
        (
          ({{jurisdiction}} == "Prince Edward Island")  & 
          (
            {{cancer}} %in% c("All 8 cancers", "Age 15-64", "Age 65-74", "Age 75-84", "Age 85-99", "Men", "Women")
          )
        ) ~ 1,
        (
          ({{jurisdiction}} == "Saskatchewan")  & 
            (
              {{cancer}} %in% c("All 8 cancers", "Age 15-64", "Age 65-74", "Age 75-84", "Age 85-99", "Men", "Women")
            )
        ) ~ 1,
        (
          ({{jurisdiction}} %in% c("Manitoba", "Victoria", "Wales", "Norway"))  & 
            (
              {{cancer}} %in% c("Colon, stage III", "Colon", "Rectal", "Liver")
            )
        ) ~ 1,
        # Specific problem with chemo ovarian data
        (
          ({{jurisdiction}} == "Saskatchewan")  & 
            (
              {{cancer}} == "Ovarian" 
            )
        ) ~ 1,
        # filter out all cancers-combined due to oral chemo issue
        (
          ({{jurisdiction}} %in% c("Manitoba", "Victoria", "Wales", "Norway"))  & 
            (
              {{cancer}} %in% c("All 8 cancers", "Age 15-64", "Age 65-74", "Age 75-84", "Age 85-99", "Men", "Women")
            )
        ) ~ 1,
        TRUE ~ 0
      )
    ) 
  
}


# Identify exclusions (comparisons data) ----------------------------------
# Applies exclusion rules for the comparisons data
chemo_exclusions_comparisons <- function(data) {
  
  data |>
    mutate(
      exclude = case_when(
        jurisdiction == "Newfoundland & Labrador" ~ 1,
        TRUE ~ 0
      )
    ) 
  
}


# Standard data cleaning tasks --------------------------------------------
# applies data cleaning tasks that are consistent across the different analyses
chemo_dat_clean <- function(data) {
  
  data |>
    mutate(trt_order   = factor(trt, levels = c("Chemo", "Radio"))) |>
    mutate(ctry_order  = factor(jurisdiction, levels = country_order_chm)) |>
    mutate(ordering    = as.numeric(ctry_order)) |>
    mutate(ordering    = ifelse(ordering >=  6, ordering+1, ordering)) |>
    mutate(ordering    = ifelse(ordering >= 15, ordering+1, ordering)) |>
    tibble() |> 
    assign_ctry_class() 
  
}


# expit function ----------------------------------------------------------
# inverse logistic function for use with proportions data
expit <- function(x) {
  
  x = exp(x)/(1+exp(x))
  
}


# Meta-analysis data formatting -------------------------------------------
# formatting of the poly_metan dataset that is consistent across analyses
common_meta_format <- function(data) {
  
  data |>
    pivot_longer(cols = starts_with("prop"), names_to = "prop", values_to = "x") |> 
    mutate(y = ifelse(prop == "prop0" | prop == "prop2", 0, ifelse(prop == "prop1", 0.3, ifelse(prop == "prop3", -0.3, NA)))) |>
    mutate(trt_order   = factor(trt, levels = c("Chemo", "Radio"))) 
  
}


# Chemotherapy very common plot options -----------------------------------
theme_icbp <- function() {
  
  theme(
    strip.text.x=element_text(size=8, color="black", hjust=0 ), 
    axis.text.x=element_text(size=6, color="black"),
    axis.text.y=element_text(size=5, color="black"),
    axis.ticks.length = unit(0, "cm"),
    panel.spacing.x = unit(0.5, "lines"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    panel.spacing.y = unit(0.3,"line") ,
    legend.position = "none"
  ) 

}

# Colour scheme
icbp_colours <- c("#524A6F", "#003566", "#18ACBD", "#BC1C80")

icbp_colour_manual <- function() {
  scale_color_manual(values = icbp_colours)
}

icbp_fill_manual <- function() {
  scale_fill_manual(values = icbp_colours)
}

# Standard y-axis for jurisdictions
jdn_yscale <- function() {
  scale_y_continuous(
    name = "", 
    #limits = c(0. 20), 
    breaks = c(1:5, 7:14, 16:17, 19), 
    minor_breaks = NULL, 
    labels = country_order_chm, 
    #labels = blank_order,
    trans="reverse"
  ) 
}

# y-axis for cancer sites
site_yscale_use <- function() {
  scale_y_continuous(
    name = "",  
    ylim(0, 10),  
    breaks=1:10, 
    minor_breaks = NULL, 
    labels = c("All 8 cancers", cancer_chemo_order),  
    trans="reverse"
  ) 
}

site_yscale_time <- function() {
  scale_y_continuous(
    name = "",  
    ylim(0, 10),  
    breaks=1:9, 
    minor_breaks = NULL, 
    labels = site,  
    trans="reverse"
  ) 
}

# x-axis for proportions
prop_xscale <- function() {
  scale_x_continuous(
    limits = c(-0.02,1), 
    breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
    labels = c("0%", "20%", "40%", "60%", "80%"), 
    minor_breaks = NULL, 
    expand = c(0,0), 
    name = "Receiving chemotherapy (%)"
  )
}

# x-axis for odds
odds_xscale <- function() {
  scale_x_continuous(
    name = "Odds ratio", 
    limits = c(1/128, 2.5), 
    breaks = c(1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1, 2), 
    labels = xlabels, 
    minor_breaks = NULL,
    trans = "log",
    expand = expansion(mult = 0, add = 0)
  ) 
}

# x-axis for days
days_xscale <- function() {
  scale_x_continuous(
    limits=c(-10, 250), 
    breaks = c(0, 50, 100, 150, 200, 250), 
    labels = c("0", "50", "100", "150", "200", ""), 
    minor_breaks = NULL, 
    expand = c(0,0), 
    name = "Days from diagnosis to treatment (median and IQR)"
  )
}

daydiff_xscale <- function() {
  scale_x_continuous(
    name = "Difference in median days-to-treatment",  
    limits = c(-33.2, 15), 
    breaks = c(-28, -14, 0, 14), 
    minor_breaks = NULL
  )
}

# Included Jdns -----------------------------------------------------------
# Plot jurisdiction-specific results for the jdns used in the meta-analysis
incl_error <- function(yvar) {
  geom_errorbarh(
    data = dat |>
      filter(!exclude),
    height = .1,
    aes(y = {{yvar}}, xmin = lower, xmax = upper, colour = ctry_class)
  )
}
incl_point <- function(yvar, xvar) {
  geom_point(
    data = dat |>
      filter(!exclude), 
    aes(y = {{yvar}}, x = {{xvar}}, colour = ctry_class)
  )
}


# Excluded jdns -----------------------------------------------------------
# Plot jurisdiction-specific results for the jdns NOT used in the meta-analysis
excl_error <- function(yvar) {
  geom_errorbarh(
    data = dat |>
      filter( exclude != 0),
    height=.5,
    colour = "#dcdcdc",
    aes(y = {{yvar}}, xmin = lower, xmax = upper)
  )
}
excl_point <- function(yvar, xvar) {
  geom_point(
    data = dat |>
      filter( exclude != 0),
    shape = 21,
    colour = "#dcdcdc",
    fill = "white",
    aes(y = {{yvar}}, x = {{xvar}})
  )
}


# Meta analysis diamond ---------------------------------------------------
# Plot the meta-analysis summaries
meta_error <- function(y_shift) {
  geom_errorbarh(
    data = poly_metan |> 
      filter( y == 0) |>
      mutate(y = y + {{y_shift}}) |> 
      mutate(
        jurisdiction = "RE meta-analysis",
        ctry_order = factor(jurisdiction, levels = country_order_chm)
      ),
    height=.3,
    colour = "darkgrey",
    linetype = "solid",
    aes(y = y, xmin = pred_lb, xmax = pred_ub)
  ) 
}

meta_diamond <- function(y_shift, groupvar) {
  geom_polygon(
    data = poly_metan |>
      mutate(y = y + {{y_shift}}) |>
      mutate(
        jurisdiction = "RE meta-analysis",
        ctry_order = factor(jurisdiction, levels = country_order_chm)
      ),
    colour = "darkgrey",
    fill = "darkgrey",
    aes(
      x = x, 
      y = y,
      group = trt,
      subgroup = {{groupvar}}
    )
  )
}



# GT appendix common ------------------------------------------------------
# Common options for appendix tables
gt_appdx_common <- function(tabdata) {
  tabdata |>
  cols_merge(
    columns = c(lower, upper),
    pattern = "({1}, {2})"
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
}


# GT meta-analysis common -------------------------------------------------
# Common options for meta-analysis tables
gt_metan_common <- function(tabdata) {
  tabdata |>
    cols_merge(
      columns = c(lower, upper),
      pattern = "({1}, {2})"
    ) |>
    cols_merge(
      columns = c(pred_lb, pred_ub),
      pattern = md("{1}-{2}")
    ) |>
    cols_merge(
      columns = c(minimum, maximum),
      pattern = md("{1}-{2}")
    ) |>
    fmt_number(
      columns = c(i2),
      decimals = 1
    ) |>
    fmt_number(
      columns = c(tau),
      decimals = 3
    ) |>
    cols_align(align = "center") |>
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
}

# Metamedian wrapper functions --------------------------------------------
# These allow us to 'map()' metamedian, which otherwise throws a fit
# because it relies on vectors :/

# Function for CIs etc
qe_function1 <- function(x) {
  min = x$pct_0
  q1  = x$pct_25
  med = x$pct_50
  q3  = x$pct_75
  max = x$pct_100
  n   = x$n_trt
  
  qe(
    #min.g1 = min,
    q1.g1  = q1 ,
    med.g1 = med,
    q3.g1  = q3 ,
    #max.g1 = max,
    n.g1   = n
  ) |>
    tidy(conf.int = TRUE) |>
    select(prop = estimate, lower = conf.low, upper = conf.high, se = std.error)
}

# Function for summary stats
qe_function2 <- function(x) {
  min = x$pct_0
  q1  = x$pct_25
  med = x$pct_50
  q3  = x$pct_75
  max = x$pct_100
  n   = x$n_trt
  
  qe(
    #min.g1 = min,
    q1.g1  = q1 ,
    med.g1 = med,
    q3.g1  = q3 ,
    #max.g1 = max,
    n.g1   = n
  ) |>
    glance() |>
    select(i2 = i.squared, tau2 = tau.squared, nobs = nobs)
}

