### Source the data loading file
source("../icbp_m9_trt_data/meta0_loading_data.R")

# Save the necessary/useful data files produced by the data loading script
trt_q |> saveRDS(file = "Rdata/trt_q.RDS")
trt_q_perc |> saveRDS(file = "Rdata/trt_q_perc.RDS")
trt_unadj_n |> saveRDS(file = "Rdata/trt_unadj_n.RDS")
trt_yn |> saveRDS(file = "Rdata/trt_yn.RDS")

save(
  country_order, 
  country_order_chm, 
  country_order_reml, 
  country_order_reml2,
  site,
  site_chm,
  term_list,
  file = "Rdata/useful_vectors.Rdata"
  )

rm(list = ls())
