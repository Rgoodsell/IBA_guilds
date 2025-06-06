# fit hierarchical gams  ------------------------------------------------------------------------------------

library(tidyverse)
library(mgcv)

set.seed(3)

# load data -----------------------------------------------------------------------------------

# scale data
all_data <- readRDS("data/tidydata/all_data.rds")

# fit gam ------------------------------------------------------------------------------------

knots <- list(week_year = c(0, 52))
ctrl <- list(nthreads=8)

mod <- bam(n_OTU ~ 
             
             s(feeding_niche , forest_cover , bs = "re")
           + s(feeding_niche , crop_cover   , bs = "re")
           + s(feeding_niche , shrub_cover  , bs = "re")
           + s(feeding_niche , grass_cover  , bs = "re")
           
           + s(feeding_niche , bs = "re")
           
           + s(week_year , bs = "cc" , k = 6 , m=2)  
           + s(week_year ,feeding_niche , bs = "fs", xt=list(bs="cc") , k=6, m = 2)
           
           + s(mean_w_temp , bs = "tp" ,k = 6, m = 2)
           + s(mean_w_temp , feeding_niche , bs = "fs", k=6 , m = 2)
           
           + s(mean_w_prec , bs = "tp" ,k = 6, m = 2)
           + s(mean_w_prec , feeding_niche , bs = "fs", k=6 , m = 2)
           
           + s(longitude_wgs84 , latitude_wgs84 , bs = "tp")
           + s(longitude_wgs84 , latitude_wgs84 , feeding_niche, bs = "fs")
           
           + s(trap_id , feeding_niche, bs = "re")
           + offset(log(sampling_time)),
           knots = knots,
           family = "nb" , method = "REML",  data = all_data , select = TRUE)

saveRDS(mod , "data/mod_pk.rds")
