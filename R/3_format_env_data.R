# format environmental data -------------------------------------------------------------------

library(tidyverse)
library(janitor)
library(lubridate)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(terra)
library(data.table)

# load ddtplyr# load data -----------------------------------------------------------------------------------

tempmax <- stack("data/raw_data/DownscaledTmax2019.tif")
tempmin <- stack("data/raw_data/DownscaledTmin2019.tif")
prec <- stack("data/raw_data/DownscaledPrcp2019.tif")

hab_data  <- readRDS("data/raw_data/trap_habitat.rds")$km_1 |> 
             bind_rows() |> clean_names() %>%
             mutate(forest_cover = rowSums(dplyr::select(.,matches("forest")))) |> 
             dplyr::select(trap_id , 
                           forest_cover , crop_cover , shrub_cover , grass_cover) 

# SR data for combining & getting sample locations
niche_sr    <- readRDS("data/tidydata/niche_SR_data.rds") |> ungroup() |> dplyr::select(-trap_habitat)
sample_locs <- niche_sr |> dplyr::select(longitude_wgs84 , latitude_wgs84 , trap_id) |> distinct()

# mask for sweden
sweden   <- ne_countries(country = "Sweden" , scale = "large")


# format data -----------------------------------------------------------------------------------

# take mean values of 4 key habitat types - already extracted from raster from previous project
 trap_habitat <- hab_data |> 
                 group_by(trap_id) |> 
                 summarise(across(forest_cover:grass_cover , ~mean(.x , na.rm=TRUE)))
                 
# temperature ---------------------------------------------------------------------------------

# Extract daily tempmax measurements for each trap
tListMx <- list()
for(i in 1:365){
  tListMx[[i]]      <- try(raster::extract(tempmax[[i]] ,  y = sample_locs[,1:2] , buffer = 1e2 , df=TRUE)) |>
                      as.data.frame.matrix() 
  
  tListMx[[i]][,1]        <- sample_locs$trap_id
  tListMx[[i]]$date      <- as.Date(i,origin="2018-12-31")
  colnames(tListMx[[i]]) <- c("trap_id" , "tmax","date")
  print(i)
}

# Extract daily tempmin measurements for each trap
tListMn <- list()
for(i in 1:365){
  tListMn[[i]]      <- raster::extract(tempmin[[i]] ,  y = sample_locs[,1:2] , buffer = 1e2 , df=TRUE) |>
    as.data.frame.matrix() 
  
  tListMn[[i]][,1]        <- sample_locs$trap_id
  tListMn[[i]]$date      <- as.Date(i,origin="2018-12-31")
  colnames(tListMn[[i]]) <- c("trap_id" , "tmin","date")
  print(i)
}


# get weekly mean measurements.

tmax <- tListMx |> bind_rows()
tmin <- tListMn |> bind_rows() |> dplyr::select(tmin)


weekly_temp  <- bind_cols(tmax , tmin) |> 
                mutate(tmax = tmax / 100 , 
                       tmin = tmin / 100 ,
                       tmax = case_when(tmax < -300 ~ NA_real_ , TRUE ~ tmax) , 
                       tmin = case_when(tmin < -300 ~ NA_real_ , TRUE ~ tmax) , 
                       week_year = week(date)) |> 
                arrange(week_year) |> 
                rowwise() |> mutate(mean_temp = mean(c_across(matches("tmin|tmax")))) |> ungroup() |> 
                group_by(trap_id) |> 
                mutate(td = case_when(mean_temp <= 5 ~ 0 , TRUE ~ mean_temp),
                     gdd5 = cumsum(td)) |> 
               group_by(week_year , trap_id) |> 
               summarise(mean_w_temp = mean(mean_temp, na.rm = TRUE),
                         mean_w_gdd5 = mean(gdd5 , na.rm=TRUE)) |> ungroup() |> drop_na()

# precipitation -------------------------------------------------------------------------------


# Extract daily temp measurements for each trap
pList <- list()
for(i in 1:365){
        print(i)
        pList[[i]]      <- raster::extract(prec[[i]] ,  y = sample_locs[,1:2] , buffer = 1e2 , df=TRUE) |>
                           as.data.frame.matrix() 
        
        pList[[i]][,1]       <- sample_locs$trap_id
        pList[[i]]$date      <- as.Date(i,origin="2018-12-31")
        colnames(pList[[i]]) <- c("trap_id" , "prec","date")
}



# get weekly mean measurements + cumulative measurements
weekly_prec  <- pList |> bind_rows() |> 
        mutate(prec = prec / 100 , 
               prec = case_when(prec < -30 ~ NA_real_ , TRUE ~ prec) , 
               week_year = week(date)) |> 
        group_by(trap_id) |> 
        mutate(tp = case_when(prec <= 0 ~ 0 , TRUE ~ prec),
               cumulative_prec = cumsum(tp)) |> 
        group_by(trap_id , week_year) |> 
          summarise(mean_w_prec = mean(prec, na.rm = TRUE),
                    mean_w_cmprec = mean(cumulative_prec , na.rm=TRUE)) |> ungroup() |> drop_na()

# temp prediction layers ---------------------------------------------------------------------------


sListMx <- list()
# For some reason looping over bands quicker than doing this one the whole stack...
# Only extract July days
sday <- yday(as.Date("2019-07-01"))
eday <- yday(as.Date("2019-07-31"))
for(i in 1:eday){
        
        print(i)
        sweden_tmax_day  <- raster::crop(x = tempmax[[i]]  , y = sweden) %>% mask( . , sweden) 
        sweden_tmin_day  <- raster::crop(x = tempmin[[i]]  , y = sweden) %>% mask( . , sweden) 
        
        valsMx             <- getValues(sweden_tmax_day) / 100
        valsMx[valsMx < -300] <- NA
        valsMn             <- getValues(sweden_tmin_day) / 100
        valsMn[valsMn < -300] <- NA
        sweden_tmax_day  <- setValues(sweden_tmax_day , valsMx)
        sweden_tmin_day  <- setValues(sweden_tmin_day , valsMx)
        

        tmax_layer <- as.data.frame(sweden_tmax_day , xy=TRUE) 
        tmin_layer <- as.data.frame(sweden_tmin_day , xy=TRUE) 
        
        temp_prediction_layer <- cbind(tmax_layer , tmin_layer[,3]) 
        
        colnames(temp_prediction_layer)[3:4] <- c("tmax","tmin")
        colnames(temp_prediction_layer)[1:2] <- c("longitude_wgs84" , "latitude_wgs84")
        temp_prediction_layer$date <- as.Date(i,origin="2018-12-31")

        sListMx[[i]] <- temp_prediction_layer
        
}


# Use data.table for transforms
temp_layer <- sListMx |> bind_rows() |>  as.data.table() 

temp_layer_sum <- temp_layer |> 
                _[, td := ifelse(tmax <= 5,0, tmax)] |> 
                na.omit() |> 
                _[,gdd5 := cumsum(td),by=list(longitude_wgs84 , latitude_wgs84)] |> 
                _[date >= sday | date <= eday, ] |> 
                _[,mean_temp := rowMeans(.SD, na.rm = TRUE), .SDcols = c("tmin", "tmax")] |> 
                _[, .(mean_m_temp = mean(mean_temp, na.rm=TRUE),  max_m_gdd5 = max(gdd5)), by=list(longitude_wgs84 , latitude_wgs84)]

                
# Clear some space
rm(temp_layer)
rm(sListMx)
                
# prec prediction layers ----------------------------------------------------------------------------------

pList <- list()
# For some reason looping over bands quicker than doing this one the whole stack...
# Only extract July days
sday <- yday(as.Date("2019-07-01"))
eday <- yday(as.Date("2019-07-31"))
for(i in 1:eday){
        
        print(i)
        sweden_prec_day  <- raster::crop(x = prec[[i]]  , y = sweden) %>% mask( . , sweden) 
        vals             <- getValues(sweden_prec_day) / 100
        vals[vals < -300] <- NA
        sweden_prec_day  <- setValues(sweden_prec_day , vals)
        
        prec_prediction_layer <- as.data.frame(sweden_prec_day , xy=TRUE) 
        colnames(prec_prediction_layer)[3] <- "prec"
        colnames(prec_prediction_layer)[1:2] <- c("longitude_wgs84" , "latitude_wgs84")
        prec_prediction_layer$date <- as.Date(i,origin="2018-12-31")
        
        pList[[i]] <- prec_prediction_layer
        
}

# Use data.table for transforms
prec_layer <- pList |> bind_rows() |>  as.data.table() 

prec_layer_sum <- prec_layer |> 
                    na.omit() |> 
                  _[,cum_prec := cumsum(prec),by=list(longitude_wgs84 , latitude_wgs84)] |> 
                  _[date >= sday | date <= eday, ] |> 
                  _[, .(mean_m_prec = mean(prec,na.rm=TRUE),  max_m_prec = max(cum_prec)), by=list(longitude_wgs84 , latitude_wgs84)]



# Clear some space
rm(prec_layer)
rm(pList)


# get prediction layers ------------------------------------------------------------------------

# The Copernicus habitat layers need downloading from the copernicus web platform. 

# Forest layer
tree_cover    <- raster("PROBAV_LC100_global_v3.0.1_2018-conso_Tree-CoverFraction-layer_EPSG-4326.tif")
tree_cover    <- raster::crop(x = tree_cover  , y = sweden) %>% mask( . , sweden) 
cover1k       <- raster::aggregate(x=tree_cover , fact = 10 , fun = mean) 
forest        <- as.data.frame(cover1k , xy =TRUE) 

# crops
crop_cover    <- raster("PROBAV_LC100_global_v3.0.1_2018-conso_Crops-CoverFraction-layer_EPSG-4326.tif")
crop_cover    <- raster::crop(x = crop_cover  , y = sweden) %>% mask( . , sweden) 
crop1k        <- raster::aggregate(x=crop_cover , fact = 10 , fun = mean) 
crop          <- as.data.frame(crop1k) 

# grass
grass_cover    <- raster("PROBAV_LC100_global_v3.0.1_2018-conso_Grass-CoverFraction-layer_EPSG-4326.tif")
grass_cover    <- raster::crop(x = grass_cover  , y = sweden) %>% mask( . , sweden) 
grass1k        <- raster::aggregate(x=grass_cover , fact = 10 , fun = mean) 
grass         <- as.data.frame(grass1k) 

# shrub
shrub_cover    <- raster("PROBAV_LC100_global_v3.0.1_2018-conso_Shrub-CoverFraction-layer_EPSG-4326.tif")
shrub_cover    <- raster::crop(x = shrub_cover  , y = sweden) %>% mask( . , sweden) 
shrub1k         <- raster::aggregate(x=shrub_cover , fact = 10 , fun = mean) 
shrub           <- as.data.frame(shrub1k) 

# water bodies - used to mask out lakes in sweden
water_cover    <- raster("PROBAV_LC100_global_v3.0.1_2018-conso_PermanentWater-CoverFraction-layer_EPSG-4326.tif")
water_cover    <- raster::crop(x = water_cover  , y = sweden) %>% mask( . , sweden) 
water1k         <- raster::aggregate(x=water_cover , fact = 10 , fun = mean) 
water           <- as.data.frame(water1k) 

hab_pred_layer <- cbind(forest, crop, grass , shrub,water)  |> drop_na()
colnames(hab_pred_layer) <- c("longitude_wgs84" , "latitude_wgs84" , 
                              "forest_cover" , "crop_cover", "grass_cover" , "shrub_cover" , "water_cover")


# Replace values outside of training data with NAs to prevent extrapolation.
fRange <- range(trap_habitat$forest_cover , na.rm = TRUE)
cRange <- range(trap_habitat$crop_cover , na.rm = TRUE)
gRange <- range(trap_habitat$grass_cover , na.rm=TRUE)
sRange <- range(trap_habitat$shrub_cover , na.rm = TRUE)

hab_pred_layer[hab_pred_layer$forest_cover > fRange[2],"forest_cover"] <- NA
hab_pred_layer[hab_pred_layer$crop_cover   > cRange[2],"crop_cover"]   <- NA
hab_pred_layer[hab_pred_layer$grass_cover  > gRange[2],"forest_cover"] <- NA
hab_pred_layer[hab_pred_layer$shrub_cover  > sRange[2],"shrub_cover"]  <- NA




# bind layers ---------------------------------------------------------------------------------

hab_sweden <- hab_pred_layer %>% 
        st_as_sf(coords = c("longitude_wgs84" , "latitude_wgs84") , crs = crs(sweden))


env_sweden <- temp_layer_sum %>% 
              mutate(mean_m_prec = prec_layer_sum$mean_m_prec,
                     max_m_prec= prec_layer_sum$max_m_prec) |> 
              st_as_sf(coords = c("longitude_wgs84" , "latitude_wgs84") , crs = crs(sweden))



env_pred_layer <- env_sweden %>% drop_na() |> 
        st_join(hab_sweden , st_nearest_feature) %>% 
        mutate(longitude_wgs84 = st_coordinates(.)[,1] , 
               latitude_wgs84  = st_coordinates(.)[,2]) %>% 
        st_drop_geometry() |> 
        filter(water_cover < 100)
        

# save ----------------------------------------------------------------------------------------

saveRDS(weekly_temp , "data/tidydata/weekly_temp.rds") 
saveRDS(weekly_prec , "data/tidydata/weekly_prec.rds") 
saveRDS(trap_habitat , "data/tidydata/trap_habitat.rds")
saveRDS(hab_pred_layer , "data/tidydata/hab_pred_layer.rds")
saveRDS(env_pred_layer , "data/tidydata/env_pred_layer.rds")


# Join all data & save
all_data <- full_join(weekly_temp , weekly_prec , by = c("trap_id" , "week_year")) |> 
                full_join(trap_habitat) |> 
                full_join(niche_sr , by = c("trap_id" , "week_year")) |>
                mutate(trap_id = factor(trap_id)) |> 
                drop_na(n_OTU) |> droplevels()


saveRDS(all_data , "data/tidydata/all_data.rds")
saveRDS(all_data , "dardel/all_data.rds")

