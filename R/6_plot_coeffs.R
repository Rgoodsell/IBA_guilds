# plot results from model ---------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(mgcv)
library(viridis)
library(gratia)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)

set.seed(3)

source("R/functions.R")
# load data -----------------------------------------------------------------------------------

raw_data <- readRDS("data/tidydata/all_data.rds") |> drop_na()

all_data <- raw_data 

trap_data <- all_data |>
  group_by(longitude_wgs84 , latitude_wgs84 , trap_id) |>
  summarise(total_richness = sum(n_OTU))

mod <- readRDS("data/mod_pk.rds")

# plot maps -----------------------------------------------------------------------------------

# function to render ggplots independently 
# -------------------------------------------------------- #
plot_fun <- function(guild_data , trap_data , sweden){
  
  
  minr   <- 0
  maxr   <- round(max(guild_data$species_richness,na.rm=TRUE),-1)
  medr   <- mean(c(minr,maxr))
  breaks <- c(minr,medr,maxr)
  limits <- c(breaks[1] , breaks[3])
  
  
  ggplot()+
    geom_sf(data=sweden, alpha = 0)+
    geom_raster(data=guild_data , 
                aes(longitude_wgs84 , latitude_wgs84, 
                    fill = (species_richness)))+
    geom_point(data = trap_data ,
               aes(longitude_wgs84 , latitude_wgs84), 
               show.legend = FALSE)+
    scale_fill_viridis_c(option="turbo" , label = function(x) sprintf("%.1f", x) , 
                         breaks = breaks , limits= limits)+
    scale_x_continuous(breaks = c(14,18,22))+
    theme_linedraw(base_size = 25)+
    theme(legend.position = "bottom" , 
          legend.key.width = unit(1 , "cm"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))+
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))+
    facet_wrap(~feeding_niche) + 
    coord_sf()+
    labs(x = "Longitude" , y ="Latitude" , fill = "Predicted richness")
  
}
# -------------------------------------------------------- #
# Get outline of sweden
sweden <- ne_countries(country = "sweden" , scale = "large" , returnclass = "sf")

# Set up new data
newData <- readRDS("data/tidydata/env_pred_layer.rds") |>
  mutate(mean_w_temp    = mean_m_temp , 
         mean_w_prec    = mean_m_prec,
         mean_w_gdd5    = max_m_gdd5 , 
         mean_w_cmprec = max_m_prec) |> 
  mutate(week_year     = 30 , 
         sampling_time = 10, 
         trap_id       = all_data$trap_id[1], 
         enf_forest    = forest_cover , dbf_forest =0) 

# Replicate for each guild            
predData <- replicate(6,newData , simplify = FALSE) |>
  bind_rows() |>
  mutate(feeding_niche = rep(unique(all_data$feeding_niche) , each = n()/6))  


# predict and plot
modsmooths <- smooths(mod)
ex         <- grep("trap_id",smooths(mod)) 
predData$species_richness <- predict(mod, predData , type="response" , exclude = c(modsmooths[ex])) 

saveRDS(predData , "data/sim_dist_sr.rds")
library(beepr) ; beep(5)
# plot ----------------------------------------------------------------------------------------

predData <- readRDS("data/sim_dist_sr.rds")

# make plots
sp_plots <- predData |> drop_na() |> 
  group_split(feeding_niche) |> 
  {\(.) map( . , ~ plot_fun(.,trap_data,sweden))}()

# build patchwork
p_geo <- (sp_plots[[1]] + sp_plots[[2]] + sp_plots[[3]]) /
  (sp_plots[[4]] + sp_plots[[5]] + sp_plots[[6]]) 


# Render figure
tiff(filename = "figures/guild_distributions.tiff" , width = 1500 , height = 1500 , compression = "lzw")
p_geo
dev.off()

# Check
browseURL("figures/guild_distributions.tiff")


# predict from best model to get change in phenology over time --------------------------------

# Get model
mod <- readRDS("data/mod_pk.rds")

# Build new data 
n_obs <- 2e3
t_trend <- get_env_s(raw_data$week_year , 5 , raw_data$mean_w_temp , n_obs = n_obs)

longitude_wgs84 <- mean(all_data$longitude_wgs84)
latitude_wgs84  <- mean(all_data$latitude_wgs84)
sampling_time   <- mean(all_data$sampling_time)
trap_id         <- all_data$trap_id[[1]]
feeding_niche   <- unique(all_data$feeding_niche)

newData  <- data.frame(week_year       =  seq(10,52,l=n_obs), 
                    mean_w_temp     = t_trend,
                    mean_w_prec     = 0,
                    longitude_wgs84 = longitude_wgs84,
                    latitude_wgs84  = latitude_wgs84 , 
                    trap_id         = trap_id,
                    forest_cover    = 0,
                    crop_cover      = 0, 
                    grass_cover     = 0,
                    shrub_cover     = 0,
                    sampling_time   = sampling_time) 
  
  predData <- replicate(6,newData , simplify = FALSE) |> 
    bind_rows() |>
    mutate(feeding_niche = rep(feeding_niche , each = n()/6)) 
  
# Predict
predData$species_richness <- predict(mod , newdata = predData , exclude = "s(trap_id)" , type = "response")
predData <- predData |> 
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) 

saveRDS(predData , "data/sim_temp_phen.rds")


# simulate from MAP  to get max SR estimates for each guild * trend ---------------------------
guilds <- unique(all_data$feeding_niche)
peakList <- list()
tList    <- list()

# Loop over guilds to estimate peaks of each guild * trend
for(i in seq_along(guilds)){
    gData      <- filter(predData , feeding_niche == guilds[i])
    peakList[[i]] <- simulate_MAP_max(gData , mod , n_sim = 100 , qint = c(0.1 , 0.9)) |> 
                  bind_rows() |>
                  mutate(feeding_niche = guilds[i])
}


peakData <- bind_rows(peakList) |> 
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) 

# plot ----------------------------------------------------------------------------------------

predData <- readRDS("data/sim_temp_phen.rds")
facet_labels <- c("Phyto" = "Phytophage", "Sapro" = "Saprophage", "Predator" = "Predator")
predData$host_parasitoid <- factor(predData$host_parasitoid , levels = c("Host" , "Parasitoid"))
peakData$host_parasitoid <- factor(peakData$host_parasitoid , levels = c("Parasitoid" , "Host"))


p_time <- predData |> ggplot(aes(week_year , species_richness , group = host_parasitoid))+
  geom_line(aes(colour=host_parasitoid),lwd = 2,show.legend = FALSE)+
  facet_wrap(~host_pair ,labeller = labeller(host_pair = facet_labels))+
  #scale_x_continuous(limits = c(20,40))+
  scale_colour_viridis_d(option = "F" , end = .8)+
  theme_linedraw(base_size = 40)+
  labs(x = "Week of the year" , y = "Species richness" , colour = "Level")


p_peaks <- peakData |>
  ggplot(aes(peak , host_parasitoid))+
  geom_errorbar(aes(xmin = lwr , xmax = upr) , width = .1 , lwd=2)+
  geom_point(aes(colour=host_parasitoid),size = 8)+
  facet_wrap(~host_pair , scales = "free",labeller = labeller(host_pair = facet_labels))+
  scale_colour_viridis_d(option = "F" , end = .8,direction=-1)+
  scale_y_discrete(expand = c(.5,.5))+
  scale_x_continuous(limits=c(27,33))+
  theme_linedraw(base_size = 40)+
  labs(x = "Week of the year" , y = "" , colour = "Level")+
  theme(legend.position = "bottom",aspect.ratio = .4)

tiff("figures/guild_phenology.tiff" , width = 2000 , height = 1500 , compression = "lzw")
p_time / p_peaks + plot_layout(guides = "collect",heights = c(1, .5)) + plot_annotation(tag_levels = 'A') &
  theme(legend.position = "bottom" , legend.box.margin = margin(-1,0,0,0 , unit="cm"))
dev.off()

browseURL("figures/guild_phenology.tiff")
