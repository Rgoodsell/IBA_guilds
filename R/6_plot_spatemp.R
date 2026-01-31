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
library(ggtern)

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
    theme_linedraw(base_size = 30)+
    theme(legend.position = "bottom" , 
          legend.key.width = unit(1 , "cm"),
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 22))+
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
         mean_w_prec    = mean_m_prec) |> 
  mutate(week_year     = 30 , 
         sampling_time = 10, 
         trap_id       = all_data$trap_id[1]) 

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
predData$feeding_niche <- factor(predData$feeding_niche , 
                                 levels = c("Phytophagous","Saprophagous","Predator",
                                            "Phytophage-parasitoid" , "Saprophage-parasitoid","Predator-parasitoid"),
                                 labels = c("Phytophage","Saprophage","Predator",
                                            "Phyto-parasitoid" , "Sapro-parasitoid","Pred-parasitoid"))

# make plots
sp_plots <- predData |> drop_na() |> 
  group_split(feeding_niche) |> 
  {\(.) map( . , ~ plot_fun(.,trap_data,sweden))}()

# build patchwork
p_geo <- (sp_plots[[1]] + sp_plots[[2]] + sp_plots[[3]]) /
  (sp_plots[[4]] + sp_plots[[5]] + sp_plots[[6]]) 


# Render figure
tiff(filename = "figures/guild_distributions.tiff" , width = 1500, height = 1500 , compression = "lzw")
p_geo 
dev.off()

# Check
browseURL("figures/guild_distributions.tiff")


# predict from best model to get change in phenology over time --------------------------------



# Build new data 
n_obs <- 2e3
t_trend <- get_env_s(raw_data$week_year , 5 , raw_data$mean_w_temp , n_obs = n_obs)

longitude_wgs84 <- mean(all_data$longitude_wgs84)
latitude_wgs84  <- mean(all_data$latitude_wgs84)
sampling_time   <- mean(all_data$sampling_time)
trap_id         <- all_data$trap_id[[1]]
feeding_niche   <- unique(all_data$feeding_niche)

newData  <- data.frame(week_year       =  seq(5,52,l=n_obs), 
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


# compositional data --------------------------------------------------------------------------

host_compositions <- predData |>
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) |> 
                      filter(host_parasitoid=="Host") |> 
                      select(week_year , feeding_niche,species_richness ,trap_id) |> 
                      mutate(rowid = row_number()) |> 
pivot_wider(names_from = "feeding_niche" , values_from = species_richness,values_fill = 0) |> 
  group_by(week_year ,trap_id) |> 
summarise(Phytophagous = sum(Phytophagous) , 
          Saprophagous = sum(Saprophagous) , 
          Predator = sum(Predator)) |> 
  rowwise() |> 
                      mutate(total_richness = sum(c_across(Phytophagous:Predator),na.rm=TRUE)) |> 
                      mutate(Phytophagous = Phytophagous / total_richness,
                             Saprophagous = Saprophagous / total_richness,
                             Predator = Predator / total_richness) %>% ungroup()

parasitoid_compositions <- predData |>
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) |> 
  filter(host_parasitoid=="Parasitoid") |> 
  select(week_year , feeding_niche,species_richness ,trap_id) |> 
  mutate(rowid = row_number()) |> 
  pivot_wider(names_from = "feeding_niche" , values_from = species_richness,values_fill = 0) |> 
  group_by(week_year ,trap_id) |> 
  summarise(Phytophage_parasitoid = sum(`Phytophage-parasitoid`) , 
            Saprophage_parasitoid = sum(`Saprophage-parasitoid`) , 
            Predator_parasitoid = sum(`Predator-parasitoid`)) |> 
  rowwise() |> 
  mutate(total_richness = sum(c_across(Phytophage_parasitoid:Predator_parasitoid),na.rm=TRUE)) |> 
  mutate(Phytophage_parasitoid = Phytophage_parasitoid / total_richness,
         Saprophage_parasitoid = Saprophage_parasitoid / total_richness,
         Predator_parasitoid = Predator_parasitoid / total_richness) %>% ungroup()


all_compositions <- predData |>
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) |> 
  select(week_year , feeding_niche,species_richness ,trap_id) |> 
  mutate(rowid = row_number()) |> 
  pivot_wider(names_from = "feeding_niche" , values_from = species_richness,values_fill = 0) |> 
  group_by(week_year ,trap_id) |> 
  
  summarise(      Phytophagous = sum(Phytophagous) , 
                  Saprophagous = sum(Saprophagous) , 
                  Predator     = sum(Predator),
    Phytophage_parasitoid      = sum(`Phytophage-parasitoid`) , 
         Saprophage_parasitoid = sum(`Saprophage-parasitoid`) , 
           Predator_parasitoid = sum(`Predator-parasitoid`)) |> 
  rowwise() |> 
  mutate(total_richness = sum(c_across(Phytophagous:Predator_parasitoid),na.rm=TRUE)) |> 
  mutate(Phytophagous = Phytophagous / total_richness,
         Saprophagous = Saprophagous / total_richness,
         Predator = Predator / total_richness,
    Phytophage_parasitoid = Phytophage_parasitoid / total_richness,
         Saprophage_parasitoid = Saprophage_parasitoid / total_richness,
         Predator_parasitoid = Predator_parasitoid / total_richness) %>% ungroup()

# Ternary plots
host_compositions %>% 
  ggtern(aes(x = Phytophagous,y=Saprophagous,z=Predator))+
  geom_point(aes(colour=week_year))

parasitoid_compositions |> 
  ggtern(aes(x = Phytophage_parasitoid,y=Saprophage_parasitoid,z=Predator_parasitoid))+
  geom_point(aes(colour=week_year))


# Summary bar charts
week_labs <- c(`10` = "Week 10" , `26`="Week 26" , `40`= "Week 40")

p_acomp <- all_compositions %>% 
  mutate(week_year = round(week_year,1)) %>% 
  filter(week_year %in% c(10,26,40)) %>% 
  pivot_longer(cols = 3:8 , names_to = "guild" , values_to = "proportion") %>% 
  summarise(proportion = mean(proportion), .by=c("week_year" , "guild")) %>% 
  mutate(host_parasitoid = case_when(str_detect(guild , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host")) %>% 
  mutate(host_pair  = str_extract(guild , "Phyto|Predator|Sapro")) %>% 
  mutate(host_parasitoid = fct_relevel(host_parasitoid , "Parasitoid" , "Host")) %>% 
  ggplot()+
  geom_bar(aes(host_pair , proportion, fill = host_parasitoid) , position="stack", stat="identity",colour="Black",linewidth=2)+
  facet_wrap(~week_year,labeller = as_labeller(week_labs))+
  theme_linedraw(base_size = 40)+
  theme(axis.text.x = element_blank(),legend.position = 'bottom')+
  labs(y = "Proportion" , x = "" , fill = "Trophic level")+
  scale_fill_viridis_d(option = "F" , end = .8,direction=-1)


p_hcomp <- host_compositions %>% 
  mutate(week_year = round(week_year,1)) %>% 
  filter(week_year %in% c(10,26,40)) %>% 
  pivot_longer(cols = 3:5 , names_to = "guild" , values_to = "proportion") %>% 
  summarise(proportion = mean(proportion), .by=c("week_year" , "guild")) %>% 
  ggplot()+
  geom_bar(aes(guild , proportion),fill="#03051AFF" , stat="identity",colour="Black" ,linewidth=2)+
  facet_wrap(~week_year,labeller = as_labeller(week_labs))+
  theme_linedraw(base_size = 40)+
  theme(axis.text.x = element_blank())+
  labs(y = "Proportion" , x = "")+
  scale_colour_viridis_d(option = "F" , end = .8,direction=-1)
  

p_pcomp <- parasitoid_compositions %>% 
  mutate(week_year = round(week_year,1)) %>% 
  filter(week_year %in% c(10,26,40)) %>% 
  pivot_longer(cols = 3:5 , names_to = "guild" , values_to = "proportion") %>% 
  summarise(proportion = mean(proportion), .by=c("week_year" , "guild")) %>% 
  mutate(guild = str_remove(guild , "_parasitoid")) %>% 
  ggplot()+
  geom_bar(aes(guild , proportion),fill="#F69C73FF" , stat="identity",colour="Black" ,linewidth=2)+
  facet_wrap(~week_year,labeller = as_labeller(week_labs))+
  theme_linedraw(base_size = 40)+
  theme(axis.text.x = element_text(angle = 90))+
  labs(y = "Proportion" , x = "Guild")+
  scale_fill_viridis_d(option = "F" , end = .8,direction=-1)
  


tiff("figures/guild_composition.tiff" , width = 1500 , height = 2000 , compression = "lzw")
p_acomp / p_hcomp / p_pcomp + plot_layout(guides='collect')+ plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')
dev.off()


browseURL("figures/guild_composition.tiff")


# plot ----------------------------------------------------------------------------------------

predData <- readRDS("data/sim_temp_phen.rds")
facet_labels <- c("Phyto" = "Phytophage", "Sapro" = "Saprophage", "Predator" = "Predator")
predData$host_parasitoid <- factor(predData$host_parasitoid , levels = c("Host" , "Parasitoid"))
peakData$host_parasitoid <- factor(peakData$host_parasitoid , levels = c("Parasitoid" , "Host"))
peakData$phase <- factor(peakData$phase , levels = c("start" , "peak", "end"), labels = c("Start","Peak","End"))
predData$host_pair <- factor(predData$host_pair , levels = c("Phyto","Sapro","Predator"))
peakData$host_pair <- factor(peakData$host_pair , levels = c("Phyto","Sapro","Predator"))



p_time <- predData |> ggplot(aes(week_year , species_richness , group = host_parasitoid))+
  geom_line(aes(colour=host_parasitoid),lwd = 2,show.legend = FALSE)+
  facet_wrap(~host_pair ,labeller = labeller(host_pair = facet_labels))+
  scale_colour_viridis_d(option = "F" , end = .8)+
  theme_linedraw(base_size = 40)+
  labs(x = "Week of the year" , y = "Species richness" , colour = "Level")


p_peaks <- peakData |>
  ggplot(aes(est , host_parasitoid))+
  geom_errorbar(aes(xmin = lwr , xmax = upr) , width = .1 , lwd=2)+
  geom_point(aes(colour=host_parasitoid,shape=phase),size = 8)+
  facet_grid(~host_pair , scales = "free",labeller = labeller(host_pair = facet_labels))+
  scale_colour_viridis_d(option = "F" , end = .8,direction=-1)+
  scale_y_discrete(expand = c(.5,.5))+
  scale_x_continuous(limits=c(10,50))+
  theme_linedraw(base_size = 40)+
  labs(x = "Week of the year" , y = "" , colour = "Level", shape = "Phase")+
  theme(legend.position = "bottom",aspect.ratio = .4)

tiff("figures/guild_phenology.tiff" , width = 1800 , height = 1250 , compression = "lzw")
p_time / p_peaks + plot_layout(guides = "collect",heights = c(1, .5)) + plot_annotation(tag_levels = 'A') &
  theme(legend.position = "bottom" , legend.box.margin = margin(-1,0,0,0 , unit="cm"))
dev.off()

 browseURL("figures/guild_phenology.tiff")

# extra -------------------------------------------------------------------

pcomp <- predData %>% 
   mutate(week_year = round(week_year,1)) %>% 
   filter(week_year %in% c(10,26,40)) %>% 
   summarise(mrichness = mean(species_richness), .by=c("week_year" , "feeding_niche")) %>% 
   arrange(week_year, feeding_niche)
   
   
 hcomp <- host_compositions %>% 
   mutate(week_year = round(week_year,1)) %>% 
   filter(week_year %in% c(10,26,40)) 
 