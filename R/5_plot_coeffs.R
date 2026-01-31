# plot coefficients from model ----------------------------------------------------------------

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

# Environmental data
all_data <- readRDS("data/tidydata/all_data.rds") 
mod      <- readRDS("data/mod_pk.rds")

# options for plotting
bs <- 40 # base size
ylimsS <-c(-2,1.1) # smooth terms ylims


# plot coefficients ---------------------------------------------------------------------------

# Check distributions for sensible values for simulating marginal effects
plot(density(all_data$mean_w_gdd5, na.rm = TRUE))
plot(density(all_data$mean_w_cmprec , na.rm = TRUE))
plot(density(all_data$mean_w_temp , na.rm = TRUE))
plot(density(all_data$mean_w_prec , na.rm = TRUE))
# Seasonality ---------------------------------------------------------------------------------


newDataS <- expand.grid(week_year        = 5:52 , 
                        mean_w_temp       = mean(all_data$mean_w_temp , na.rm=TRUE) , 
                        mean_w_gdd5       = mean(all_data$mean_w_gdd5, na.rm=TRUE),
                        mean_w_prec       = mean(all_data$mean_w_prec,na.rm=TRUE) ,
                        mean_w_cmprec     = mean(all_data$mean_w_cmprec,na.rm=TRUE),
                        forest_cover      = 10 , 
                        grass_cover       = 10 , 
                        shrub_cover       = 10 ,
                        crop_cover        = 10 , 
                        feeding_niche     =  unique(all_data$feeding_niche), 
                        trap_id           = all_data$trap_id[1], 
                        latitude_wgs84    = 64.1 , 
                        longitude_wgs84   = 20.2 , 
                        sampling_time     = 10 ) |> 
           mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
                  host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) 




# Global and group level predictions

lp           <- predict(mod , newdata = newDataS,exclude = "s(trap_id,feeding_niche)" , type = "terms")
newDataS$lp   <- lp[,"s(week_year)"] + lp[,"s(week_year,feeding_niche)"]
newDataS$lpg  <- lp[,"s(week_year)"] 

facet_labels <- c("Phyto" = "Phytophage", "Sapro" = "Saprophage", "Predator" = "Predator")
newDataS$host_parasitoid <- factor(newDataS$host_parasitoid  , levels = c("Host" , "Parasitoid"))
newDataS$host_pair <- factor(newDataS$host_pair  , levels = c("Phyto","Sapro" , "Predator"))

pS <-  ggplot(newDataS ,
              aes(week_year , lp , 
                  colour = host_parasitoid, group=feeding_niche))+
  geom_line(lwd = 3 , lty=1)+
  scale_colour_viridis_d(option = "F" , end = .8)+
  theme_linedraw(base_size = bs)+
  labs(x = "Week of the year" , y = "Partial effect" , colour = "Feeding guild")+
  ggtitle(label = expression(f(S[g])))+
  scale_y_continuous(limits = ylimsS)+
  theme(plot.title = element_text(hjust = 0.5) , legend.position = "bottom")+
  facet_wrap(~host_pair,labeller=as_labeller(facet_labels))



# temperature ---------------------------------------------------------------------------------

newDataT <- expand.grid(week_year        = 25 , 
                        mean_w_temp       = seq(-5,30,l=100) , 
                        mean_w_gdd5       = mean(all_data$mean_w_gdd5, na.rm=TRUE),
                        mean_w_prec       = mean(all_data$mean_w_prec,na.rm=TRUE) ,
                        mean_w_cmprec     = mean(all_data$mean_w_cmprec,na.rm=TRUE),
                        forest_cover = 10,
                        grass_cover  = 10 , 
                        shrub_cover  = 10 ,
                        crop_cover   = 10 , 
                        feeding_niche   =  unique(all_data$feeding_niche), 
                        trap_id           = all_data$trap_id[1], 
                        latitude_wgs84    = 64.1 , 
                        longitude_wgs84   = 20.2 , 
                        sampling_time     = 10 )|> 
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) 



# Global and group level predictions
lp           <- predict(mod , newdata = newDataT, exclude = "s(trap_id,feeding_niche)" , type = "terms")
newDataT$lp   <- lp[,"s(mean_w_temp)"] + lp[,"s(mean_w_temp,feeding_niche)"]
newDataT$lpg  <- lp[,"s(mean_w_temp)"] 

facet_labels <- c("Phyto" = "Phytophage", "Sapro" = "Saprophage", "Predator" = "Predator")
newDataT$host_parasitoid <- factor(newDataT$host_parasitoid  , levels = c("Host" , "Parasitoid"))
newDataT$host_pair <- factor(newDataT$host_pair  , levels  = c("Phyto","Sapro" , "Predator"))

pT <-  ggplot(newDataT ,
              aes(mean_w_temp , lp , 
                  colour = host_parasitoid, group=feeding_niche))+
  geom_line(lwd = 3 , lty=1)+
  scale_colour_viridis_d(option = "F" , end = .8)+
  theme_linedraw(base_size = bs)+
  scale_y_continuous(limits = ylimsS)+
  labs(x = "Average weekly temperature" , y = "Partial effect" , colour = "Feeding guild")+
  ggtitle(label = expression(f(T[g])))+
  theme(plot.title = element_text(hjust = 0.5) , legend.position = "bottom")+
  facet_wrap(~host_pair,labeller=as_labeller(facet_labels))


# precipitation -------------------------------------------------------------------------------

newDataP <- expand.grid(week_year        = 25 , 
                        mean_w_temp       = mean(all_data$mean_w_temp,na.rm=TRUE) ,    
                        mean_w_gdd5       = mean(all_data$mean_w_gdd5, na.rm=TRUE),
                        mean_w_prec       = seq(0 , 10 , l = 100) ,
                        mean_w_cmprec     = mean(all_data$mean_w_cmprec,na.rm=TRUE),
                        forest_cover = 10 , 
                        grass_cover  = 10 , 
                        shrub_cover  = 10 ,
                        crop_cover   = 10 , 
                        feeding_niche   =  unique(all_data$feeding_niche), 
                        trap_id           = all_data$trap_id[1], 
                        latitude_wgs84    = 64.1 , 
                        longitude_wgs84   = 20.2 , 
                        sampling_time     = 10 ) |> 
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro")) 



# Global and group level predictions
lp           <- predict(mod , newdata = newDataP,exclude = "s(trap_id,feeding_niche)" , type = "terms")
newDataP$lp   <- lp[,"s(mean_w_prec)"] + lp[,"s(mean_w_prec,feeding_niche)"]
newDataP$lpg  <- lp[,"s(mean_w_prec)"] 


facet_labels <- c("Phyto" = "Phytophage", "Sapro" = "Saprophage", "Predator" = "Predator")
newDataP$host_parasitoid <- factor(newDataP$host_parasitoid  , levels = c("Host" , "Parasitoid"))
newDataP$host_pair <- factor(newDataP$host_pair  , levels = c("Phyto","Sapro" , "Predator"))


pP <- ggplot(newDataP ,
             aes(mean_w_prec , lp , 
                 colour = host_parasitoid, group=feeding_niche))+
  geom_line(lwd = 3 , lty=1)+
  scale_colour_viridis_d(option = "F" , end = .8)+
  theme_linedraw(base_size = bs)+
  scale_y_continuous(limits = ylimsS)+
  labs(x = "Average weekly precipitation" , y = "Partial effect" , colour = "Feeding guild")+
  ggtitle(label = expression(f(P[g])))+
  theme(plot.title = element_text(hjust = 0.5) , legend.position = "bottom")+
  facet_wrap(~host_pair,labeller=as_labeller(facet_labels))


# habitat coefficients ------------------------------------------------------------------------
newDataH <- expand.grid(week_year        = 25 , 
                        mean_w_temp       = mean(all_data$mean_w_temp , na.rm=TRUE) , 
                        mean_w_gdd5       = mean(all_data$mean_w_gdd5, na.rm=TRUE),
                        mean_w_prec       = mean(all_data$mean_w_prec,na.rm=TRUE) ,
                        mean_w_cmprec     = mean(all_data$mean_w_cmprec,na.rm=TRUE),
                        forest_cover = 10 , 
                        grass_cover  = 10 , 
                        shrub_cover  = 10 ,
                        crop_cover   = 10 , 
                        feeding_niche   =  unique(all_data$feeding_niche), 
                        trap_id           = all_data$trap_id[1], 
                        latitude_wgs84    = 64.1 , 
                        longitude_wgs84   = 20.2 , 
                        sampling_time     = 10 )

habP   <- predict(mod , newdata = newDataH,exclude = "s(trap_id,feeding_niche)" , type = "terms" , se.fit=TRUE)
lph    <- habP$fit[,1:4]
lph.se <- habP$se.fit[,1:4]
colnames(lph) <- paste0(c("forest_cover" , "crop_cover" , "shrub_cover" , "grass_cover"),".val")
colnames(lph.se) <- paste0(c("forest_cover" , "crop_cover" , "shrub_cover" , "grass_cover"),".se")

habData <- cbind(lph , lph.se) |> as.data.frame() |> 
  transform(feeding_niche = newDataH$feeding_niche) |> 
  pivot_longer(-feeding_niche,
               names_to = c("habitat", ".value"), 
               names_sep="\\." ) |> 
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host"),
         host_pair       = str_extract(feeding_niche , "Phyto|Predator|Sapro"),
         habitat = str_remove(habitat, "_cover$")) 


habData$feeding_niche <- factor(habData$feeding_niche , 
                                 levels = c("Phytophagous","Saprophagous","Predator",
                                            "Phytophage-parasitoid" , "Saprophage-parasitoid","Predator-parasitoid"),
                                 labels = c("Phytophage","Saprophage","Predator",
                                            "Phytophage-parasitoid" , "Saprophage-parasitoid","Predator-parasitoid"))


pH <- ggplot(habData, aes(habitat , val))+ 
  geom_hline(yintercept = 0 , lty = 3 , size = 2)+
  geom_errorbar(aes(ymin = val - se , ymax = val + se) , width = 0 , size = 2)+
  geom_point(aes(colour = host_parasitoid) , show.legend = TRUE , size = 7)+
  scale_colour_viridis_d(option = "F" , end = .8)+
  facet_wrap(~feeding_niche)+
  theme_linedraw(base_size = bs)+
  theme(axis.text.x =  element_text(angle = 90,vjust = 0.3), legend.position="bottom")+
  labs(x = "Habitat cover" , y = "Partial effect" , colour = "Level")

# sensitivity analysis for habitat variables --------------------------------------------------

hab <- all_data %>% select(forest_cover, crop_cover, grass_cover, shrub_cover) %>% drop_na()
pca <- prcomp(hab, scale. = TRUE)
summary(pca)

quant <- quantile(pca$x[,1],probs = seq(0,1,l=20))
xbreaks <- names(quant) %>% str_remove("%") %>% as.numeric() %>% round() 
scenarios <- map_dfr(quant, function(q) {
  all_data %>%
    slice(which.min(abs(pca$x[,1] - q))) %>%
    select(crop_cover, shrub_cover, forest_cover, grass_cover)
})

p1hs <- scenarios %>% 
  pivot_longer(everything()) %>% 
  mutate(sample = rep(1:20 , each = 4)) %>% 
  ggplot(aes(sample , value , group = name))+
  geom_line(aes(colour = name),lwd=1.5)+
  scale_colour_brewer(palette = 2, type = "qual")+
  theme_linedraw()+
  theme(legend.position = "none")+
  labs(y="Cover value" , x= "Quantile along PC1")+
  scale_x_continuous(breaks = xbreaks)


# choose some representative scenarios
scenarios <- scenarios %>% slice(c(1,3,5,18)) 

p2hs <- scenarios %>% 
  pivot_longer(everything()) %>% 
  mutate(scenario = 
           rep(factor(c("Grass-dominated" , 
                        "Agricultural" , 
                        "Mixed mosaic",  
                        "High forest"),
                      levels = c("Grass-dominated" , 
                                 "Agricultural" , 
                                 "Mixed mosaic",  
                                 "High forest")), each = 4)) %>% 
  ggplot()+
  geom_bar(aes(name, value,fill = name),stat = "identity")+
  facet_wrap(~scenario,ncol=4)+
  theme_linedraw()+
  scale_fill_brewer(palette = 2, type = "qual")+
  theme(legend.position = "bottom", axis.text = element_text(angle = 90))+
  labs(x = "Cover type" , y = "Proportion")


p1hs +p2hs  


# predictions -------------------------------------------------------------



longitude_wgs84 <- mean(all_data$longitude_wgs84,na.rm=TRUE)
latitude_wgs84  <- mean(all_data$latitude_wgs84,na.rm=TRUE)
trap_id         <- all_data$trap_id[[1]]
feeding_niche   <- unique(all_data$feeding_niche)
mean_w_temp     <- mean(all_data$mean_w_temp,na.rm=TRUE)
mean_w_prec     <- mean(all_data$mean_w_prec,na.rm=TRUE)



newData  <- data.frame(week_year       =  30, 
                       mean_w_temp     = mean_w_temp ,
                       mean_w_prec     = mean_w_prec,
                       longitude_wgs84 = longitude_wgs84,
                       latitude_wgs84  = latitude_wgs84 , 
                       trap_id         = trap_id,
                       forest_cover    = scenarios$forest_cover[],
                       crop_cover      = scenarios$crop_cover, 
                       grass_cover     = scenarios$grass_cover,
                       shrub_cover     = scenarios$shrub_cover,
                       sampling_time   = 7,
                       scenario = factor(c("GD" , "AM" , "MM",  "HF"),
                                         levels = c("GD" , "AM" , "MM",  "HF"))) 

# Replicate for each guild            
predData <- replicate(6,newData , simplify = FALSE) |>
  bind_rows() |>
  mutate(feeding_niche = rep(unique(all_data$feeding_niche) , each = n()/6))  

# predict and plot
modsmooths <- smooths(mod)
ex         <- grep("trap_id",smooths(mod)) 
sr_preds <- predict(mod, predData , type="response" , exclude = c(modsmooths[ex]),se.fit = TRUE) 
predData$species_richness <- sr_preds$fit
predData$ci_fit <- sr_preds$se.fit

pHab_sens <- predData %>% 
  group_by(feeding_niche) %>%  
  mutate(host_parasitoid = case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , TRUE ~ "Host")) %>% 
  ggplot(aes(scenario , species_richness))+
  geom_errorbar(aes(ymin = species_richness - ci_fit , ymax = species_richness+ci_fit), lwd = 2, width = 0.1)+
  geom_point(pch = 21 , aes(fill = host_parasitoid),size = 7)+
  facet_wrap(~feeding_niche ,scales="free")+
  theme_linedraw(base_size = bs)+
  theme(axis.text.x  = element_text(angle = 90))+
  theme(axis.text.x =element_text(angle = 90), legend.position = 'bottom')+
  labs(y = "Proportion" , x = "" , fill = "Guild")+
  scale_fill_viridis_d(option = "F" , end = .8,direction=1)


# plot --------------------------------------------------------------------


tiff("figures/habitat_coeffs.tiff", width = 1200 , height = 1000 , compression = "lzw")
pH + plot_layout(guides = "collect") & theme(legend.position = 'bottom')   
dev.off()

browseURL("figures/habitat_coeffs.tiff")


tiff("figures/habitat_sens.tiff", width = 1500 , height = 1000 , compression = "lzw")
pHab_sens  & theme(legend.position = 'bottom')   
dev.off()


browseURL("figures/habitat_sens.tiff")

tiff("figures/smooth_coefficient_plots.tiff" , width = 2000 , height = 2000 , compression = "lzw")
 (pT / pP / pS) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')      
dev.off()


browseURL("figures/smooth_coefficient_plots.tiff")


