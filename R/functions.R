
get_env_s <- function(week_year,k,envar , n_obs){
  # all_data = all observations by week from study sites
  # K = knots
  # envar
  
  mod     <- gam( envar ~  s(week_year , bs = "cc" , k = 5))
  newdata <- data.frame(week_year = seq(min(week_year), max(week_year) , l = n_obs))
  trend   <- predict(mod , newdata = newdata , type = "response")
  }


# Function to build new data to simulate changes in environment
#! Not sure if the Gaussian approximation of the posterior is working well here when there are very few observations !
#! For some groups it ends up with the model making high SR estimates very early in the the year !
#! For now I have limited the prediction period to be from week 10 to avoid this when estimating SR peaks !
#! Ideally we should implement something more sensible !

#! As a quick and dirty option for some preliminary results, this uses a gam fit from our observations as the seasonal temp + prec components
#! Need to fix this so that we get a better observation of temperature to simulate from. Maybe fit from all data from COPERNICUS?

build_new_data <- function(fit_data, m_obj ,n_obs, t_temp , t_prec , delta_temp , delta_prec) {
  # fit_data   = data that model was fit to
  # m_obj      = fitted gam 
  # n_obs      = number of new observations to simulate
  # t_x        = current trend in environmental variables
  
  # delta_X    = average change to environmental variables
  
  
  # Extract seasonal trends of temp & precipitation from data 
  week_year  <- seq(min(fit_data$week_year), max(fit_data$week_year) , l = n_obs)
  mean_temp  <- t_temp + delta_temp # Modify average trends
  mean_prec  <- t_prec + delta_prec
  
  # Set data
  longitude_wgs84 <- mean(fit_data$longitude_wgs84)
  latitude_wgs84  <- mean(fit_data$latitude_wgs84)
  sampling_time   <- mean(fit_data$sampling_time)
  trap_id         <- fit_data$trap_id[[1]]
  trap_habitat    <- "forest"
  
  
  phen_dat <- data.frame( week_year       = week_year , 
                          mean_temp       = mean_temp ,
                          mean_prec       = mean_prec ,
                          longitude_wgs84 = longitude_wgs84,
                          latitude_wgs84  = latitude_wgs84 , 
                          trap_id = trap_id,
                          trap_habitat = trap_habitat,
                          sampling_time)
  
  return(phen_dat)
  
}

# Function to simulate from MAP posterior of fitted gam and calculate peak species richness
simulate_MAP_max <- function(phen_dat , m_obj , n_sim, qint = c(0.1, 0.9)){
  # phen_dat = new data
  # m_obj    = fitted gam object
  # n_sim    = number of simulations
  # qint     = length 2 vector of quantile interval
  
  pred_r   <- predict(m_obj, phen_dat, type = "response") # Predicted response
  Xp       <- predict(m_obj, phen_dat, type = "lpmatrix") # Coefficents
  m.offset <- attr(Xp , "model.offset")                   # Model offset
  
  beta     <- coef(m_obj) # Posterior mean and covariance of regression coefficients
  Vb       <- vcov(m_obj) 
  Vbdist    <- MASS::mvrnorm(n_sim, beta, Vb) # Random draws from joint distribution
  
  # Loop over mv distribution and calculate the maximum SR
  opt   <- rep(NA, n_sim)
  pred  <- vector(mode = "list" , length = n_sim)
  ilink <- family(m_obj)$linkinv 
  
  for (i in seq_len(n_sim)) { 
    pred[[i]]   <- ilink((Xp %*% Vbdist[i, ]) +  m.offset)  # save predicted curve 
    opt[i]      <- phen_dat$week_year[which.max(pred[[i]])] # save peak SR
  }
  
  peak   <- phen_dat$week_year[which.max(pred_r)]
  ci     <- quantile(opt, qint)
  
  # data.frame of different curves from simulated coefs
  # sim_df <- tibble(week_year = rep(phen_dat$week_year,times=n_sim), 
  #                      n_OTU = do.call(rbind,pred)[,1]) |>
  #           mutate(sim = rep(1:n_sim, each=nrow(phen_dat)))
  

  out <-  data.frame(peak = peak , upr = ci[2] , lwr = ci[1]) 
  return(out)
}




