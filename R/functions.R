
get_env_s <- function(week_year,k,envar , n_obs){
  # all_data = all observations by week from study sites
  # K = knots
  # envar
  
  mod     <- gam( envar ~  s(week_year , bs = "cc" , k = 5))
  newdata <- data.frame(week_year = seq(min(week_year), max(week_year) , l = n_obs))
  trend   <- predict(mod , newdata = newdata , type = "response")
  }


# Function to build new data to simulate changes in environment
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
  opt      <- rep(NA, n_sim)
  pstart   <- rep(NA, n_sim)
  pend     <- rep(NA, n_sim)
  
  pred  <- vector(mode = "list" , length = n_sim)
  ilink <- family(m_obj)$linkinv 
  
  for (i in seq_len(n_sim)) { 
    pred[[i]]   <- ilink((Xp %*% Vbdist[i, ]) +  m.offset)  # save predicted curve 
    opt[i]      <- phen_dat$week_year[which.max(pred[[i]])] # save peak SR
    
    phix       <- min(pred[[i]]) + max(pred[[i]]) /100 * 5 # get 5% of maximum species richness as an indicator of the start / end of season phenostates
    phase_ix   <- phen_dat$week_year[which(pred_r<phix)] # 
    pstart[i]  <- max(phase_ix[phase_ix<(26)])
    pend[i]    <- min(phase_ix[phase_ix>(26)])
  }

  peak   <- median(opt)
  p_ci   <- quantile(opt, qint)
  start  <- median(pstart,na.rm = TRUE)
  s_ci   <- quantile(pstart, qint)
  end  <- median(pend,na.rm=TRUE)
  e_ci   <- quantile(pend, qint)
  
  
  
  # data.frame of different curves from simulated coefs
  # sim_df <- tibble(week_year = rep(phen_dat$week_year,times=n_sim), 
  #                      n_OTU = do.call(rbind,pred)[,1]) |>
  #           mutate(sim = rep(1:n_sim, each=nrow(phen_dat)))
  

  out <-  data.frame(peak_est = peak   , peak_upr = p_ci[2] , peak_lwr = p_ci[1],
                     start_est = start , start_upr = s_ci[2] , start_lwr = s_ci[1],
                     end_est = end     , end_upr = e_ci[2] , end_lwr = e_ci[1]) |> 
    pivot_longer(
      everything(),
      names_to = c("phase", ".value"),
      names_pattern = "(.*)_(.*)"
    ) 
  return(out)
}




