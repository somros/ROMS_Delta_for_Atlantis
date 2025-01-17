# Alberto Rovellini
# The monthly approach will not work here, we need to go to the finest time step
# Which, unfortunately, happens to be 12-hours
# So, ymonmean and ydaymean from CDO are no help here
# For maximum flexibility, use a hybrid approach of CDO and R
# In CDO: merge the annual NETCDF for historical and hindcast into one nc file
# in R, read in the netcdf file and transform it into a dataframe (it WILL be large)
# get time step index (be careful with leap years)
# Group by box, layer, julian day, and summarize: mean and sd for hindcast and historical run
# turn to array
# read in projection files one at a time
# apply deltas to projection

library(tidyverse)
library(ncdf4)
library(tidync)
library(lubridate)

# merge hindcast files
# shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime data\\hindcast\\temp\\*.nc data\\hindcast\\temp\\hindcast_merged.nc")

# merge the historical
# shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime data\\historical\\temp\\*.nc data\\historical\\temp\\historical_merged.nc")

# what variable are we handling?
this_variable <- "salinity"

# files
hindcast_file <- "data/hindcast/salt/hindcast_merged.nc"
historical_file <- "data/historical/salt/historical_merged.nc"

# function that produces the delta arrays
# there are 4: mean and sd hindcast and historical run per time step for the reference period
# need two versions of each: one with and one without leap years
# produce one for all data including leap years, then drop julian day 60 to make a set of arrays for non-leap years (else array operations will not work on all projection files)
# make sure that the 732 time steps refer to the correct days (i.e., we have issues with midnight of Jan 1st missing)
# how does that propagate?

# function
make_delta_array <- function(variable, sim_period = "hindcast", mean=TRUE, leap=FALSE){
  
  # variable = "temperature"
  # sim_period = "projection"
  # mean=TRUE
  # leap=FALSE
  
  # set file
  if(sim_period == "hindcast"){
    this_file <- hindcast_file
  } else if(sim_period == "historical") {
    this_file <- historical_file
  } else if(sim_period == "projection") {
    this_file <- projection_file
  } else {
    stop("Check simulation period")
  }
  
  # pull variables with tidync
  salt_nc <- tidync(this_file) 
  
  # list variables
  these_vars <- hyper_grids(salt_nc) %>% # all available grids in the ROMS ncdf
    pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
    purrr::map_df(function(x){
      salt_nc %>% activate(x) %>% hyper_vars() %>% 
        mutate(grd=x)
    })
  
  # pull time origin
  salt_con <- nc_open(this_file)
  time_units <- ncatt_get(salt_con, "t", "units")$value # get time
  nc_close(salt_con)
  
  file_origin <- as.POSIXct(gsub("seconds since ", "", time_units))
  
  # pull variable
  grids <- these_vars %>% filter(name==variable) %>% pluck('grd')
  
  dat_salt <- salt_nc %>% activate(grids) %>% hyper_tibble()
  
  # numbering of b and z starts from 1 - change
  # and replace time steps with dates
  dat_salt <- dat_salt %>% 
    mutate(t = t + file_origin) %>%
    rename(b = y, z = x)
  
  # handle time
  # the time dimension in the transfromed ROMS forcings has a lot of obnoxious properties
  # The time series begins at 12:00 of Jan 1st
  # leap years cause problems, remove Feb 29 from the data
  dat_salt <- dat_salt %>%
    mutate(year = year(t),
           month = month(t),
           day = day(t),
           tod  = format(t, "%H:%M:%S")) 
  
  # add a first time step for the correct indexing of the data
  filler <- dat_salt %>% 
    filter(year == 1991, month == 1, day == 1) %>%
    mutate(t = t - hours(12),
           tod = format(t, "%H:%M:%S"))
  
  dat_salt <- rbind(filler, dat_salt) %>%
    filter(t < max(t)) # remove very last time step
  
  # remove Feb 29
  dat_salt <- dat_salt %>%
    filter(!(month == 2 & day == 29))
  
  # group by year, arrange by date, assign a time step
  dat_salt <- dat_salt %>%
    group_by(year, b, z) %>%
    mutate(ts = row_number()) %>%
    ungroup()
  
  # Group by box, layer, time step, and summarize (mean and sd as separate data frames)
  if(mean){
    aggregated_df <- dat_salt %>%
      group_by(b, z, ts) %>%
      summarise(!!sym(variable) := mean(!!sym(variable))) 
  } else {
    aggregated_df <- dat_salt %>%
      group_by(b, z, ts) %>%
      summarise(!!sym(variable) := sd(!!sym(variable))) 
  }
  
  # for the purpose of applying the correct delta to the projections, 
  # the delta for Jan 1 00:00 (jday = 1, tod = 00:00) must be the last time step, because annual files start at noon of Jan 1 and end at midnight of Jan 1 of the following year
  aggregated_df <- aggregated_df %>%
    mutate(ts = ifelse(ts == 1, (max(ts)+1), ts)) %>% # put the first time step last
    mutate(ts = ts - 1) %>%
    arrange(ts, b, z) %>%
    dplyr::select(ts, b, z, all_of(this_variable))
  
  # now turn the data frame into an array with dimensions z,b,t
  # the correct dimensioning is important because we will use element-wise operations for the delta correction of the projection files
  # first restore missing levels
  all_cells <- expand.grid(ts = 1:max(aggregated_df$ts), 
                           b = 1:max(aggregated_df$b), 
                           z = 1:max(aggregated_df$z)) # 109 boxes, 7 layers
  
  aggregated_df_complete <- all_cells %>%
    left_join(aggregated_df) 
  
  # for leap years, duplicate Feb 28
  
  if(leap){
    tmp_1 <- aggregated_df_complete %>%
      filter(ts < 118)
    
    tmp_2 <- aggregated_df_complete %>%
      filter(ts %in% c(116,117)) %>%
      mutate(ts = ts + 2) # add 2 time steps for Feb 29
    
    tmp_3 <- aggregated_df_complete %>%
      filter(ts >= 118) %>%
      mutate(ts = ts + 2) # add 2 time steps 
    
    aggregated_df_complete <- rbind(tmp_1, tmp_2, tmp_3)
    
  }
  
  # there is an annoying blip at the last 2 time steps
  # It seems to be mostly due to the hindcast
  # I do not fully understand why it is there, but to smooth it out use the delta from before the blip
  # moving the first time step to the end of the time series may be part of the cause, but this does not explain the second-to-last point
  # Future revisions should address it, but this is a close enough approximation for the purpose of delta-correcting
  tmp_1 <- aggregated_df_complete %>%
    filter(ts < 728)
  tmp_2 <- aggregated_df_complete %>%
    filter(ts == 727)
  
  if(leap){reps <- 5} else {reps <- 3} # number of time steps to pad depending on whether this is a leap year
  tmp_4 <- data.frame()
  for(i in 1:reps) {
    tmp_3 <- tmp_2
    tmp_3$ts <- 727 + i
    tmp_4 <- rbind(tmp_4, tmp_3)
  }
  
  aggregated_df_complete <- rbind(tmp_1, tmp_4)
  
  # Reshape the data frame to a wider format for each 'ts'
  ts_list <- aggregated_df_complete %>%
    group_by(ts) %>%
    pivot_wider(names_from = b, values_from = this_variable) %>%
    ungroup() %>%
    split(.$ts)
  
  # Convert each wide data frame to a matrix and stack them
  arr <- array(dim = c(length(unique(aggregated_df_complete$z)), 
                       length(unique(aggregated_df_complete$b)), 
                       length(unique(aggregated_df_complete$ts))))
  
  for (i in seq_along(ts_list)) {
    arr[,,i] <- as.matrix(ts_list[[i]][, -c(1,2)])  # Exclude the 'ts' and 'Z' columns
  }
  
  return(arr)
  
}

# apply the function to get mean and sd of the reference period
# TODO: set up a table with all combinations of conditions and use purrr::map() or similar

mean_hind <- make_delta_array(variable = this_variable,
                              sim_period = "hindcast",
                              mean = T,
                              leap = F)

# sd_hind <- make_delta_array(variable = this_variable,
#                               sim_period = "hindcast",
#                               mean = F,
#                               leap = F)

mean_hist <- make_delta_array(variable = this_variable,
                              sim_period = "historical",
                              mean = T,
                              leap = F)

# sd_hist <- make_delta_array(variable = this_variable,
#                               sim_period = "historical",
#                               mean = F,
#                               leap = F)

# for leap years
mean_hind_leap <- make_delta_array(variable = this_variable,
                                   sim_period = "hindcast",
                                   mean = T,
                                   leap = T)

# sd_hind_leap <- make_delta_array(variable = this_variable,
#                             sim_period = "hindcast",
#                             mean = F,
#                             leap = T)

mean_hist_leap <- make_delta_array(variable = this_variable,
                                   sim_period = "historical",
                                   mean = T,
                                   leap = T)

# sd_hist_leap <- make_delta_array(variable = this_variable,
#                             sim_period = "historical",
#                             mean = F,
#                             leap = T)

########################################################################################

# now loop over the projection files, and add the deltas with R's array operations
# This approach will be time-step agnostic
# be careful with NA's
# this works on annual forcing files
# In ACLIM, the formula for the bias correction is:
# X_proj' = X_hind + ((sigma_hind/sigma_hist) * (X_proj - X_hist))
# while the formula for the delta correction is
# X_proj' = X_hind + ((1/1) * (X_proj - X_hist))
# Scaling by the ratio of the variances does not work, not on 12-hourly files (some large differences in varainces can make the scalar very big)

proj_files <- list.files("data/ssp585/salt/projection/", full.names = T)

for(i in 1:length(proj_files)){
  
  # point at pre-correction file
  original_proj_file <- proj_files[i]
  # just file name to re-write later
  fn <- gsub("data/ssp585/salt/projection/","",original_proj_file)
  
  # open pre-correction file
  original_proj_nc <- nc_open(original_proj_file)
  original_proj_values <- ncvar_get(original_proj_nc, this_variable)
  
  # calculate values for the delta-corrected file
  if(dim(original_proj_values)[3]==730){
    
    # corrected_proj_values <- mean_hind + ((sd_hind/sd_hist) * (original_proj_values - mean_hist))
    corrected_proj_values <- mean_hind + ((1/1) * (original_proj_values - mean_hist))
    
  } else if (dim(original_proj_values)[3]==732) {
    
    # corrected_proj_values <- mean_hind_leap + ((sd_hind_leap/sd_hist_leap) * (original_proj_values - mean_hist_leap))
    corrected_proj_values <- mean_hind_leap + ((1/1) * (original_proj_values - mean_hist_leap))
    
  } else {
    
    stop("The projection nc array is not 730 or 732") # this would indicate a bug upstream in the translation
    
  }
  
  # now re-pack this to nc forcing files in a different folder (bias-corrected)
  # create a new nc file with the dimensions, variables, and attributes of the original one
  # the easiest way to do this is to copy the original file to a new location, open it, modify it, and close it
  corrected_proj_file <- paste0("data/ssp585/salt/projection_corrected/",fn)
  
  # create the folder
  if(!dir.exists("data/ssp585/salt/projection_corrected/")){
    dir.create("data/ssp585/salt/projection_corrected/")
  }
  
  # copy the original necdf file that we are correcting
  file.copy(original_proj_file, corrected_proj_file, recursive = T)
  
  # open the new file
  corrected_proj_nc <- nc_open(corrected_proj_file, write = T)
  
  # write the delta-corrected variable to it
  ncvar_put(corrected_proj_nc, varid = this_variable, vals = corrected_proj_values)
  
  # close both nc files
  nc_close(original_proj_nc)
  nc_close(corrected_proj_nc)
  
}
