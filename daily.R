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
# shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime data\\hindcast\\*.nc data\\hindcast\\hindcast_merged.nc")

# merge the historical
# shell("\"C:\\Users\\Alberto Rovellini\\CDO\\cdo.exe\" mergetime data\\historical\\*.nc data\\historical\\historical_merged.nc")

# what variable are we handling?
this_variable <- "temperature"

# files
hindcast_file <- "data/hindcast/hindcast_merged.nc"
historical_file <- "data/historical/historical_merged.nc"

# function that produces the delta arrays
# there are 4: mean and sd hindcast and historical run per time step for the reference period
# need two versions of each: one with and one without leap years
# produce one for all data including leap years, then drop julian day 60 to make a set of arrays for non-leap years (else array operations will not work on all projection files)
# make sure that the 732 time steps refer to the correct days (i.e., we have issues with midnight of Jan 1st missing)
# how does that propagate?

# function
make_delta_array <- function(variable, hindcast=TRUE, mean=TRUE, leap=FALSE){
  
  # set file
  if(hindcast){
    this_file <- hindcast_file
  } else {
    this_file <- historical_file
  }
  
  # pull variables with tidync
  temp_nc <- tidync(this_file) 

  # list variables
  these_vars <- hyper_grids(temp_nc) %>% # all available grids in the ROMS ncdf
    pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
    purrr::map_df(function(x){
      temp_nc %>% activate(x) %>% hyper_vars() %>% 
        mutate(grd=x)
    })
  
  # pull time origin
  temp_con <- nc_open(this_file)
  time_units <- ncatt_get(temp_con, "t", "units")$value # get time
  nc_close(temp_con)
  
  file_origin <- as.POSIXct(gsub("seconds since ", "", time_units))
  
  # pull variable
  grids <- these_vars %>% filter(name==variable) %>% pluck('grd')
  
  dat_temp <- temp_nc %>% activate(grids) %>% hyper_tibble()

  # numbering of b and z starts from 1 - change
  # and replace time steps with dates
  dat_temp <- dat_temp %>% 
    mutate(t = t + file_origin) %>%
    rename(b = y, z = x)
  
  # handle time: go from date to julian day and time of the day (12:00 vs 00:00)
  # daylight saving causes issues with grouping later - change 01:00 to 00:00 and 13:00 to 12:00
  dat_temp <- dat_temp %>%
    mutate(jday = yday(t),
           tod  = format(t, "%H:%M:%S")) %>%
    mutate(tod = str_replace_all(tod, "13:00:00", "12:00:00"),
           tod = str_replace_all(tod, "01:00:00", "00:00:00"))
  
  # Group by box, layer, julian day, and summarize (mean and sd as separate data frames)
  if(mean){
    aggregated_df <- dat_temp %>%
      group_by(b, z, jday, tod) %>%
      summarise(!!sym(variable) := mean(!!sym(variable))) 
  } else {
    aggregated_df <- dat_temp %>%
      group_by(b, z, jday, tod) %>%
      summarise(!!sym(variable) := sd(!!sym(variable))) 
  }
  
  # add time step index and clean up
  # In addition: for the purpose of applying the correct delta to the projections, 
  # the delta for Jan 1 00:00 (jday = 1, tod = 00:00) must be the last time step
  aggregated_df <- aggregated_df %>%
    group_by(jday, tod) %>%
    mutate(ts = cur_group_id()) %>%
    ungroup() %>%
    mutate(ts = ifelse(ts == 1, (max(ts)+1), ts)) %>% # put the first time step last
    mutate(ts = ts -1) %>%
    arrange(ts, b, z) %>%
    dplyr::select(ts, b, z, temperature)
  
  # now turn the data frame into an array with dimensions z,b,t
  # the correct dimensioning is important because we will use element-wise operations for the delta correction of the projection files
  # first restore missing levels
  all_cells <- expand.grid(ts = 1:max(aggregated_df$ts), 
                           b = 1:max(aggregated_df$b), 
                           z = 1:max(aggregated_df$z)) # 109 boxes, 7 layers
  
  aggregated_df_complete <- all_cells %>%
    left_join(aggregated_df) 
  
  # for leap years, eliminate ts for Feb 29
  # considering that ts starts on Jan 1 12:00, we have:
  # (31*2)-1 = 61 ts for jan
  # 28*2 = 56 ts for Feb
  # so, the first ts to drop is Feb 29 00:00, which is:
  # 61+56+1=(118,119)
  
  if(!leap){
    aggregated_df_complete <- aggregated_df_complete %>%
      filter(ts != 118, ts != 119)
    }
  
  # Reshape the data frame to a wider format for each 'ts'
  ts_list <- aggregated_df_complete %>%
    group_by(ts) %>%
    pivot_wider(names_from = b, values_from = temperature) %>%
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
                              hindcast = T,
                              mean = T,
                              leap = F)

sd_hind <- make_delta_array(variable = this_variable,
                              hindcast = T,
                              mean = F,
                              leap = F)

mean_hist <- make_delta_array(variable = this_variable,
                              hindcast = F,
                              mean = T,
                              leap = F)

sd_hist <- make_delta_array(variable = this_variable,
                              hindcast = F,
                              mean = F,
                              leap = F)

# for leap years
mean_hind_leap <- make_delta_array(variable = this_variable,
                              hindcast = T,
                              mean = T,
                              leap = T)

sd_hind_leap <- make_delta_array(variable = this_variable,
                            hindcast = T,
                            mean = F,
                            leap = T)

mean_hist_leap <- make_delta_array(variable = this_variable,
                              hindcast = F,
                              mean = T,
                              leap = T)

sd_hist_leap <- make_delta_array(variable = this_variable,
                            hindcast = F,
                            mean = F,
                            leap = T)

# xx
check <- mean_hind - mean_hist
check <- check[,,1:250]
check_df <- data.frame("ts" = 1:(dim(check)[3]), "mean_delta" = NA)
for(i in check_df$ts){
  check_df[i,2] <- mean(check[,,i], na.rm = T)
}

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

proj_files <- list.files("data/projection/", full.names = T)


for(i in 1:length(proj_files)){
  
  # point at pre-correction file
  original_proj_file <- proj_files[i]
  # just file name to re-write later
  fn <- gsub("data/projection/","",original_proj_file)
  
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
  corrected_proj_file <- paste0("data/projection_corrected/",fn)
  
  # copy the original necdf file that we are correcting
  file.copy(original_proj_file, corrected_proj_file)
  
  # open the new file
  corrected_proj_nc <- nc_open(corrected_proj_file, write = T)
  
  # write the delta-corrected variable to it
  ncvar_put(corrected_proj_nc, varid = this_variable, vals = corrected_proj_values)
  
  # close both nc files
  nc_close(original_proj_nc)
  nc_close(corrected_proj_nc)
  
}